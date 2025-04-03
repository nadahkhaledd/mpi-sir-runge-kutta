#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <cmath>

#include <mpi.h>
#include <opencv2/core.hpp>
#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>

struct StateData {
    std::string name;
    double population;
    std::string date;
    double lat, lon;
    int confirmed, deaths, recovered, active;
    int blockId; // New field
};

// Function to read the CSV file
std::vector<StateData> readCSV(const std::string& filename) {
    std::vector<StateData> data;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return data;
    }

    std::getline(file, line); // Skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        StateData entry;
        std::string temp;

        std::getline(ss, entry.name, ',');
        std::getline(ss, temp, ','); entry.population = std::stod(temp);
        std::getline(ss, entry.date, ',');
        std::getline(ss, temp, ','); entry.lat = std::stod(temp);
        std::getline(ss, temp, ','); entry.lon = std::stod(temp);
        std::getline(ss, temp, ','); entry.confirmed = std::stoi(temp);
        std::getline(ss, temp, ','); entry.deaths = std::stoi(temp);
        std::getline(ss, temp, ','); entry.recovered = std::stoi(temp);
        std::getline(ss, temp, ','); entry.active = std::stoi(temp);
        
        data.push_back(entry);
    }
    return data;
}

// Function to assign block IDs using OpenCV k-means
void assignBlockIds(std::vector<StateData>& data, int numBlocks) {
    // Prepare the data for OpenCV k-means
    std::vector<cv::Point2f> locations;
    for (const auto& entry : data) {
        locations.emplace_back(entry.lat, entry.lon);
    }

    // Convert the data to OpenCV's cv::Mat (points matrix)
    cv::Mat points(locations.size(), 1, CV_32FC2);  // CV_32FC2 means 2D float matrix
    for (size_t i = 0; i < locations.size(); ++i) {
        points.at<cv::Vec2f>(i) = cv::Vec2f(locations[i].x, locations[i].y);
    }

    // Apply k-means clustering using OpenCV
    cv::Mat labels;
    cv::Mat centers;
    int attempts = 10; // Number of times the algorithm is executed with different initializations
    cv::kmeans(points, numBlocks, labels, cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 100, 0.2), attempts, cv::KMEANS_PP_CENTERS, centers);

    // Assign block IDs based on the cluster labels
    for (size_t i = 0; i < data.size(); ++i) {
        data[i].blockId = labels.at<int>(i);  // Set blockId from the k-means cluster label
    }
}

void saveCSV(const std::vector<StateData>& data, const std::string& filename) {
    std::ofstream file(filename);
    file << "Province_State,Population,Date,Lat,Long,Confirmed,Deaths,Recovered,Active,BlockId\n";
    
    for (const auto& entry : data) {
        file << entry.name << "," << entry.population << "," << entry.date << "," 
             << entry.lat << "," << entry.lon << "," << entry.confirmed << ","
             << entry.deaths << "," << entry.recovered << "," << entry.active << ","
             << entry.blockId << "\n";
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<StateData> dataset;
    if (rank == 0) {
        // Read the dataset on rank 0 and distribute to all other ranks
        dataset = readCSV("covid_data.csv");

        // Distribute the data among MPI processes (you can use MPI_Scatter here)
        int dataPerProcess = dataset.size() / size;
        std::vector<StateData> localData(dataPerProcess);

        // Scatter the dataset across processes
        MPI_Scatter(dataset.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                    localData.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                    0, MPI_COMM_WORLD);

        // Each process runs the k-means algorithm on its local data
        assignBlockIds(localData, dataPerProcess / 4);

        // Gather the results back to rank 0
        MPI_Gather(localData.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                   dataset.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                   0, MPI_COMM_WORLD);

        // Rank 0 saves the updated data
        saveCSV(dataset, "updated_covid_data.csv");
    } else {
        // Non-rank 0 processes receive their part of the data
        int dataPerProcess = dataset.size() / size;
        std::vector<StateData> localData(dataPerProcess);

        MPI_Scatter(nullptr, 0, MPI_BYTE, localData.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                    0, MPI_COMM_WORLD);

        // Perform clustering locally
        assignBlockIds(localData, dataPerProcess / 4);

        // Send the results back to rank 0
        MPI_Gather(localData.data(), dataPerProcess * sizeof(StateData), MPI_BYTE, 
                   nullptr, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
