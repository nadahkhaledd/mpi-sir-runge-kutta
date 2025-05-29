import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def analyze_test_results(base_path):
    # Create results directory
    results_dir = os.path.join(base_path, 'test_results')
    os.makedirs(results_dir, exist_ok=True)

    # Load all test results
    results = []
    for file in os.listdir(base_path):
        if file.endswith('_results.csv'):
            df = pd.read_csv(os.path.join(base_path, file))
            test_name = file.replace('_results.csv', '')
            df['Test'] = test_name
            results.append(df)

    if not results:
        print("No test results found!")
        return

    all_results = pd.concat(results)

    # Generate plots
    # 1. Performance comparison
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=all_results, x='Processes', y='ExecTime', hue='Test')
    plt.title('Execution Time by Number of Processes')
    plt.savefig(os.path.join(results_dir, 'performance_comparison.png'))
    plt.close()

    # 2. SIR curves for each test
    for test in all_results['Test'].unique():
        test_data = all_results[all_results['Test'] == test]
        plt.figure(figsize=(10, 6))
        plt.plot(test_data['Time'], test_data['S'], label='Susceptible')
        plt.plot(test_data['Time'], test_data['I'], label='Infected')
        plt.plot(test_data['Time'], test_data['R'], label='Recovered')
        plt.title(f'SIR Curves - {test}')
        plt.xlabel('Time')
        plt.ylabel('Population Fraction')
        plt.legend()
        plt.savefig(os.path.join(results_dir, f'sir_curves_{test}.png'))
        plt.close()

    # Generate summary statistics
    summary = all_results.groupby(['Test', 'Processes'])['ExecTime'].agg(['mean', 'std', 'min', 'max'])
    summary.to_csv(os.path.join(results_dir, 'performance_summary.csv'))

if __name__ == "__main__":
    analyze_test_results('./data')
