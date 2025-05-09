import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Load Simulation CSV
# -----------------------------
def load_simulation_data(filename):
    """
    Load the simulation results CSV into a DataFrame.
    """
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print(f"Error: File not found - {filename}")
        return None


# -----------------------------
# Plot 1: SIR Line Plot per Rank
# -----------------------------
def plot_sir_per_rank(df):
    """
    Plot S, I, R lines separately for each MPI rank.
    """
    plt.figure(figsize=(10, 6))
    for rank in df["Rank"].unique():
        sub = df[df["Rank"] == rank]
        plt.plot(sub["Time"], sub["S_avg"], color="blue", alpha=0.5)
        plt.plot(sub["Time"], sub["I_avg"], color="orange", alpha=0.5)
        plt.plot(sub["Time"], sub["R_avg"], color="green", alpha=0.5)

    plt.xlabel("Time")
    plt.ylabel("Proportion of Population")
    plt.title("SIR Simulation Per Rank")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(["Susceptible (S)", "Infected (I)", "Recovered (R)"])
    plt.tight_layout()
    plt.savefig("./plots/sir_per_rank_lines.png")
    print("Saved: sir_per_rank_lines.png")


# -----------------------------
# Plot 2: Stacked Area Plot per Rank
# -----------------------------
def plot_stacked_area_per_rank(df):
    """
    Generate one stacked area chart showing rank-averaged composition.
    """
    plt.figure(figsize=(10, 6))
    for rank in df["Rank"].unique():
        sub = df[df["Rank"] == rank]
        plt.stackplot(sub["Time"], sub["S_avg"], sub["I_avg"], sub["R_avg"],
                      alpha=0.3)

    plt.xlabel("Time")
    plt.ylabel("Proportion of Population")
    plt.title("Stacked SIR Area (Per Rank Overlay)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig("./plots/sir_stacked_area_per_rank.png")
    print("Saved: sir_stacked_area_per_rank.png")


# -----------------------------
# Plot 3: Global Line Plot (Averaged Across Ranks)
# -----------------------------
def plot_sir_global_line(df):
    """
    Plot a single averaged S, I, R line over time.
    """
    df_avg = df.groupby("Time")[["S_avg", "I_avg", "R_avg"]].mean().reset_index()

    plt.figure(figsize=(10, 6))
    plt.plot(df_avg["Time"], df_avg["S_avg"], label="Susceptible (S)", linewidth=2)
    plt.plot(df_avg["Time"], df_avg["I_avg"], label="Infected (I)", linewidth=2)
    plt.plot(df_avg["Time"], df_avg["R_avg"], label="Recovered (R)", linewidth=2)

    plt.xlabel("Time")
    plt.ylabel("Proportion of Population")
    plt.title("Global Average SIR Line Plot")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("./plots/sir_global_line_plot.png")
    print("Saved: sir_global_line_plot.png")


# -----------------------------
# Plot 4: Global Stacked Area (Averaged)
# -----------------------------
def plot_global_stacked_area(df):
    """
    Stacked area showing average S/I/R across all ranks.
    """
    df_avg = df.groupby("Time")[["S_avg", "I_avg", "R_avg"]].mean().reset_index()

    plt.figure(figsize=(10, 6))
    plt.stackplot(df_avg["Time"], df_avg["S_avg"], df_avg["I_avg"], df_avg["R_avg"],
                  labels=["Susceptible", "Infected", "Recovered"], alpha=0.8)

    plt.xlabel("Time")
    plt.ylabel("Proportion of Population")
    plt.title("Global Average SIR Stacked Area")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("./plots/sir_global_stacked_area.png")
    print("Saved: sir_global_stacked_area.png")


# -----------------------------
# Run all plots
# -----------------------------
if __name__ == "__main__":
    filename = "./data/simulation_results.csv"
    df = load_simulation_data(filename)

    if df is not None:
        plot_sir_per_rank(df)
        plot_stacked_area_per_rank(df)
        plot_sir_global_line(df)
        plot_global_stacked_area(df)
