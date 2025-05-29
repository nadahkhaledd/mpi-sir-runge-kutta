import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns
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
# Plot 5: Heatmap – Infection Intensity Over Time (Per Rank)
# -----------------------------
def plot_infection_heatmap_per_rank(df, vmin=0.0, vmax=0.25):
    """
    Creates a heatmap of the average infected proportion (I_avg)
    across time for each MPI rank.

    - The x-axis represents time steps.
    - The y-axis corresponds to MPI ranks (processes).
    - Cell values (color) show the level of infection.

    This visualization is helpful for:
    - Spotting temporal infection peaks across ranks
    - Identifying which regions were hit harder or earlier
    - Debugging load balance or propagation dynamics

    Parameters:
    - df: Pandas DataFrame containing simulation results.
    - vmin: Minimum value for the heatmap color scale (default: 0.0).
    - vmax: Maximum value for the heatmap color scale (default: 0.25).
      Adjust these to control sensitivity/contrast in the plot.

    NOTE:
    - Assumes the `Rank`, `Time`, and `I_avg` columns are present.
    - Grouped by Rank → Time for matrix-style heatmap.
    """

    # Pivot to get a matrix: Rows = Ranks, Columns = Time, Values = I_avg
    pivoted = df.pivot(index="Rank", columns="Time", values="I_avg")

    plt.figure(figsize=(12, 6))
    sns.heatmap(
        pivoted,
        cmap="YlOrRd",
        vmin=vmin,
        vmax=vmax,
        cbar_kws={"label": "Infected Proportion"}
    )
    plt.title("Heatmap of Infection Intensity Over Time (Per Rank)")
    plt.xlabel("Time")
    plt.ylabel("MPI Rank")
    plt.tight_layout()
    plt.savefig("./plots/infection_heatmap_per_rank.png")
    print("Saved: infection_heatmap_per_rank.png")




# -----------------------------
# Plot 6: Animated Infection Spread (Per Rank Over Time)
# -----------------------------
def animate_infection(df):
    """
    Generates an animated bar chart showing the spread of infection (I_avg)
    across MPI ranks over time.

    - Each frame of the animation corresponds to a timestep.
    - Bar height represents the infection proportion for that rank.
    - Useful for visually tracking temporal progression of the epidemic
      and how it impacts different ranks differently.

    Parameters:
    - df: Pandas DataFrame with columns 'Rank', 'Time', 'I_avg'.

    Output:
    - Saves an animated .gif of the infection evolution to ./plots/infection_animation.gif.
    """

    ranks = df["Rank"].unique()
    times = sorted(df["Time"].unique())
    fig, ax = plt.subplots(figsize=(10, 6))
    bar = ax.bar(ranks, [0] * len(ranks), color="orange")
    max_i = df["I_avg"].max()
    ax.set_ylim(0, max_i * 1.2)  # Add headroom

    def update(t):
        time_data = df[df["Time"] == t]
        i_vals = [time_data[time_data["Rank"] == r]["I_avg"].values[0]
                  if r in time_data["Rank"].values else 0 for r in ranks]
        for rect, h in zip(bar, i_vals):
            rect.set_height(h)
        ax.set_title(f"Infection Spread at t = {t:.1f}")
        return bar

    ani = animation.FuncAnimation(fig, update, frames=times, blit=False, repeat=False)
    ani.save("./plots/infection_animation.gif", writer="pillow", fps=5)
    print("Saved: infection_animation.gif")

# -----------------------------
# Plot 7: Peak Infection Time per MPI Rank
# -----------------------------
def plot_peak_infection_times(df):
    """
    Plots a bar chart showing the time at which each MPI rank experienced
    its peak infection level.

    - X-axis: MPI Rank
    - Y-axis: Time of highest infection (I_avg)
    - Highlights temporal heterogeneity of epidemic spread

    Parameters:
    - df: Pandas DataFrame containing columns 'Rank', 'Time', and 'I_avg'.

    Output:
    - Saves the figure to ./plots/peak_infection_time_per_rank.png
    """

    peak_times = df.loc[df.groupby("Rank")["I_avg"].idxmax()][["Rank", "Time", "I_avg"]]

    plt.figure(figsize=(10, 5))
    plt.bar(peak_times["Rank"], peak_times["Time"], color="purple", alpha=0.7)
    plt.xlabel("MPI Rank")
    plt.ylabel("Time of Peak Infection")
    plt.title("Peak Infection Timing per Rank")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("./plots/peak_infection_time_per_rank.png")
    print("Saved: peak_infection_time_per_rank.png")

# -----------------------------
# Plot 8: Timing Comparison Across Simulation Phases
# -----------------------------
def plot_timing_comparison_phases(timing_df, output_dir="./plots"):
    """
    Generates a bar plot comparing Min, Max, and Avg times
    across simulation phases.

    Parameters:
    - timing_df: DataFrame with 'PhaseName', 'Statistic', and 'Value' columns.
    - output_dir: Directory where the plot will be saved.
    """
    df_general = timing_df[timing_df['Statistic'].isin(['Min', 'Max', 'Avg'])]

    plt.figure(figsize=(12, 6))
    sns.barplot(data=df_general, x="PhaseName", y="Value", hue="Statistic")
    plt.title("Timing Comparison Across Simulation Phases")
    plt.ylabel("Time (s)")
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y')
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "timing_comparison_phases.png")
    plt.savefig(output_path)
    print(f"Saved: {output_path}")
    plt.close()


# -----------------------------
# Plot 9: Rank_0_Time by Phase
# -----------------------------
def plot_rank0_time_phases(timing_df, output_dir="./plots"):
    """
    Bar chart showing Rank_0_Time across simulation phases.
    """
    df_rank0 = timing_df[timing_df['Statistic'] == 'Rank_0_Time']

    plt.figure(figsize=(10, 5))
    sns.barplot(data=df_rank0, x="PhaseName", y="Value", color="skyblue")
    plt.title("Rank_0_Time Across MPI and Computation Phases")
    plt.ylabel("Time (s)")
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y')
    plt.tight_layout()
    output_path = os.path.join(output_dir, "rank0_time_phases.png")
    plt.savefig(output_path)
    print(f"Saved: {output_path}")
    plt.close()


# -----------------------------
# Run all plots
# -----------------------------
if __name__ == "__main__":

    # Define file paths
    sim_data_file = "./data/output/simulation_results.csv"
    timing_log_file = "./data/output/timing_log.csv"
    output_dir = "./plots"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # --- Load and plot simulation data ---
    df = load_simulation_data(sim_data_file)
    if df is not None:
        print("Generating SIR and infection plots...")
        plot_sir_per_rank(df)
        plot_stacked_area_per_rank(df)
        plot_sir_global_line(df)
        plot_global_stacked_area(df)
        plot_infection_heatmap_per_rank(df)
        animate_infection(df)
        plot_peak_infection_times(df)
        print("Simulation plots generated successfully.")
    else:
        print("Simulation data not found. Skipping simulation plots.")

    # --- Load and plot timing data ---
    if os.path.exists(timing_log_file):
        timing_df = pd.read_csv(timing_log_file)
        print("Generating timing analysis plots...")
        plot_timing_comparison_phases(timing_df, output_dir)
        plot_rank0_time_phases(timing_df, output_dir)
        print("Timing plots generated successfully.")
    else:
        print("Timing log not found. Skipping timing plots.")
