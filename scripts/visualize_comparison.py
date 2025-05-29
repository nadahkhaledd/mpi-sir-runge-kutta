import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from matplotlib.gridspec import GridSpec

def create_comparison_plots(test_results_dir, output_dir):
    """Create comprehensive comparison visualizations"""
    
    # Load all test results
    results = []
    for file in os.listdir(test_results_dir):
        if file.endswith('_results.csv'):
            df = pd.read_csv(os.path.join(test_results_dir, file))
            process_count = int(file.split('_p')[1].split('_')[0])
            test_name = file.split('_p')[0]
            df['Processes'] = process_count
            df['Test'] = test_name
            results.append(df)
    
    all_results = pd.concat(results)
    
    # 1. Process scalability plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=all_results, x='Processes', y='Runtime')
    plt.title('Runtime vs Number of Processes')
    plt.savefig(os.path.join(output_dir, 'scalability.png'))
    plt.close()
    
    # 2. Parameter sensitivity heatmap
    sensitivity_data = all_results.pivot_table(
        values='max_infection_rate', 
        index='beta',
        columns='gamma'
    )
    plt.figure(figsize=(10, 8))
    sns.heatmap(sensitivity_data, annot=True, cmap='YlOrRd')
    plt.title('Parameter Sensitivity Analysis')
    plt.savefig(os.path.join(output_dir, 'parameter_sensitivity.png'))
    plt.close()
    
    # 3. Multi-test comparison
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 2, figure=fig)
    
    # SIR curves
    ax1 = fig.add_subplot(gs[0, :])
    for test in all_results['Test'].unique():
        test_data = all_results[all_results['Test'] == test]
        ax1.plot(test_data['Time'], test_data['S'], label=f'{test} S')
        ax1.plot(test_data['Time'], test_data['I'], label=f'{test} I')
        ax1.plot(test_data['Time'], test_data['R'], label=f'{test} R')
    ax1.set_title('SIR Curves Comparison')
    ax1.legend()
    
    # Runtime distribution
    ax2 = fig.add_subplot(gs[1, 0])
    sns.violinplot(data=all_results, x='Test', y='Runtime', ax=ax2)
    ax2.set_title('Runtime Distribution')
    
    # Error metrics
    ax3 = fig.add_subplot(gs[1, 1])
    error_metrics = calculate_error_metrics(all_results)
    sns.barplot(data=error_metrics, ax=ax3)
    ax3.set_title('Error Metrics by Test')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'comprehensive_comparison.png'))
    plt.close()

def calculate_error_metrics(results_df, baseline_df=None):
    """Calculate various error metrics for comparison"""
    metrics = []
    
    if baseline_df is not None:
        # Compare with baseline data
        for test in results_df['Test'].unique():
            test_data = results_df[results_df['Test'] == test]
            mse = np.mean((test_data['I'] - baseline_df['I'])**2)
            mae = np.mean(np.abs(test_data['I'] - baseline_df['I']))
            metrics.append({
                'Test': test,
                'MSE': mse,
                'MAE': mae
            })
    
    return pd.DataFrame(metrics)

if __name__ == "__main__":
    create_comparison_plots(
        './data/test_results',
        './data/analysis'
    )
