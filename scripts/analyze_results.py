import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
import numpy as np

def analyze_results(test_dir="./data/test_results", 
                   main_results="./data/output/simulation_results.csv",
                   output_dir="./data/analysis"):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print("\nStarting analysis...")
    
    # Structure analysis outputs
    subdirs = ['comparisons', 'errors', 'summaries']
    for subdir in subdirs:
        Path(os.path.join(output_dir, subdir)).mkdir(parents=True, exist_ok=True)

    # Load and validate reference data
    print("\nLoading and validating data...")
    reference_df = pd.read_csv(main_results)
    print(f"Reference data shape: {reference_df.shape}")
    print(f"Reference data columns: {reference_df.columns.tolist()}")

    # Organize tests by type
    test_categories = {}
    
    for file in os.listdir(test_dir):
        if not file.endswith('_results.csv'):
            continue
            
        test_name = file.split('_p')[0]
        test_df = pd.read_csv(os.path.join(test_dir, file))
        
        # Determine test category from filename
        if 'sensitivity' in test_name:
            category = 'sensitivity'
        else:
            category = 'temporal'
            
        if category not in test_categories:
            test_categories[category] = {}
            
        test_categories[category][test_name] = test_df
        print(f"Loaded {category} test: {test_name}")

    # Create comparison plots by category
    for category, tests in test_categories.items():
        print(f"\nProcessing {category} tests...")
        
        # Create category directory
        category_dir = os.path.join(output_dir, category)
        Path(category_dir).mkdir(parents=True, exist_ok=True)
        
        # Plot all tests in this category together
        plt.figure(figsize=(12, 8))
        
        # Reference data
        plt.plot(reference_df['Time'], reference_df['S_avg'], 'k-', 
                label='Reference', linewidth=2)
        
        # Test data with different colors for each test
        colors = plt.cm.rainbow(np.linspace(0, 1, len(tests)))
        for (test_name, test_df), color in zip(tests.items(), colors):
            plt.plot(test_df['Time'], test_df['S_avg'], '--', 
                    color=color, label=test_name, linewidth=1.5)
        
        plt.title(f'{category.title()} Tests Comparison')
        plt.xlabel('Time')
        plt.ylabel('Population Fraction')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        # Save category plot
        save_path = os.path.join(category_dir, f'{category}_comparison.png')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {save_path}")

    # Basic validation
    for category, tests in test_categories.items():
        for test_name, test_df in tests.items():
            validate_data(test_df, test_name)
            
            # Calculate and save error metrics
            calculate_errors(reference_df, test_df, test_name, output_dir)

    print("\nAnalysis complete! Check the following directories:")
    for subdir in subdirs:
        print(f"- {os.path.join(output_dir, subdir)}")

def validate_data(df, name):
    """Validates SIR constraints and data quality."""
    print(f"\nValidating {name}:")
    print(f"- Shape: {df.shape}")
    print(f"- Time range: [{df['Time'].min():.2f}, {df['Time'].max():.2f}]")
    
    # Check SIR constraints
    sir_sum = df['S_avg'] + df['I_avg'] + df['R_avg']
    if not np.allclose(sir_sum, 1.0, rtol=1e-5):
        print(f"Warning: S+I+R != 1 (mean={sir_sum.mean():.6f})")

def plot_comparison(ref_df, test_df, test_name, output_dir):
    """Creates detailed comparison plots."""
    plt.figure(figsize=(12, 8))
    
    # Reference data (solid lines)
    plt.plot(ref_df['Time'], ref_df['S_avg'], 'b-', 
            label='Training - Susceptible', linewidth=2)
    plt.plot(ref_df['Time'], ref_df['I_avg'], 'r-', 
            label='Training - Infected', linewidth=2)
    plt.plot(ref_df['Time'], ref_df['R_avg'], 'g-', 
            label='Training - Recovered', linewidth=2)

    # Test data (dashed lines)
    plt.plot(test_df['Time'], test_df['S_avg'], 'b--', 
            label=f'Test - Susceptible', linewidth=2)
    plt.plot(test_df['Time'], test_df['I_avg'], 'r--', 
            label=f'Test - Infected', linewidth=2)
    plt.plot(test_df['Time'], test_df['R_avg'], 'g--', 
            label=f'Test - Recovered', linewidth=2)

    plt.title(f'SIR Model: Training vs {test_name}')
    plt.xlabel('Time')
    plt.ylabel('Population Fraction')
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Save individual test comparison
    save_path = os.path.join(output_dir, 'comparisons', f'{test_name}_comparison.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved comparison plot: {save_path}")

def calculate_errors(ref_df, test_df, test_name, output_dir):
    """Calculates detailed error metrics."""
    errors = {
        'Susceptible': np.abs(test_df['S_avg'] - ref_df['S_avg']),
        'Infected': np.abs(test_df['I_avg'] - ref_df['I_avg']),
        'Recovered': np.abs(test_df['R_avg'] - ref_df['R_avg'])
    }
    
    error_summary = pd.DataFrame(errors)
    error_summary['Time'] = test_df['Time']
    
    # Save error metrics to CSV
    save_path = os.path.join(output_dir, 'errors', f'{test_name}_error_metrics.csv')
    error_summary.to_csv(save_path, index=False)
    print(f"Saved error metrics: {save_path}")

    # Create error plot
    plt.figure(figsize=(12, 6))
    
    for label, error in errors.items():
        plt.plot(test_df['Time'], error, label=label, linewidth=2)
    
    plt.title(f'Absolute Error: {test_name} vs Training')
    plt.xlabel('Time')
    plt.ylabel('|Test - Training|')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save error plot
    save_path = os.path.join(output_dir, 'errors', f'{test_name}_error.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved error plot: {save_path}")

if __name__ == "__main__":
    analyze_results()
