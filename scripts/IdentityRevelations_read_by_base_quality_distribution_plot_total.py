#!/usr/bin/env python3

"""
Generate comparison violin plots for quality score distributions across 4 sequencing platforms.
Optimized for large datasets (20M+ reads) using chunked reading with multiprocessing.

Creates 12 subplots (4 error categories × 3 quality bins):
- Rows: 1_error, 2_errors, 3_errors, gt3_errors
- Columns: Q0-19 (Low), Q20-29 (Medium), Q30+ (High)
- Each subplot shows violin plots for 4 platforms:
  1. Element AVITI
  2. GeneMind SURFSeq 5000
  3. Illumina NovaSeq X
  4. MGI T1+
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
import multiprocessing as mp
from functools import partial
from matplotlib.ticker import LogLocator, FuncFormatter


# Platform names in order
PLATFORMS = ["Element AVITI", "GeneMind SURFSeq", "Illumina NovaSeq X", "MGI T1+"]
COLORS = ["#FFE7B3", "#6BFFD6", "#55A9FF", "#F0E0F0"]


def load_tsv_chunked(tsv_path, chunksize=500000):
    """
    Load TSV file in chunks (memory-efficient for large files).
    
    Args:
        tsv_path: path to TSV file (.tsv or .tsv.gz)
        chunksize: number of rows to read per chunk
    
    Yields:
        pd.DataFrame chunks
    """
    if tsv_path.endswith('.gz'):
        return pd.read_csv(tsv_path, sep='\t', compression='gzip', chunksize=chunksize)
    else:
        return pd.read_csv(tsv_path, sep='\t', chunksize=chunksize)


def process_single_file(args):
    """
    Process a single TSV file in chunks.
    
    Args:
        args: tuple (tsv_path, platform_name, chunksize)
    
    Returns:
        dict with aggregated data for each error class
    """
    tsv_path, platform_name, chunksize = args
    
    error_order = ["1_error", "2_errors", "3_errors", "gt3_errors"]
    
    # Initialize data structure for each error class
    data_by_error = {err: {'q0_19': [], 'q20_29': [], 'q30plus': []} for err in error_order}
    
    print(f"[{platform_name}] Reading {tsv_path}...")
    
    chunk_count = 0
    total_rows = 0
    filtered_rows = 0
    
    for chunk in load_tsv_chunked(tsv_path, chunksize=chunksize):
        chunk_count += 1
        total_rows += len(chunk)
        
        # Exclude reads with 0 errors
        chunk = chunk[chunk['error_class'] != '0_errors']
        
        # Filter to valid error classes
        chunk = chunk[chunk['error_class'].isin(error_order)]
        filtered_rows += len(chunk)
        
        # Aggregate data by error class
        for error_class in error_order:
            subset = chunk[chunk['error_class'] == error_class]
            if len(subset) > 0:
                data_by_error[error_class]['q0_19'].extend(subset['q0_19'].values.tolist())
                data_by_error[error_class]['q20_29'].extend(subset['q20_29'].values.tolist())
                data_by_error[error_class]['q30plus'].extend(subset['q30plus'].values.tolist())
    
    print(f"[{platform_name}] ✔ Finished: {total_rows:,} rows total, {filtered_rows:,} rows after filtering")
    
    return {
        'platform': platform_name,
        'data': data_by_error
    }


def process_all_files_parallel(tsv_paths, threads=4, chunksize=500000):
    """
    Process all TSV files in parallel.
    
    Args:
        tsv_paths: list of 4 TSV file paths [Element, GeneMind, Illumina, MGI]
        threads: number of parallel workers
        chunksize: rows per chunk
    
    Returns:
        dict with data for each platform
    """
    print(f"\nProcessing {len(tsv_paths)} platforms in parallel ({threads} threads)...\n")
    
    args_list = [
        (tsv_paths[i], PLATFORMS[i], chunksize)
        for i in range(len(tsv_paths))
    ]
    
    with mp.Pool(threads) as pool:
        results = pool.map(process_single_file, args_list)
    
    # Reorganize results: {platform: {error_class: {q_bin: [values]}}}
    all_data = {}
    for result in results:
        all_data[result['platform']] = result['data']
    
    print(f"\n✔ All files processed\n")
    return all_data


def create_comparison_plots(all_data, output_prefix):
    """
    Create 12-subplot comparison figure (4 error categories × 3 quality bins).
    
    Args:
        all_data: dict with {platform: {error_class: {q_bin: [values]}}}
        output_prefix: prefix for output file names
    """
    
    error_order = ["1_error", "2_errors", "3_errors", "gt3_errors"]
    q_bins = ["q0_19", "q20_29", "q30plus"]
    q_titles = ["Low Qscore (Q0-19)", "Medium Qscore (Q20-29)", "High Qscore (Q30+)"]
    
    # Set style
    sns.set_style("whitegrid")
    
    # Create figure with 4 rows (error categories) × 3 columns (quality bins)
    fig, axes = plt.subplots(4, 3, figsize=(20, 25))
    fig.suptitle('Quality Score Distribution Comparison Across Sequencing Platforms\n(Reads with 1+ Errors)', 
                 fontsize=18, fontweight='bold', y=0.995)
    
    # Iterate through error categories (rows)
    for row_idx, error_class in enumerate(error_order):
        
        # Iterate through quality bins (columns)
        for col_idx, (q_bin, q_title) in enumerate(zip(q_bins, q_titles)):
            
            ax = axes[row_idx, col_idx]
            
            # Prepare data for this subplot
            plot_data_list = []
            
            for platform in PLATFORMS:
                if platform in all_data and len(all_data[platform][error_class][q_bin]) > 0:
                    values = all_data[platform][error_class][q_bin]
                    platform_data = pd.DataFrame({
                        'platform': [platform] * len(values),
                        'value': values
                    })
                    plot_data_list.append(platform_data)
            
            if plot_data_list:
                df_plot = pd.concat(plot_data_list, ignore_index=True)
                
                # # Create violin plot
                # sns.violinplot(
                #     data=df_plot,
                #     x='platform',
                #     y='value',
                #     ax=ax,
                #     palette=COLORS,
                #     inner='box'
                # )

                sns.boxplot(
                    data=df_plot,
                    x='platform',
                    y='value',
                    ax=ax,
                    palette=COLORS,
                    showfliers=False,
                    linewidth=2
                )
                if q_bin=="q30plus":
                    ax.set_ylim(120, 150)
                elif q_bin=="q20_29":
                    ax.set_ylim(0, 15)
                else:
                    ax.set_ylim(0, 15)


                # Set labels and title
                if row_idx == 0:  # Top row
                    ax.set_title(q_title, fontsize=13, fontweight='bold', pad=10)
                
                if col_idx == 0:  # Left column
                    ax.set_ylabel(f'{error_class}\nNumber of Bases', fontsize=11, fontweight='bold')
                else:
                    ax.set_ylabel('')
                
                # Set x-axis labels for all rows
                ax.set_xticklabels(PLATFORMS, rotation=45, ha='right')

                # Only add xlabel to bottom row
                if row_idx == 3:
                    ax.set_xlabel('Platform', fontsize=11, fontweight='bold')
                else:
                    ax.set_xlabel('')
                
                ax.grid(axis='y', alpha=0.3)
                ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    # Save figure
    output_file = f"{output_prefix}_platform_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✔ Comparison plot saved: {output_file}")
    
    output_pdf = f"{output_prefix}_platform_comparison.pdf"
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    print(f"✔ Comparison plot (PDF) saved: {output_pdf}")
    
    plt.close()


def generate_summary_stats(all_data):
    """Generate summary statistics for each platform and error category."""
    
    error_order = ["1_error", "2_errors", "3_errors", "gt3_errors"]
    
    print("\n" + "="*100)
    print("SUMMARY STATISTICS BY PLATFORM AND ERROR CATEGORY")
    print("="*100)
    
    for platform in PLATFORMS:
        print(f"\n{'='*100}")
        print(f"PLATFORM: {platform}")
        print(f"{'='*100}")
        
        if platform not in all_data:
            print(f"  No data found for {platform}")
            continue
        
        for error_class in error_order:
            q0_19 = np.array(all_data[platform][error_class]['q0_19'])
            q20_29 = np.array(all_data[platform][error_class]['q20_29'])
            q30plus = np.array(all_data[platform][error_class]['q30plus'])
            
            if len(q0_19) == 0:
                print(f"\n  {error_class}: No reads found")
                continue
            
            print(f"\n  {error_class}:")
            print(f"    Number of reads: {len(q0_19):,}")
            
            print(f"\n    Q0-19 bases (Low Quality):")
            print(f"      Mean:   {np.mean(q0_19):8.2f}")
            print(f"      Median: {np.median(q0_19):8.2f}")
            print(f"      Std:    {np.std(q0_19):8.2f}")
            print(f"      Min:    {np.min(q0_19):8.0f}")
            print(f"      Max:    {np.max(q0_19):8.0f}")
            
            print(f"\n    Q20-29 bases (Medium Quality):")
            print(f"      Mean:   {np.mean(q20_29):8.2f}")
            print(f"      Median: {np.median(q20_29):8.2f}")
            print(f"      Std:    {np.std(q20_29):8.2f}")
            print(f"      Min:    {np.min(q20_29):8.0f}")
            print(f"      Max:    {np.max(q20_29):8.0f}")
            
            print(f"\n    Q30+ bases (High Quality):")
            print(f"      Mean:   {np.mean(q30plus):8.2f}")
            print(f"      Median: {np.median(q30plus):8.2f}")
            print(f"      Std:    {np.std(q30plus):8.2f}")
            print(f"      Min:    {np.min(q30plus):8.0f}")
            print(f"      Max:    {np.max(q30plus):8.0f}")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("USAGE: python plot_quality_platforms_comparison.py output_prefix threads element.tsv genemind.tsv illumina.tsv mgi.tsv")
        print("\nArguments:")
        print("  output_prefix: prefix for output files")
        print("  threads: number of parallel workers (e.g., 4)")
        print("  element.tsv: TSV file for Element AVITI")
        print("  genemind.tsv: TSV file for GeneMind SURFSeq")
        print("  illumina.tsv: TSV file for Illumina NovaSeq X")
        print("  mgi.tsv: TSV file for MGI T1+")
        print("\nExample:")
        print("  python plot_quality_platforms_comparison.py results 4 element.tsv.gz genemind.tsv.gz illumina.tsv.gz mgi.tsv.gz")
        sys.exit(1)
    
    output_prefix = sys.argv[1]
    threads = int(sys.argv[2])
    tsv_paths = sys.argv[3:7]
    
    print(f"\n{'='*100}")
    print("SEQUENCING PLATFORM COMPARISON - QUALITY DISTRIBUTION ANALYSIS")
    print(f"{'='*100}\n")
    print(f"Output prefix: {output_prefix}")
    print(f"Parallel threads: {threads}")
    print(f"\nPlatforms:")
    for i, (platform, path) in enumerate(zip(PLATFORMS, tsv_paths)):
        print(f"  {i+1}. {platform:25s} → {path}")
    
    # Process all files in parallel
    all_data = process_all_files_parallel(tsv_paths, threads=threads)
    
    # Generate comparison plots
    create_comparison_plots(all_data, output_prefix)
    
    # Print summary statistics
    # generate_summary_stats(all_data)
    
    print("\n" + "="*100)
    print("✔ Analysis complete!")
    print("="*100 + "\n")
