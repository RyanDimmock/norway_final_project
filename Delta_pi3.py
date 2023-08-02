import pandas as pd
import matplotlib.pyplot as plt

# Load in data
SOR_df = pd.read_csv('50SNPs_final_allsites_SORvsNonspring.csv')
VES_df = pd.read_csv('50SNPs_final_allsites_VESvsNonspring.csv')
common_genes = pd.read_csv("common_spring_Gpos_with_scaf_startend_ordered.tsv", sep="\t")

# Calculate change in nucleotide diversity vs control
SOR_df['delta_pi'] = SOR_df['SL2_SOR_pi'] - SOR_df['north_nonspring_pi']
VES_df['delta_pi'] = VES_df['SL1_VES_pi'] - VES_df['north_nonspring_pi']

# Extract scaffold number
SOR_df['scaffold'] = SOR_df['scaffold'].str.extract('(\d+)', expand=False).astype(int)
VES_df['scaffold'] = VES_df['scaffold'].str.extract('(\d+)', expand=False).astype(int)

# Set window size
window_size = 1000000

# Determine scaffold ranges
max_position = SOR_df.groupby('scaffold')['midpoint'].max()
scaf_ranges = {scaf: (sum(max_position[:i]), sum(max_position[:i+1])) for i, scaf in enumerate(max_position.index)}

# Adjust the midpoint values
SOR_df['adjusted_position'] = SOR_df.apply(lambda row: scaf_ranges[row['scaffold']][0] + row['midpoint'], axis=1)
VES_df['adjusted_position'] = VES_df.apply(lambda row: scaf_ranges[row['scaffold']][0] + row['midpoint'], axis=1)
common_genes['adjusted_position'] = common_genes.apply(lambda row: scaf_ranges[row['scaffold']][0] + row['midpoint'], axis=1)

# Adjust the start_pos and end_pos values
common_genes['adjusted_start_pos'] = common_genes.apply(lambda row: scaf_ranges[row['scaffold']][0] + row['start_pos'], axis=1)
common_genes['adjusted_end_pos'] = common_genes.apply(lambda row: scaf_ranges[row['scaffold']][0] + row['end_pos'], axis=1)

# Sort data by adjusted midpoint
SOR_df_sorted = SOR_df.sort_values('adjusted_position')
VES_df_sorted = VES_df.sort_values('adjusted_position')
common_genes_sorted = common_genes.sort_values('adjusted_position')

# plot window function
def plot_window(window_start, window_end, gene_list, SOR_df_sorted, VES_df_sorted):
    data_window = SOR_df_sorted[(SOR_df_sorted['adjusted_position'] >= window_start) & (SOR_df_sorted['adjusted_position'] <= window_end)]
    new_data_window = VES_df_sorted[(VES_df_sorted['adjusted_position'] >= window_start) & (VES_df_sorted['adjusted_position'] <= window_end)]

    # Plot the data for this window
    plt.figure(figsize=(10, 5))
    plt.plot(data_window['adjusted_position'], data_window['delta_pi'], color='blue', linewidth=0.8, alpha=0.7, zorder=3)
    plt.plot(new_data_window['adjusted_position'], -new_data_window['delta_pi'], color='orange', linewidth=0.8, alpha=0.7, zorder=3)

    # Plot common gene line
    for i, gene_name in gene_list:
        gene = common_genes_sorted.loc[i]
        plt.axvspan(gene['adjusted_start_pos'], gene['adjusted_end_pos'], color='lightgreen', alpha=0.5, linewidth=0.5, zorder=2)

    # Add gene names
    gene_text = '\n'.join([f'{i+1}. {name}' for i, name in gene_list])
    plt.text(window_end, 0.023, gene_text, horizontalalignment='right', verticalalignment='top', fontsize=9, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5'))

    # Set titles
    plt.title(f'Scaffold {common_genes_sorted.loc[gene_list[0][0]]["scaffold"]}')  # using the scaffold as the title
    plt.xlabel('Position')
    plt.ylabel('Nucleotide Diversity (pi) vs non-spring')

    # set tciks and limits
    plt.xticks([window_start, window_end], [f'{int(window_start):,}', f'{int(window_end):,}'])
    plt.xlim(window_start, window_end)
    plt.ylim(-0.05, 0.05)
    plt.yticks([-0.02, 0, 0.02])
    plt.grid(axis='y', which='both', linewidth=0.5)

    # Save
    filename = f'scaffold_plot_{plot_number}.pdf'
    plt.savefig(filename)

    plt.close()

# Empty list to store the previous gene end position and scaffold
previous_gene_end = [None, None]
gene_list = []  
window_start, window_end = None, None  
plot_number = 1  

# Loop over each common gene
for i, gene in common_genes_sorted.iterrows():
    if previous_gene_end[1] == gene['scaffold'] and abs(previous_gene_end[0] - gene['adjusted_start_pos']) < window_size:
        window_end = max(window_end, gene['adjusted_end_pos'] + window_size / 2)
        gene_list.append((i, gene['gene_name']))  
    else:  
        if window_start is not None and window_end is not None:  
            plot_window(window_start, window_end, gene_list, SOR_df_sorted, VES_df_sorted)
            plot_number += 1
        window_start = max(0, gene['adjusted_position'] - window_size / 2)
        window_end = gene['adjusted_position'] + window_size / 2
        gene_list = [(i, gene['gene_name'])]  

    # Update previous position and scaf
    previous_gene_end = [gene['adjusted_end_pos'], gene['scaffold']]

# Plot last window after the loop
plot_window(window_start, window_end, gene_list, SOR_df_sorted, VES_df_sorted)
