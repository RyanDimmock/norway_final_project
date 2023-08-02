import pandas as pd
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

# Function to extract chromosome number from scaffold column
def extract_chrom_number(scaffold):
    match = re.search(r'\d+$', scaffold)
    if match:
        return int(match.group())
    else:
        return None

# Load in data
SOR_df = pd.read_csv('50SNPs_final_allsites_SORvsNonspring.csv')
VES_df = pd.read_csv('50SNPs_final_allsites_VESvsNonspring.csv')
common_genes = pd.read_csv("common_spring_Gpos_with_scaf_unique.tsv", sep="\t")

# Extract chromosome number from the scaffold
SOR_df['CHROM_Number'] = SOR_df['scaffold'].apply(extract_chrom_number)
VES_df['CHROM_Number'] = VES_df['scaffold'].apply(extract_chrom_number)
common_genes['CHROM_Number'] = common_genes['Cexcelsa_scaf_'].apply(extract_chrom_number)

# Calculate max position for each chromosome and define ranges
max_positions = SOR_df.groupby('CHROM_Number')['midpoint'].max()
chrom_ranges = {chrom: (sum(max_positions[:chrom-1]), sum(max_positions[:chrom])) for chrom in range(1, max(SOR_df['CHROM_Number']) + 1)}

# Adjust the midpoint values
SOR_df['adjusted_midpoint'] = SOR_df.apply(lambda row: chrom_ranges[row['CHROM_Number']][0] + row['midpoint'], axis=1)
VES_df['adjusted_midpoint'] = VES_df.apply(lambda row: chrom_ranges[row['CHROM_Number']][0] + row['midpoint'], axis=1)
common_genes['adjusted_midpoint'] = common_genes.apply(lambda row: chrom_ranges[row['CHROM_Number']][0] + row['midpoint'], axis=1)

# Compute the threshold for the top 1% of Fst values
fst_threshold = SOR_df['Fst'].quantile(0.99)
new_fst_threshold = VES_df['Fst'].quantile(0.99)

# Convert fst to negative to create mirror image
VES_df['Fst'] = -VES_df['Fst']

# Set y limits and plot
y_max = 0.6
y_min = -0.4
plt.figure(figsize=(15, 8))

# Define colors for alternating between chromosomes
original_colors = ['blue', 'darkblue']
mirror_colors = ['orange', 'darkorange']

# Plot the original data with new colors and label
for chrom in range(1, max(SOR_df['CHROM_Number']) + 1):
    SOR_data = SOR_df[SOR_df['CHROM_Number'] == chrom]
    plt.scatter(SOR_data['adjusted_midpoint'], SOR_data['Fst'], s=1, color=original_colors[chrom % 2], label='SOR (SL2)' if chrom==1 else "", marker='o', zorder = 3)

plt.axhline(y=fst_threshold, color='purple', linestyle='--', linewidth = 0.7, label=f'Top 1% Fst values (SOR) ({fst_threshold:.4f})')

# Plot the new data as mirror image with new colors and label
for chrom in range(1, max(VES_df['CHROM_Number']) + 1):
    VES_data = VES_df[VES_df['CHROM_Number'] == chrom]
    plt.scatter(VES_data['adjusted_midpoint'], VES_data['Fst'], s=1, color=mirror_colors[chrom % 2], label='VES (SL1)' if chrom==1 else "", marker='o', zorder = 3)

plt.axhline(y=-new_fst_threshold, color='purple', linestyle='--', linewidth = 0.7, label=f'Top 1% Fst values (VES) ({new_fst_threshold:.4f})')

# Plot the common genes
for midpoint in common_genes['adjusted_midpoint']:
    plt.axvline(x=midpoint, color='lightgreen', alpha=0.3, linewidth=2.1, zorder = 2)

plt.axhline(y=0, color='black', linestyle='-', linewidth=1)


# grids and ticks
ticks = [chrom_ranges[chrom][0] + (chrom_ranges[chrom][1] - chrom_ranges[chrom][0]) / 2 for chrom in range(1, max(SOR_df['CHROM_Number']) + 1)]
gridlines = [chrom_ranges[chrom][i] for chrom in range(1, max(SOR_df['CHROM_Number']) + 1) for i in (0, 1)]
for line in gridlines:
    plt.axvline(x=line, color='gray', linestyle='--', linewidth=0.5)
labels = [str(chrom) for chrom in range(1, max(SOR_df['CHROM_Number']) + 1)]

# Set the ticks and labels
plt.xticks(ticks, labels)
plt.xlabel('Scaffold')
plt.ylabel('Fst')
plt.ylim(y_min, y_max)

# Using the absolute value for the y ticks and rounding to 1 decimal place
locs, _ = plt.yticks()
plt.yticks(locs, [f'{abs(loc):.1f}' for loc in locs])

plt.grid(False)

# Show the plot
plt.savefig("final_FST_Spring_vsNON_grn.pdf")
plt.show()
