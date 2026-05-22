import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import sys

if len(sys.argv) != 7:
    print("Usage: python motif_cooc_...py <motif_overlap_matrix> <fimo_summary_table_w_deneIDs> <iterations>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
motif_overlap_matrix_path = sys.argv[1]
summary_table_path = sys.argv[2]
repeats = int(sys.argv[3])
summary_excel_path = sys.argv[4]
plot_title_sys = sys.argv[5]
plot_file_sys = sys.argv[6]

plot_title_sys = plot_title_sys.replace("\\n", "\n")

# read files
motif_overlap_matrix = pd.read_csv(motif_overlap_matrix_path)
motif_overlap_matrix.set_index('present_motifs', inplace=True)

summary_table = pd.read_excel(summary_table_path)
# first 5 columns are 'chrom' 'start' 	'end' 'peak_code' and 'gene_ids' - not taken into account
all_motifs = list(summary_table.columns)[5:-1]  # all motifs that were attempted to be mapped

# get columns that were not mapped at all or empy columns - have sum of 0 over all rows
column_sums = summary_table[all_motifs].sum()
mapped_motifs_sum_zero = column_sums[column_sums == 0].index.tolist()

# drop these columns from the dataframe
summary_table_59 = summary_table.drop(columns=mapped_motifs_sum_zero)

# assign peak codes as indices
summary_table_59.set_index('peak_code', inplace=True)
mapped_motifs = list(summary_table_59.columns)[4:-1]

# remove all unnecessary columns
summary_table_59_short = summary_table_59.loc[:, mapped_motifs]

# transpose the table to cluster TF motifs with peaks as features and not peaks with TF motifs as features
summary_table_59_T = summary_table_59_short.transpose()
peak_codes = list(summary_table_59_T.columns)


# calculating total number of peaks that contain both (any pair of TFs???):
dimensions = len(mapped_motifs)
peaks_overlap_m = np.zeros((dimensions, dimensions))

peaks_overlap_initial = np.zeros((dimensions, dimensions))
p_val_m = np.zeros((dimensions, dimensions))


for row in range(peaks_overlap_m.shape[0]):
    cooc_dict = {}  # dictionary that would contain all calculated co-oc values. Will be used for the p-val calculation
    motif_1 = mapped_motifs[row]
    motif_1_occurrence = (summary_table[motif_1] > 0).sum()
    # collect initial values
    for col in range(peaks_overlap_m.shape[1]):
        motif_2 = mapped_motifs[col]
        cooccur_val = ((summary_table[motif_1] > 0) & (summary_table[motif_2] > 0)).sum()
        peaks_overlap_initial[row, col] = cooccur_val

    for _ in range(repeats):
        randomized_sub_df = summary_table.sample(n=motif_1_occurrence)
        for col in range(peaks_overlap_m.shape[1]):
            motif_2 = mapped_motifs[col]
            if motif_2 in cooc_dict:
                cooc_dict[motif_2].append((randomized_sub_df[motif_2] > 0).sum())
            else:
                cooc_dict[motif_2] = [(randomized_sub_df[motif_2] > 0).sum()]
    for col in range(peaks_overlap_m.shape[1]):
        key = list(cooc_dict.keys())[col]
        p_value = sum(1 for peaks in cooc_dict[key] if peaks >= (peaks_overlap_initial[row, col])) / len(cooc_dict[key])
        p_val_m[row, col] = p_value

peaks_overlap_initial_df = pd.DataFrame(peaks_overlap_initial, columns=mapped_motifs, index=mapped_motifs)
p_val_df = pd.DataFrame(p_val_m, columns=mapped_motifs, index=mapped_motifs)

# save results
p_val_df.to_excel(summary_excel_path, index=False)


# Reorder  motif_overlap_matrix_reord to match the structure of p_val_m (same index and columns)
motif_overlap_matrix_reord = motif_overlap_matrix.loc[p_val_df.index, p_val_df.columns]

# Create a mask where motif_overlap_matrix_reord > 5 -> more than 5% of the coordinates overlap between TFs
mask = motif_overlap_matrix_reord > 5

# Create a custom colormap that adds gray for the masked values
cmap = sns.color_palette("viridis_r", as_cmap=True)

# Set a colormap for the masked values (gray)
gray_cmap = ListedColormap(['gray'])

# Plot heatmap with a mask and color scheme
fig = plt.figure(figsize=(59, 59))
ax = sns.heatmap(p_val_df, cmap=cmap, mask=mask, cbar_kws={"label": "Values in p_val_df"}, annot=True, fmt=".2f", vmin=0,
                 vmax=0.075, linewidths=0, rasterized=True)

# Overlay grayed-out cells
sns.heatmap(p_val_df, mask=~mask, cmap=gray_cmap, cbar=False, linewidths=0, rasterized=True)

# Display the plot
plt.title(plot_title_sys, pad=30, fontsize=38)

fig.savefig(plot_file_sys, dpi=fig.dpi, bbox_inches='tight')
