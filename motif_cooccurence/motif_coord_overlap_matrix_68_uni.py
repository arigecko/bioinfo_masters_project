import pandas as pd
import time
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 7:
    print("Usage: python script.py <fimo_path> <genome_fasta> <resulting_csv_path> <plot_title> "
          "<plot_path_grid> <plot_path_no_grid>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
fimo_path_sys = sys.argv[1]
genome_f_path_sys = sys.argv[2]
summary_csv_path_sys = sys.argv[3]
plot_title_sys = sys.argv[4]
plot_path_grid_sys = sys.argv[5]
plot_path_nogrid_sys = sys.argv[6]
# species = sys.argv[7]
#
# if species == 's':
#     f_scaffold = 'sca1'
# elif species == 'c':
#     f_scaffold = 'Cowc_Chr01'
# elif species == 'a':
#     f_scaffold = 'Abeoforma_whisleri_P_RNA_scaffold_1'

plot_title_sys = plot_title_sys.replace("\\n", "\n")


mapped_motifs_68 = ["MA0670.2", "MA0798.3", "MA0754.3", "MA0807.1", "MA1991.2", "MA0664.2", "MA0151.1", "MA0752.2",
                    "MA0052.5", "MA0914.2", "MA0657.2", "MA0596.1", "MA1644.2", "MA0874.2", "MA0691.1", "MA0849.1",
                    "MA0002.3", "MA0863.1", "MA0784.3", "MA0799.3", "MA0852.3", "MA0846.2", "MA0095.4", "MA0838.1",
                    "MA0613.1", "MA0058.4", "MA0831.3", "MA0090.4", "MA1625.2", "MA1990.2", "MA1968.2", "MA1632.2",
                    "MA0639.2", "MA0626.2", "MA0704.2", "MA0839.2", "MA1466.2", "MA2328.1", "MA0876.2", "MA1122.2",
                    "MA0143.5", "MA0147.4", "MA1511.2", "MA0868.3", "MA0469.4", "MA0823.1", "MA0521.3", "MA1153.2",
                    "MA1607.2", "MA0051.2", "MA0861.2", "MA0685.2", "MA0505.3", "MA0770.1", "MA0659.4", "MA0506.3",
                    "MA0037.5", "MA0098.4", "MA0162.5", "MA0063.3", "MA1627.2", "MA0083.3", "MA0502.3", "MA1513.2",
                    "MA0105.4", "MA0593.2", "MA0708.3", "MA0875.2"]


# create dictionary to later translate coordinates:
def new_coord_calc(genome_f_path):
    genome = SeqIO.to_dict(SeqIO.parse(genome_f_path, "fasta"))
    scaffolds = list(genome.keys())
    scaffold_lengths_added = [0]
    scaffoldlength_dict = {scaffolds[0]: 0}  # contains summed up lengths of all previous scaffolds (based on the key - sca3 : length of sca1+sca2)
    for i in range(1, len(scaffolds)):
        prev_len = len(genome[scaffolds[i-1]].seq)
        scaffold_lengths_added.append(scaffold_lengths_added[i-1] + prev_len)
        scaffoldlength_dict[scaffolds[i]] = scaffold_lengths_added[i-1] + prev_len
    # print(scaffoldlength_dict)
    genome_len = len(genome[scaffolds[0]].seq)
    scaffold_lengths_added = [0]
    scaffoldlength_dict = {scaffolds[0]: 0}  # contains summed up lengths of all previous scaffolds (based on the key - sca3 : length of sca1+sca2)
    for i in range(1, len(scaffolds)):
        prev_len = len(genome[scaffolds[i-1]].seq)
        current_len = len(genome[scaffolds[i]].seq)
        genome_len += current_len
        scaffold_lengths_added.append(scaffold_lengths_added[i-1] + prev_len)
        scaffoldlength_dict[scaffolds[i]] = scaffold_lengths_added[i-1] + prev_len

    return scaffoldlength_dict, genome_len


# access fimo tsv files and create one dataframe

def read_update_coord(mapped_motifs, fimo_path, genome_f_path):
    # get the values as to how coordinates should be adjusted based on the scafold
    scaffoldlength_dict, genome_len = new_coord_calc(genome_f_path)
    # future columns of the df summary for all the mapped motifs
    motifs_summary = {
        'scafold': [],
        'motif_id': [],
        'motif_name': [],
        'start': [],
        'end': [],
        'new_start': [],
        'new_end': []
    }

    # read corresponding tsv files and prepare data for the dic
    for motif_id in mapped_motifs:
        fimo_path_full = fimo_path + '/' + motif_id + "/fimo.tsv"
        motif_tsv = pd.read_csv(fimo_path_full, sep='\t')
        # remove 3 last rows from the df - don't contain any relevant information
        motif_tsv.drop(motif_tsv.tail(3).index, inplace=True)

        if motif_tsv.empty:
            # return peaks_bed_extended, peaks_dict
            continue

        for ind in motif_tsv.index:
            m_scaffold = motif_tsv['sequence_name'][ind]
            m_start, m_end = int(motif_tsv['start'][ind] - 1), int(motif_tsv['stop'][ind] - 1)
            scaffold_len = scaffoldlength_dict[m_scaffold]

            motifs_summary['scafold'].append(m_scaffold)
            motifs_summary['motif_id'].append(motif_id)
            motifs_summary['motif_name'].append(motif_tsv['motif_alt_id'][ind])
            motifs_summary['start'].append(m_start)
            motifs_summary['end'].append(m_end)
            motifs_summary['new_start'].append(m_start + scaffold_len)
            motifs_summary['new_end'].append(m_end + scaffold_len)

    # create summary df
    motifs_summary_df = pd.DataFrame.from_dict(motifs_summary)

    return motifs_summary_df


motifs_summary = read_update_coord(mapped_motifs_68, fimo_path_sys, genome_f_path_sys)

motifs_overlap_val = {}
# for each key(motif1) in the dictionary following values for all other values will be stored -
# proportion of mappings of motif1 that overlap with the motif2 (iterates through all available motifs)

present_motifs = motifs_summary['motif_name'].unique().tolist()


# Create a dictionary of DataFrames for each motif for faster access
motif_dfs = {motif: motifs_summary[motifs_summary['motif_name'] == motif] for motif in present_motifs}

start = time.time()
for motif_1 in present_motifs:
    motif_1_df = motif_dfs[motif_1]
    total_motif_1 = len(motif_1_df)
    overlap_fractions = []

    for motif_2 in present_motifs:
        if motif_1 == motif_2:
            overlap_fraction = 100
        else:
            motif_2_df = motif_dfs[motif_2]

            # Vectorized conditions to find overlaps
            condition1 = (motif_2_df['new_start'].values[:, None] <= motif_1_df['new_start'].values) & \
                         (motif_2_df['new_end'].values[:, None] >= motif_1_df['new_start'].values)
            condition2 = (motif_2_df['new_start'].values[:, None] <= motif_1_df['new_end'].values) & \
                         (motif_2_df['new_end'].values[:, None] >= motif_1_df['new_end'].values)

            overlap_count = (condition1 | condition2).any(axis=0).sum()
            overlap_fraction = overlap_count / total_motif_1 * 100

        overlap_fractions.append(overlap_fraction)

    motifs_overlap_val[motif_1] = overlap_fractions


end1 = time.time()
print("motif coordinates comparison: ", "%.2f" % (end1 - start), " sec")


# created a df from the collected data
resulting_df = pd.DataFrame.from_dict(motifs_overlap_val, orient='index')
resulting_df.columns = present_motifs

# save the dataframe to a csv file
resulting_df_csv = pd.DataFrame.from_dict(motifs_overlap_val, orient='index')
resulting_df_csv.columns = present_motifs
resulting_df_csv.insert(0, 'present_motifs', present_motifs)
resulting_df_csv.to_csv(summary_csv_path_sys, sep=",", index=False)

# save two plots: with and without grid

fig = plt.figure(figsize=(59, 59))
plt.title(plot_title_sys, fontsize=30)
sns.heatmap(resulting_df, annot=True, fmt=".2f", vmin=0, cmap='viridis_r', linewidths=0.005, linecolor='gray')
fig.savefig(plot_path_grid_sys, dpi=fig.dpi, bbox_inches='tight')

fig1 = plt.figure(figsize=(59, 59))
plt.title(plot_title_sys, fontsize=30)
sns.heatmap(resulting_df, annot=True, fmt=".2f", vmin=0, cmap='viridis_r')
fig1.savefig(plot_path_nogrid_sys, dpi=fig1.dpi, bbox_inches='tight')
