import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import gc
import sys
import os

if len(sys.argv) != 10:
    print("Usage: python motif_cooc_analysis_updated_wip.py <motif_translation> <fimo output> "
          "<peaks bed file>  <iterations> <summary_excel_path> <plot title> <plot file> "
          "<peak calling approach (CR/MACS)> <species>")
    sys.exit(1)

# paths to the required files
motifs_translation_path = sys.argv[1]
fimo_path = sys.argv[2]
peaks_bed_path = sys.argv[3]
iterations = int(sys.argv[4])
summary_excel_path = sys.argv[5]
plot_title_sys = sys.argv[6]
plot_file_sys = sys.argv[7]
peak_calling_approach = sys.argv[8]
species = sys.argv[9]

plot_title_sys = plot_title_sys.replace("\\n", "\n")
motifs_transl = pd.read_excel(motifs_translation_path)

# number of comment rows that needs to be skipped in peaks.bed file
if peak_calling_approach == 'CR':
    if species in ['abeo', 'abeoforma']:
        skipr = 636
    elif species in ['capsa', 'capsaspora']:
        skipr = 30
    elif species in ['sub', 'suberites', 'sponge']:
        skipr = 29
elif peak_calling_approach == 'MACS':
    skipr = 0

# jaspar IDs of mapped TFs that were selected for the analysis
mapped_motifs_68 = ["MA0147.4_comb", "MA0147.4_all", "MA0147.4_exact", "MA0147.4", "MA0595.1", "MA0596.1", "MA0613.1",
                    "MA0486.2", "MA0663.1", "MA0669.1", "MA0083.3", "MA0691.1", "MA0737.1", "MA0758.1", "MA0770.1",
                    "MA0105.4", "MA0795.1", "MA0807.1", "MA0823.1", "MA0838.1", "MA0849.1", "MA0863.1", "MA1489.1",
                    "MA1570.1", "MA0685.2", "MA0798.3", "MA0831.3", "MA1511.2", "MA2325.1", "MA2328.1", "MA0619.2",
                    "MA0853.2", "MA0874.2", "MA1632.2", "MA1466.2", "MA1636.2", "MA0018.5", "MA0839.2", "MA0609.3",
                    "MA0754.3", "MA0639.2", "MA0469.4", "MA0154.5", "MA0162.5", "MA0598.4", "MA0828.3", "MA0098.4",
                    "MA0645.2", "MA0156.4", "MA0492.2", "MA0491.3", "MA0846.2", "MA0851.2", "MA0852.3", "MA1607.2",
                    "MA0593.2", "MA0037.5", "MA0143.5", "MA1990.2", "MA0647.2", "MA0131.3", "MA1991.2", "MA0046.3",
                    "MA0485.3", "MA0050.4", "MA0051.2", "MA0493.3", "MA0657.2", "MA1513.2", "MA1515.2", "MA1516.2",
                    "MA0039.5", "MA0768.3", "MA1518.3", "MA0659.4", "MA0052.5", "MA0497.2", "MA0620.4", "MA0664.2",
                    "MA1642.2", "MA0502.3", "MA1644.2", "MA0063.3", "MA0122.4", "MA0505.3", "MA0506.3", "MA0067.3",
                    "MA0070.2", "MA0782.3", "MA0784.3", "MA0799.3", "MA0002.3", "MA1118.2", "MA1153.2", "MA1562.2",
                    "MA0868.3", "MA0829.3", "MA1625.2", "MA0804.2", "MA0688.2", "MA1648.2", "MA0090.4", "MA1968.2",
                    "MA1122.2", "MA0861.2", "MA0526.5", "MA1627.2", "MA0095.4", "MA0749.2", "MA0752.2", "MA0147.NC_1",
                    "MA0147.NC_3", "MA0599.1", "MA1517.2", "MA1959.2", "MA1107.3", "MA1512.2", "MA0742.2", "MA0740.2",
                    "MA0741.1", "MA1514.2"]


# read fimo.tsv file and check if it is empty
def fimo_tsv(fimo_path, motif_id):
    dtype_mapping = {
        'motif_alt_id': 'category',
        'sequence_name': 'category',
        'start': 'uint32',
        'stop': 'uint32'
    }

    fimo_path = fimo_path + '/' + motif_id + "/fimo.tsv"
    file_size = os.path.getsize(fimo_path)
    if file_size < 500:
        return
    else:  # read file and skip the first column
        motif_tsv = pd.read_csv(fimo_path, sep='\t', usecols=range(1, 5), skipfooter=3,
                                engine='python', dtype=dtype_mapping)

        motif_tsv = motif_tsv.rename(columns={'sequence_name': 'chrom', 'start': 'start_m'})

        motif_tsv['chrom'] = motif_tsv["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')

        return motif_tsv


# get motif names based on the motif ids
def motif_name(motif_id, motifs_transl):
    motif_name = motifs_transl.loc[motifs_transl['matrix_id'] == motif_id, 'name']
    motif_name = ''.join(motif_name.astype(str))
    return motif_name


# get motif names and ids only for motifs that were successfully mapped by FIMO
def motif_id_name(mapped_motifs, motifs_transl, fimo_path):
    motif_ids = []
    motif_names = []

    for iid in range(len(mapped_motifs)):
        motif_id = mapped_motifs[iid]

        motif_n = motif_name(motif_id, motifs_transl)
        fimo_file = os.path.join(fimo_path, motif_id, "fimo.tsv")
        if os.path.getsize(fimo_file) < 500:
            continue
        else:
            motif_ids.append(motif_id)
            motif_names.append(motif_n)

    return motif_ids, motif_names


def motif_summary_update(peaks_bed, updated_df, motif_tsv, motif_name_t):

    merged = pd.merge(peaks_bed, motif_tsv, on='chrom')

    # merge and drop unnecessary columns
    merged = merged.loc[
        (merged['start_m'] >= merged['start']) &
        (merged['stop'] <= merged['end']),
        ['peak_code', 'motif_alt_id']
    ]

    # group df by peak_code and count overlapping elements
    peak_element_counts = merged.groupby(['peak_code', 'motif_alt_id']).size().reset_index(name=motif_name_t)
    updated_df_n = pd.merge(updated_df, peak_element_counts[['peak_code', motif_name_t]], on='peak_code', how='left')
    updated_df_n.iloc[:, -1] = updated_df_n.iloc[:, -1].fillna(0).astype('int32')

    # free up memory
    del merged, peak_element_counts
    gc.collect()
    return updated_df_n


# clean up tsv of the second motif by removing mapped motifs overlapping with any in the 1st TF tsv
def remove_overlap(tsv_fixed, tsv_check):
    tsv_check['code_check'] = np.arange(1, len(tsv_check) + 1, dtype=np.uint32)
    tsv_fixed['code_fix'] = np.arange(1, len(tsv_fixed) + 1, dtype=np.uint32)

    merge_cols_check = ['chrom', 'start_m', 'stop', 'code_check']
    merge_cols_fix = ['chrom', 'start_m', 'stop', 'code_fix']
    merged = pd.merge(
        tsv_fixed[merge_cols_fix].rename(columns={'start_m': 'start_m_fix', 'stop': 'stop_fix'}),
        tsv_check[merge_cols_check].rename(columns={'start_m': 'start_m_check', 'stop': 'stop_check'}),
        on='chrom',
        how='inner'
    )
    overlap_condition = ((merged['start_m_check'] <= merged['stop_fix']) &
                         (merged['stop_check'] >= merged['start_m_fix']))

    overlapping_rows = merged[overlap_condition]
    del merged

    tsv_check = tsv_check[~tsv_check['code_check'].isin(overlapping_rows['code_check'].unique())]
    tsv_check = tsv_check.drop(columns='code_check')
    return tsv_check


peaks_bed = pd.read_csv(peaks_bed_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skipr)
# adjust coordinates for indexing
peaks_bed['start'] = peaks_bed['start'] + 1
peaks_bed['end'] = peaks_bed['end'] + 1
#######
peaks_bed['chrom'] = peaks_bed["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')

peak_codes = [i for i in range(1, len(peaks_bed) + 1)]  # 1-based indexing
peaks_bed.insert(len(peaks_bed.columns), 'peak_code', peak_codes)

mapped_ids, mapped_names = motif_id_name(mapped_motifs_68, motifs_transl, fimo_path)

dimensions = len(mapped_names)
peaks_overlap_initial = np.zeros((dimensions, dimensions))
p_val_m = np.zeros((dimensions, dimensions))


for row in range(len(mapped_ids)):
    cooc_dict = {}  # dictionary that would contain all calculated co-oc values. Will be used for the p-val calculation
    # get the motif name and id and the tsv dataframe
    motif_1_id = mapped_ids[row]
    motif_1 = mapped_names[row]
    motif_1_tsv = fimo_tsv(fimo_path, motif_1_id)
    motif_1_tsv['code_fix'] = np.arange(1, len(motif_1_tsv) + 1, dtype=np.uint32)
    motif_1_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=peaks_bed,
                                           motif_tsv=motif_1_tsv, motif_name_t=motif_1)

    motif2_n = motif_1 + '_check'

    for col in range(len(mapped_ids)):
        motif_2_id = mapped_ids[col]
        motif_2 = mapped_names[col]
        motif_2_tsv = fimo_tsv(fimo_path, motif_2_id)
        if motif_2_tsv.empty:
            peaks_overlap_initial[row, col] = 0
            continue

        if motif_1 == motif_2:
            motif_2 = motif2_n
            motif_2_summary = motif_1_summary.copy()
            motif_2_summary[motif_2] = np.zeros(len(motif_1_summary), dtype='int32')
            motif_1_occurrence = 0

        else:
            # remove the overlap
            motif_2_tsv_clean = remove_overlap(tsv_fixed=motif_1_tsv, tsv_check=motif_2_tsv)

            motif_1_occurrence = (motif_1_summary[motif_1] > 0).sum()

            motif_2_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=motif_1_summary,
                                                   motif_tsv=motif_2_tsv_clean, motif_name_t=motif_2)
            del motif_2_tsv
            del motif_2_tsv_clean
            gc.collect()

        cooccur_val = ((motif_2_summary[motif_1] > 0) & (motif_2_summary[motif_2] > 0)).sum()
        peaks_overlap_initial[row, col] = cooccur_val

        # run random sampling
        for _ in range(iterations):
            randomized_sub_df = motif_2_summary.sample(n=motif_1_occurrence)

            if motif_2 in cooc_dict:
                cooc_dict[motif_2].append((randomized_sub_df[motif_2] > 0).sum())
            else:
                cooc_dict[motif_2] = [(randomized_sub_df[motif_2] > 0).sum()]

        p_value = sum(1 for peaks in cooc_dict[motif_2] if
                      peaks >= (peaks_overlap_initial[row, col])) / len(cooc_dict[motif_2])
        p_val_m[row, col] = p_value


peaks_overlap_initial_df = pd.DataFrame(peaks_overlap_initial, columns=mapped_names, index=mapped_names)
p_val_df = pd.DataFrame(p_val_m, columns=mapped_names, index=mapped_names)
# save results table
p_val_df.to_excel(summary_excel_path, index=False)


# Plotting and saving the heatmap
fig = plt.figure(figsize=(59, 59))
cmap = sns.color_palette("viridis_r", as_cmap=True)
ax = sns.heatmap(p_val_df, cmap=cmap, cbar_kws={"label": "Values in p_val_df"}, annot=True, fmt=".2f", vmin=0,
                 vmax=0.075, linewidths=0, rasterized=True)

plt.title(plot_title_sys, pad=30, fontsize=38)

fig.savefig(plot_file_sys, dpi=fig.dpi, bbox_inches='tight')
