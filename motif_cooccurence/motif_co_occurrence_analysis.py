import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import gc
import sys
import os

if len(sys.argv) != 11:
    print("Usage: python motif_cooc_analysis_updated_wip.py <motif_translation> <fimo output> "
          "<peaks bed file>  <iterations> <summary_excel_path> <plot title> <plot file> "
          "<peak calling approach (CR/MACS)> <species> <<pwm_ids_txt>>")
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
selected_motifs_path = sys.arg[10]

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

# read jaspar IDs of mapped TFs that were selected for the analysis
with open(selected_motifs_path, 'r') as file:
    mapped_motifs = [line.strip() for line in file]


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

mapped_ids, mapped_names = motif_id_name(mapped_motifs, motifs_transl, fimo_path)

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
fig = plt.figure(figsize=(89, 84))
cmap = sns.color_palette("viridis_r", as_cmap=True)
ax = sns.heatmap(p_val_df, cmap=cmap, cbar_kws={"label": "Values in p_val_df"}, vmin=0,
                 vmax=0.005, linewidths=0, rasterized=True)

plt.title(plot_title_sys, pad=30, fontsize=50)

fig.savefig(plot_file_sys, dpi=fig.dpi, bbox_inches='tight')
