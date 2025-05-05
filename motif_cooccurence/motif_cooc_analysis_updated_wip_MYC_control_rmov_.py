import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# import pyranges as pr
import timeit
import sys
import os

if len(sys.argv) != 11:
    print("Usage: python motif_cooc_analysis_updated_wip.py <motif_translation> <fimo output> "
          "<peaks bed file> <fimo_summary_table> <iterations> <summary_excel_path> <plot title> <plot file> "
          "<peak calling approach (CR/MACS)> <species>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
motifs_translation_path = sys.argv[1]
fimo_path = sys.argv[2]
peaks_bed_path = sys.argv[3]
summary_table_path = sys.argv[4]
iterations = int(sys.argv[5])
summary_excel_path = sys.argv[6]
plot_title_sys = sys.argv[7]
plot_file_sys = sys.argv[8]
peak_calling_approach = sys.argv[9]
species = sys.argv[10]

plot_title_sys = plot_title_sys.replace("\\n", "\n")
motifs_transl = pd.read_excel(motifs_translation_path)

if peak_calling_approach == 'CR':
    if species in ['abeo', 'abeoforma']:
        skipr = 636
    elif species in ['capsa', 'capsaspora']:
        skipr = 30
    elif species in ['sub', 'suberites', 'sponge']:
        skipr = 114
elif peak_calling_approach == 'MACS':
    skipr = 0


mapped_motifs_68 = ["MA0670.2", "MA0798.3", "MA0754.3", "MA0807.1", "MA1991.2", "MA0664.2", "MA0151.1", "MA0752.2",
                    "MA0052.5", "MA0914.2", "MA0657.2", "MA0596.1", "MA1644.2", "MA0874.2", "MA0691.1", "MA0849.1",
                    "MA0002.3", "MA0863.1", "MA0784.3", "MA0799.3", "MA0852.3", "MA0846.2", "MA0095.4", "MA0838.1",
                    "MA0613.1", "MA0058.4", "MA0831.3", "MA0090.4", "MA1625.2", "MA1990.2", "MA1968.2", "MA1632.2",
                    "MA0639.2", "MA0626.2", "MA0704.2", "MA0839.2", "MA1466.2", "MA2328.1", "MA0876.2", "MA1122.2",
                    "MA0143.5", "MA0147.4", "MA1511.2", "MA0868.3", "MA0469.4", "MA0823.1", "MA0521.3", "MA1153.2",
                    "MA1607.2", "MA0051.2", "MA0861.2", "MA0685.2", "MA0505.3", "MA0770.1", "MA0659.4", "MA0506.3",
                    "MA0037.5", "MA0098.4", "MA0162.5", "MA0063.3", "MA1627.2", "MA0083.3", "MA0502.3", "MA1513.2",
                    "MA0105.4", "MA0593.2", "MA0708.3", "MA0875.2", "MA0147.4_comb", "MA0147.4_mut", "MA0147.4_mut2",
                    "MA0147.4_ex", "MA000.2", "MA0122.4", "MA0124.3"] #, "MA000.2"]

def fimo_tsv(fimo_path, motif_id):

    dtype_mapping = {
        # 'motif_id': 'category',
        'motif_alt_id': 'category',
        'sequence_name': 'category',
        'start': 'uint32',
        'stop': 'uint32'
    }

    # fimo_path = fimo_path + '\\' + motif_id + "\\fimo.tsv"
    fimo_path = fimo_path + '/' + motif_id + "/fimo.tsv"
    file_size = os.path.getsize(fimo_path)
    if file_size < 500:
        motif_tsv = pd.read_csv(fimo_path, sep='\t', skipfooter=3, engine='python')
        return motif_tsv
    else:  ### Remove the first column
        motif_tsv = pd.read_csv(fimo_path, sep='\t', usecols=range(1,5), skipfooter=3,
                                engine='python', dtype=dtype_mapping)

        motif_tsv = motif_tsv.rename(columns={'sequence_name': 'chrom', 'start': 'start_m'})
        #######
        motif_tsv['chrom'] == motif_tsv["chrom"].str.replace("sca", "", regex=True).astype('uint16')

        return motif_tsv


def motif_name(motif_id, motifs_transl):
    query_condition = "matrix_id=='" + motif_id + "'"
    motif_name = motifs_transl.query(query_condition)["name"]
    motif_name = ''.join(motif_name.astype(str))
    return motif_name


def motif_id_name(mapped_motifs, motifs_transl):
    motif_ids = []
    motif_names = []

    for iid in range(len(mapped_motifs)):
        motif_id = mapped_motifs[iid]

        motif_n = motif_name(motif_id, motifs_transl)
        motif_tsv_f = fimo_tsv(fimo_path, motif_id)
        if motif_tsv_f.empty:
            del motif_tsv_f
            continue
        else:
            motif_ids.append(motif_id)
            motif_names.append(motif_n)
            del motif_tsv_f


    return motif_ids, motif_names



def motif_summary_update(peaks_bed, motif_tsv, motif_name):

    merged = pd.merge(peaks_bed, motif_tsv, on='chrom')

    merged = merged.loc[
        (merged['start_m'] >= merged['start']) &
        (merged['stop'] <= merged['end']),
        ['peak_code', 'motif_alt_id']  # Drop unnecessary cols early
    ]
    # old approach
    # merged = pd.merge(peaks_bed, motif_tsv, on='chrom')
    #
    # merged = merged[(merged['start_m'] >= merged['start']) & (merged['stop'] <= merged['end'])]

    # Group by peak and count overlapping elements
    peak_element_counts = merged.groupby(['peak_code', 'motif_alt_id']).size().reset_index(name=motif_name)
    peaks_bed = pd.merge(peaks_bed, peak_element_counts[['peak_code', motif_name]], on='peak_code', how='left')
    peaks_bed.iloc[:, -1] = peaks_bed.iloc[:, -1].fillna(0).astype('int32')

    del merged  # Free memory
    return peaks_bed


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
    overlap_condition = ((merged['start_m_check'] <= merged['stop_fix']) & (merged['stop_check'] >= merged['start_m_fix']))

    overlapping_rows = merged[overlap_condition]
    del merged

    tsv_check = tsv_check[~tsv_check['code_check'].isin(overlapping_rows['code_check'].unique())]
    tsv_fixed = tsv_fixed[~tsv_fixed['code_fix'].isin(overlapping_rows['code_fix'].unique())]

    # Drop 'code' columns to restore original structure
    tsv_check = tsv_check.drop(columns='code_check')
    tsv_fixed = tsv_fixed.drop(columns='code_fix')

    return tsv_check, tsv_fixed


peaks_bed = pd.read_csv(peaks_bed_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skipr)
#######
peaks_bed['chrom'] == peaks_bed["chrom"].str.replace("sca", "", regex=True).astype('uint16')

peak_codes = [i for i in range(1, len(peaks_bed)+1)]  # 1-based indexing
peaks_bed.insert(len(peaks_bed.columns), 'peak_code', peak_codes)

mapped_ids, mapped_names = motif_id_name(mapped_motifs_68, motifs_transl)

dimensions = len(mapped_names)
peaks_overlap_initial = np.zeros((dimensions, dimensions))
p_val_m = np.zeros((dimensions, dimensions))

start = timeit.default_timer()




for row in range(len(mapped_ids)):
    cooc_dict = {}  # dictionary that would contain all calculated co-oc values. Will be used for the p-val calculation
    # get the motif name and id and the tsv dataframe
    motif_1_id = mapped_ids[row]
    motif_1 = mapped_names[row]
    motif_1_tsv = fimo_tsv(fimo_path, motif_1_id)

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
            motif_1_summary = pd.DataFrame(0, index=range(len(peaks_bed)), columns=[motif_1, motif_2])
            motif_1_occurrence = 0

        else:
            # remove the overlap
            motif_2_tsv_clean, motif_1_tsv_clean = remove_overlap(tsv_fixed=motif_1_tsv, tsv_check=motif_2_tsv)
            # _clean - no overlapping coordinates are present

            # trying to repeat the initial approach:
            motif_1_summary = motif_summary_update(peaks_bed=peaks_bed, motif_tsv=motif_1_tsv_clean, motif_name=motif_1)
            motif_1_occurrence = (motif_1_summary[motif_1] > 0).sum()

            motif_1_summary = motif_summary_update(peaks_bed=motif_1_summary, motif_tsv=motif_2_tsv, motif_name=motif_2)

        cooccur_val = ((motif_1_summary[motif_1] > 0) & (motif_1_summary[motif_2] > 0)).sum()
        peaks_overlap_initial[row, col] = cooccur_val

        for _ in range(iterations):
            randomized_sub_df = motif_1_summary.sample(n=motif_1_occurrence)

            if motif_2 in cooc_dict:
                cooc_dict[motif_2].append((randomized_sub_df[motif_2] > 0).sum())
            else:
                cooc_dict[motif_2] = [(randomized_sub_df[motif_2] > 0).sum()]

        p_value = sum(1 for peaks in cooc_dict[motif_2] if peaks >= (peaks_overlap_initial[row, col])) / len(cooc_dict[motif_2])
        p_val_m[row, col] = p_value


stop = timeit.default_timer()
print('Time: ', stop - start)

peaks_overlap_initial_df = pd.DataFrame(peaks_overlap_initial, columns=mapped_names, index=mapped_names)
p_val_df = pd.DataFrame(p_val_m, columns=mapped_names, index=mapped_names)
# save results
p_val_df.to_excel(summary_excel_path, index=False)


# Plot heatmap with a mask and color scheme
fig = plt.figure(figsize=(59, 59))
cmap = sns.color_palette("viridis_r", as_cmap=True)
ax = sns.heatmap(p_val_df, cmap=cmap, cbar_kws={"label": "Values in p_val_df"}, annot=True, fmt=".2f", vmin=0,
                 vmax=0.075, linewidths=0, rasterized=True)

# Display the plot
plt.title(plot_title_sys, pad=30, fontsize=38)

fig.savefig(plot_file_sys, dpi=fig.dpi, bbox_inches='tight')
