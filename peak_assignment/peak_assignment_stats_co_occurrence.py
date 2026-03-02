import pandas as pd
import numpy as np
import statistics
import math
import os
import gc
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 8:
    print("Usage: peak_assigning_stat.py <full_peak_assignment_file> <motifs_translation_file> <fimo_path> "
          "<peaks_bed_path> <summary_output_file> <co-occurence table> <p-val cutoff>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
peak_assignment_sys = sys.argv[1]
motifs_transl = sys.argv[2]
fimo_path = sys.argv[3]
peaks_bed_path = sys.argv[4]
summary = sys.argv[5]
co_oc_table_path = sys.argv[6]
pval = int(sys.argv[7])



# get significant combinations from the co-occurrence table
def significant_combi(co_oc_table_path, max_p_val): # overlap_table_path,
    co_oc_table = pd.read_excel(co_oc_table_path)
    co_oc_table = co_oc_table.set_axis(list(co_oc_table.columns))
    positions = set(co_oc_table[co_oc_table <= max_p_val].stack().index.tolist())
    unique_positions = set(tuple(sorted(pos)) for pos in positions)

    return unique_positions



def fimo_tsv(fimo_path, motif_id):
    dtype_mapping = {
        'motif_alt_id': 'category',
        'sequence_name': 'category',
        'start': 'uint32',
        'stop': 'uint32'
    }

    fimo_path = fimo_path + '/' + motif_id + "/fimo.tsv"
    # file_size = os.path.getsize(fimo_path)
    # if file_size < 500:
    #     return  # motif_tsv
    # else:  # Remove the first column
    motif_tsv = pd.read_csv(fimo_path, sep='\t', usecols=range(1, 5), skipfooter=3,
                            engine='python', dtype=dtype_mapping)

    motif_tsv = motif_tsv.rename(columns={'sequence_name': 'chrom', 'start': 'start_m'})

    motif_tsv['chrom'] = motif_tsv["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')

    return motif_tsv


def motif_name(motif_id, motifs_transl):
    motif_name = motifs_transl.loc[motifs_transl['matrix_id'] == motif_id, 'name']
    motif_name = ''.join(motif_name.astype(str))
    return motif_name


def motif_id_name(mapped_motifs, motifs_transl, fimo_path, id_list=False):
    motif_ids = []
    motif_names = []
    if id_list:
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
    else:
        motif_id = mapped_motifs
        print(motif_id)
        motif_n = motif_name(motif_id, motifs_transl)
        return motif_n




def motif_summary_update(peaks_bed, updated_df, motif_tsv, motif_name_t):

    merged = pd.merge(peaks_bed, motif_tsv, on='chrom')

    merged = merged.loc[
        (merged['start_m'] >= merged['start']) &
        (merged['stop'] <= merged['end']),
        ['peak_code', 'motif_alt_id']  # Drop unnecessary cols
    ]

    # Group by peak and count overlapping elements
    peak_element_counts = merged.groupby(['peak_code', 'motif_alt_id'], observed=False).size().reset_index(name=motif_name_t)
    updated_df_n = pd.merge(updated_df, peak_element_counts[['peak_code', motif_name_t]], on='peak_code', how='left')
    updated_df_n.iloc[:, -1] = updated_df_n.iloc[:, -1].fillna(0).astype('int32')

    del merged, peak_element_counts  # Free memory
    gc.collect()
    return updated_df_n


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
    # tsv_fixed = tsv_fixed[~tsv_fixed['code_fix'].isin(overlapping_rows['code_fix'].unique())]

    # Drop 'code' columns to restore original structure
    tsv_check = tsv_check.drop(columns='code_check')

    return tsv_check


def get_signif_peaks(motif_1_id, motif_2_id, motif_2, fimo_path, peaks_bed_path, skipr, save_f, save_path,  for_stats, motifs_transl): #,motif_1,
    motifs_transl = pd.read_excel(motifs_transl)
    peaks_bed_orig = pd.read_csv(peaks_bed_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skipr)
    peak_codes = [i for i in range(1, len(peaks_bed_orig) + 1)]  # 1-based indexing
    peaks_bed_orig.insert(len(peaks_bed_orig.columns), 'peak_code', peak_codes)
    peaks_bed = peaks_bed_orig.copy()
    peaks_bed['start'] = peaks_bed['start'] + 1
    peaks_bed['end'] = peaks_bed['end'] + 1
    peaks_bed['chrom'] = peaks_bed["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')

    motif_1_tsv = fimo_tsv(fimo_path, motif_1_id)
    motif_1_tsv['code_fix'] = np.arange(1, len(motif_1_tsv) + 1, dtype=np.uint32)
    motif_1 = motif_id_name(mapped_motifs=motif_1_id, motifs_transl = motifs_transl, fimo_path=fimo_path)
    motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=peaks_bed,
                                           motif_tsv=motif_1_tsv, motif_name_t=motif_1)
    if motif_2:

        motif_2 = motif_id_name(mapped_motifs=motif_2_id, motifs_transl = motifs_transl, fimo_path=fimo_path)
        motif_2_tsv = fimo_tsv(fimo_path, motif_2_id)
        motif_2_tsv_clean= remove_overlap(tsv_fixed=motif_1_tsv, tsv_check=motif_2_tsv)

        motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=motif_summary,
                                               motif_tsv=motif_2_tsv_clean, motif_name_t=motif_2)
        cooccur_val = ((motif_summary[motif_1] > 0) & (motif_summary[motif_2] > 0)).sum()
        print('Number of peaks where motif 1 and 2 co-occur: ', cooccur_val)
        filtered_motif_summary = motif_summary[(motif_summary[motif_1] > 0) & (motif_summary[motif_2] > 0)]

    else:
        filtered_motif_summary = motif_summary[(motif_summary[motif_1] > 0) ]
    relevant_peaks = filtered_motif_summary['peak_code'].to_list()
    if save_f:
        relevant_peaks_bed = peaks_bed_orig[peaks_bed_orig['peak_code'].isin(relevant_peaks)]
        if for_stats:
            relevant_peaks_bed[['peak_code', 'chrom', 'start', 'end']].to_csv(save_path, index=False, header=True, sep='\t')
        else:
            relevant_peaks_bed[['chrom', 'start', 'end']].to_csv(save_path, index=False, header=False, sep='\t')
    print("Done!")
    return filtered_motif_summary, motif_1

def get_vals(rand_df, abs_vals, orig_val = False):
    # get absolute values
    categories = ['TSS', 'intragenic', 'proximal_intergenic', 'distal_intergenic']
    category_per = rand_df['category'].value_counts().reindex(categories).reset_index()
    category_per.columns = ['category', 'Frequency']
    counts_abs = category_per.set_index('category')['Frequency'].to_dict()

    if orig_val:
        return counts_abs['TSS'], counts_abs['intragenic'], counts_abs['proximal_intergenic'], counts_abs['distal_intergenic']
    else:
        abs_vals['TSS'].append(counts_abs['TSS'])
        abs_vals['intragenic'].append(counts_abs['intragenic'])
        abs_vals['proximal_intergenic'].append(counts_abs['proximal_intergenic'])
        abs_vals['distal_intergenic'].append(counts_abs['distal_intergenic'])
        return abs_vals

def get_p_val_above(values, orig_v, val_dict, categ):
    val_dict[categ].append(round(sum([1 for x in values if x >= orig_v])/len(values), 4))
    return val_dict


def get_p_val_below(values, orig_v, val_dict, categ):
    values = [0 if (isinstance(v, float) and math.isnan(v)) else v for v in values]
    val_dict[categ].append(round(sum([1 for x in values if x <= orig_v])/len(values), 4))
    val_dict[categ].append(round(statistics.mean(values), 2))

    val_dict[categ].append(round(statistics.stdev(values), 2))
    return val_dict

def plot_stats(peak_assignment_path, selected_peaks, iterations, summary_dict):
    #  read input files and rename the categories to fit the three category approach
    peak_assignment = pd.read_csv(peak_assignment_path, sep='\t')
    peak_assignment['category'] = peak_assignment['category'].replace(['proximal intergenic'], 'proximal_intergenic')
    peak_assignment['category'] = peak_assignment['category'].replace(['distal intergenic'], 'distal_intergenic')

    # get the number of peaks where a motif occurs or co-occurs with another motif
    peak_codes = selected_peaks['peak_code'].to_list()
    selected_peaks = peak_assignment[peak_assignment['peak_code'].isin(peak_codes)]
    peaks_count = len(selected_peaks)

    abs_vals={ 'TSS': [],
               'intragenic': [],
               'proximal_intergenic': [],
               'distal_intergenic': []}

    orig_abs_vals={'param': ['value', 'pval_above', 'pval_below', 'mean', 'stdev'],
                   'TSS': [0],
                   'intragenic': [0],
                   'proximal_intergenic': [0],
                   'distal_intergenic': [0]}

    # store the original absolute values for each category from the selected peaks
    orig_abs_vals['TSS'][0], orig_abs_vals['intragenic'][0], orig_abs_vals['proximal_intergenic'][0], orig_abs_vals['distal_intergenic'][0] = get_vals(selected_peaks, orig_abs_vals, orig_val=True)

    # perform random sampling of all assigned peaks
    for _ in range(iterations):
        subset = peak_assignment.sample(n=peaks_count)
        abs_vals = get_vals(subset, abs_vals)


    orig_abs_vals = get_p_val_above(abs_vals['TSS'], orig_abs_vals['TSS'][0], orig_abs_vals, 'TSS')
    orig_abs_vals = get_p_val_above(abs_vals['intragenic'], orig_abs_vals['intragenic'][0], orig_abs_vals, 'intragenic')
    orig_abs_vals = get_p_val_above(abs_vals['proximal_intergenic'], orig_abs_vals['proximal_intergenic'][0], orig_abs_vals, 'proximal_intergenic')
    orig_abs_vals = get_p_val_above(abs_vals['distal_intergenic'], orig_abs_vals['distal_intergenic'][0], orig_abs_vals, 'distal_intergenic')


    orig_abs_vals = get_p_val_below(abs_vals['TSS'], orig_abs_vals['TSS'][0], orig_abs_vals, 'TSS')
    orig_abs_vals = get_p_val_below(abs_vals['intragenic'], orig_abs_vals['intragenic'][0], orig_abs_vals, 'intragenic')
    orig_abs_vals = get_p_val_below(abs_vals['proximal_intergenic'], orig_abs_vals['proximal_intergenic'][0], orig_abs_vals, 'proximal_intergenic')
    orig_abs_vals = get_p_val_below(abs_vals['distal_intergenic'], orig_abs_vals['distal_intergenic'][0], orig_abs_vals, 'distal_intergenic')

    # update values
    # absolute values
    summary_dict['TSS_abs_val'].append(orig_abs_vals['TSS'][0])
    summary_dict['intragenic_abs_val'].append(orig_abs_vals['intragenic'][0])
    summary_dict['proximal_intergenic_abs_val'].append(orig_abs_vals['proximal_intergenic'][0])
    summary_dict['distal_intergenic_abs_val'].append(orig_abs_vals['distal_intergenic'][0])

    # pval above
    summary_dict['TSS_pval_above'].append(orig_abs_vals['TSS'][1])
    summary_dict['intragenic_pval_above'].append(orig_abs_vals['intragenic'][1])
    summary_dict['proximal_intergenic_pval_above'].append(orig_abs_vals['proximal_intergenic'][1])
    summary_dict['distal_intergenic_pval_above'].append(orig_abs_vals['distal_intergenic'][1])

    # pval below
    summary_dict['TSS_pval_below'].append(orig_abs_vals['TSS'][2])
    summary_dict['intragenic_pval_below'].append(orig_abs_vals['intragenic'][2])
    summary_dict['proximal_intergenic_pval_below'].append(orig_abs_vals['proximal_intergenic'][2])
    summary_dict['distal_intergenic_pval_below'].append(orig_abs_vals['distal_intergenic'][2])

    # mean
    summary_dict['TSS_mean'].append(orig_abs_vals['TSS'][3])
    summary_dict['intragenic_mean'].append(orig_abs_vals['intragenic'][3])
    summary_dict['proximal_intergenic_mean'].append(orig_abs_vals['proximal_intergenic'][3])
    summary_dict['distal_intergenic_mean'].append(orig_abs_vals['distal_intergenic'][3])

    # sd
    summary_dict['TSS_sd'].append(orig_abs_vals['TSS'][4])
    summary_dict['intragenic_sd'].append(orig_abs_vals['intragenic'][4])
    summary_dict['proximal_intergenic_sd'].append(orig_abs_vals['proximal_intergenic'][4])
    summary_dict['distal_intergenic_sd'].append(orig_abs_vals['distal_intergenic'][4])

    return summary_dict


all_motifs = ["MA0151.1", "MA0595.1", "MA0596.1", "MA0613.1", "MA0486.2", "MA0663.1", "MA0669.1", "MA0083.3",
              "MA0691.1", "MA0737.1", "MA0758.1", "MA0770.1", "MA0105.4", "MA0795.1", "MA0807.1", "MA0823.1",
              "MA0838.1", "MA0849.1", "MA0863.1", "MA1489.1", "MA1570.1", "MA0685.2", "MA0798.3", "MA0831.3",
              "MA1511.2", "MA2325.1", "MA2328.1", "MA0619.2", "MA0853.2", "MA0874.2", "MA1632.2", "MA1466.2",
              "MA0875.2", "MA0876.2", "MA1636.2", "MA0018.5", "MA0839.2", "MA0609.3", "MA0754.3", "MA0639.2",
              "MA1481.2", "MA0469.4", "MA0154.5", "MA0162.5", "MA0598.4", "MA0828.3", "MA0058.4", "MA0098.4",
              "MA0645.2", "MA0156.4", "MA0492.2", "MA0491.3", "MA0846.2", "MA0851.2", "MA0852.3", "MA1607.2",
              "MA0593.2", "MA0037.5", "MA0143.5", "MA1990.2", "MA0647.2", "MA1106.2", "MA0131.3", "MA1991.2",
              "MA0046.3", "MA0485.3", "MA0050.4", "MA0051.2", "MA0914.2", "MA0493.3", "MA0657.2", "MA1513.2",
              "MA1515.2", "MA1516.2", "MA0039.5", "MA0768.3", "MA1518.3", "MA0704.2", "MA0659.4", "MA0147.4",
              "MA0052.5", "MA0497.2", "MA0620.4", "MA0664.2", "MA0708.3", "MA1642.2", "MA0670.2", "MA0502.3",
              "MA1644.2", "MA0063.3", "MA0122.4", "MA0675.2", "MA0164.2", "MA0505.3", "MA0506.3", "MA0067.3",
              "MA0070.2", "MA0782.3", "MA0784.3", "MA1116.2", "MA0799.3", "MA0002.3", "MA1118.2", "MA1153.2",
              "MA1562.2", "MA0868.3", "MA0829.3", "MA1625.2", "MA0108.3", "MA0804.2", "MA0688.2", "MA0521.3",
              "MA1648.2", "MA0090.4", "MA1968.2", "MA1122.2", "MA0861.2", "MA0526.5", "MA1627.2", "MA0095.4",
              "MA0749.2", "MA0752.2"]

motifs_transl_exel = pd.read_excel(motifs_transl)
mapped_motif_id, mapped_motif_name = motif_id_name(all_motifs, motifs_transl=motifs_transl_exel, fimo_path=fimo_path,
                                                   id_list=True)
translation = dict(zip(mapped_motif_name, mapped_motif_id))

motifs_enrichment_summ = {
    'motif_name': [],
    'motif_id': [],
    'TSS_abs_val': [],
    'TSS_pval_above': [],
    'TSS_pval_below': [],
    'TSS_mean': [],
    'TSS_sd': [],
    'intragenic_abs_val': [],
    'intragenic_pval_above': [],
    'intragenic_pval_below': [],
    'intragenic_mean': [],
    'intragenic_sd': [],
    'proximal_intergenic_abs_val': [],
    'proximal_intergenic_pval_above': [],
    'proximal_intergenic_pval_below': [],
    'proximal_intergenic_mean': [],
    'proximal_intergenic_sd': [],
    'distal_intergenic_abs_val': [],
    'distal_intergenic_pval_above': [],
    'distal_intergenic_pval_below': [],
    'distal_intergenic_mean': [],
    'distal_intergenic_sd': [],
}

sign_cooc = significant_combi(co_oc_table_path, max_p_val=pval)
sign_cooc = list(sign_cooc)


for pair in sign_cooc:
    pair = list(pair)

    for ii in range(2):
        if ii == 0:
            motif1, motif2 = pair[0], pair[1]
        else:
            motif2, motif1 = pair[0], pair[1]
        tfbm1 = translation[motif1]
        tfbm2 = translation[motif2]
        combination = motif1 + ', ' + motif2

        tf_summary, motif_name_j = get_signif_peaks(motif_1_id=tfbm1,
                                                    motif_2=True,
                                                    motif_2_id=tfbm2,
                                                    for_stats=True,
                                                    fimo_path=fimo_path,
                                                    peaks_bed_path=peaks_bed_path,
                                                    skipr=29,
                                                    save_f=False,
                                                    save_path=False,
                                                    motifs_transl=motifs_transl)
        motifs_enrichment_summ['motif_name'].append(motif_name_j)
        motifs_enrichment_summ['motif_id'].append(combination)
        motifs_enrichment_summ = plot_stats(peak_assignment_path=peak_assignment_sys,
                                            selected_peaks=tf_summary,
                                            iterations=10000,
                                            summary_dict=motifs_enrichment_summ)


comb_stats = pd.DataFrame.from_dict(motifs_enrichment_summ)
comb_stats.to_excel(summary,  index=False)