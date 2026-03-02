import pandas as pd
import numpy as np
import gc
import sys
import os

if len(sys.argv) != 8:
    print("Usage: create_fimo_summary_gene_ID.py <peaks_bed> <motifs_translation> <fimo_output> "
          "<peak_calling_approach> <species> <output_excel_file> <peak_assignment_tsv_path>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
peaks_bed_path = sys.argv[1]
motifs_translation_path = sys.argv[2]
fimo_output_path = sys.argv[3]
peak_calling_approach = sys.argv[4]
species = sys.argv[5]
excel_file_path = sys.argv[6]
peak_encoding_path_sys = sys.argv[7]

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

# jaspar IDs of mapped TFs
mapped_motifs_set = ['MA0002.3', 'MA0063.3', 'MA0122.4', 'MA0151.1', 'MA0492.2', 'MA0595.1', 'MA0645.2', 'MA0685.2',
                     'MA0749.2', 'MA0798.3', 'MA0839.2', 'MA0874.2', 'MA1153.2', 'MA1516.2', 'MA1636.2', 'MA2328.1',
                     'MA0018.5', 'MA0067.3', 'MA0131.3', 'MA0154.5', 'MA0493.3', 'MA0596.1', 'MA0647.2', 'MA0688.2',
                     'MA0752.2', 'MA0799.3', 'MA0846.2', 'MA0875.2', 'MA1466.2', 'MA1517.2', 'MA0037.5', 'MA0070.2',
                     'MA0143.5', 'MA0156.4', 'MA0497.2', 'MA0598.4', 'MA0657.2', 'MA0754.3', 'MA0804.2', 'MA0849.1',
                     'MA0876.2', 'MA1481.2', 'MA1518.3', 'MA1644.2', 'MA0039.5', 'MA0083.3', 'MA0147.4', 'MA0162.5',
                     'MA0502.3', 'MA0599.1', 'MA0659.4', 'MA0704.2', 'MA0758.1', 'MA0807.1', 'MA0851.2', 'MA0914.2',
                     'MA1489.1', 'MA1562.2', 'MA1648.2', 'MA0090.4', 'MA0147.4_all', 'MA0164.2', 'MA0505.3', 'MA0609.3',
                     'MA0663.1', 'MA0708.3', 'MA0768.3', 'MA0823.1', 'MA0852.3', 'MA1106.2', 'MA1511.2', 'MA1570.1',
                     'MA1959.2', 'MA0050.4', 'MA0095.4', 'MA0469.4', 'MA0506.3', 'MA0613.1', 'MA0664.2', 'MA0737.1',
                     'MA0770.1', 'MA0853.2', 'MA1107.3', 'MA1512.2', 'MA1607.2', 'MA1968.2', 'MA0051.2', 'MA0098.4',
                     'MA0147.4_exact', 'MA0619.2', 'MA0669.1', 'MA0740.2', 'MA0782.3', 'MA0861.2', 'MA1116.2',
                     'MA1513.2', 'MA1625.2', 'MA1990.2', 'MA0052.5', 'MA0105.4', 'MA0486.2', 'MA0526.5', 'MA0620.4',
                     'MA0670.2', 'MA0741.1', 'MA0784.3', 'MA0831.3', 'MA0863.1', 'MA1118.2', 'MA1514.2', 'MA1627.2',
                     'MA1991.2', 'MA0058.4', 'MA0108.3', 'MA0593.2', 'MA0639.2', 'MA0675.2', 'MA0742.2', 'MA0795.1',
                     'MA0868.3', 'MA1122.2', 'MA1515.2', 'MA1632.2', 'MA2325.1']


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
        return  # motif_tsv
    else:  # Remove the first column
        motif_tsv = pd.read_csv(fimo_path, sep='\t', usecols=range(1, 5), skipfooter=3,
                                engine='python', dtype=dtype_mapping)

        motif_tsv = motif_tsv.rename(columns={'sequence_name': 'chrom', 'start': 'start_m'})

        motif_tsv['chrom'] = motif_tsv["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')

        return motif_tsv


# get motif names based on the motif ids
def motif_name(motif_id, motifs_transl):
    motif_name = motifs_transl.loc[motifs_transl['matrix_id'] == motif_id, 'name']
    # query_condition = "matrix_id=='" + motif_id + "'"
    # motif_name = motifs_transl.query(query_condition)["name"]
    motif_name = ''.join(motif_name.astype(str))
    return motif_name


# get motif names and ids only for motifs that were successfully mapped by FIMO
def motif_id_name(mapped_motifs, motifs_transl, fimo_path):
    motif_ids = []
    motif_names = []

    for iid in range(len(mapped_motifs)):
        motif_id = mapped_motifs[iid]

        motif_n = motif_name(motif_id, motifs_transl)
        # motif_tsv_f = fimo_tsv(fimo_path, motif_id)
        fimo_file = os.path.join(fimo_path, motif_id, "fimo.tsv")
        if os.path.getsize(fimo_file) < 500:
            continue
        else:
            motif_ids.append(motif_id)
            motif_names.append(motif_n)

    return motif_ids, motif_names


# update motif summary for a motif
def motif_summary_update(peaks_bed, updated_df, motif_tsv, motif_name_t):
    merged = pd.merge(peaks_bed, motif_tsv, on='chrom')

    merged = merged.loc[
        (merged['start_m'] >= merged['start']) &
        (merged['stop'] <= merged['end']),
        ['peak_code', 'motif_alt_id']  # Drop unnecessary cols
    ]

    # Group by peak and count overlapping elements
    peak_element_counts = merged.groupby(['peak_code', 'motif_alt_id']).size().reset_index(name=motif_name_t)
    updated_df_n = pd.merge(updated_df, peak_element_counts[['peak_code', motif_name_t]], on='peak_code', how='left')
    updated_df_n.iloc[:, -1] = updated_df_n.iloc[:, -1].fillna(0).astype('int32')

    del merged, peak_element_counts  # Free memory
    gc.collect()
    return updated_df_n


def add_genes_to_fimo_table(fimo_df, encoding_df):
    encoding_subset = encoding_df[['peak_code', 'gene_ids']]
    resulting_df = pd.merge(fimo_df, encoding_subset, on="peak_code")
    col_names = resulting_df.columns.tolist()
    last_col = col_names.pop()  # Remove the last column
    col_names.insert(4, last_col)  # move it to the 4th index
    # Reorder the df
    resulting_df = resulting_df[col_names]
    return resulting_df


def create_genome_summary(peaks_data_path, motifs_transl_path, mapped_motifs, fimo_path, skip_rows, peak_encoding_path,
                          excel_name_path=False):
    motifs_transl = pd.read_excel(motifs_transl_path)
    peaks_bed = pd.read_csv(peaks_data_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skip_rows)
    peaks_bed['chrom'] = peaks_bed["chrom"].str.replace("Sdo_chr", "", regex=True).astype('uint16')
    peaks_bed['start'] = peaks_bed['start'] + 1
    peaks_bed['end'] = peaks_bed['end'] + 1
    peaks_bed['peak_code'] = np.arange(1, len(peaks_bed) + 1, dtype=np.uint32)
    mapped_ids, mapped_names = motif_id_name(mapped_motifs, motifs_transl, fimo_path)

    first = 1
    for motif_id, motif_name in zip(mapped_ids, mapped_names):
        motif_tsv = fimo_tsv(fimo_path, motif_id)
        # if ii == 0:
        if first:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=peaks_bed, motif_tsv=motif_tsv,
                                                 motif_name_t=motif_name)
            first = 0

        else:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=motif_summary, motif_tsv=motif_tsv,
                                                 motif_name_t=motif_name)
        del motif_tsv
        gc.collect()

    # add gene IDs to the summary df
    encoding_df = pd.read_csv(peak_encoding_path, sep='\t')
    motif_summary_w_IDs = add_genes_to_fimo_table(motif_summary, encoding_df)

    motif_summary_w_IDs.to_excel(excel_name_path, index=False)

    return f'Done!\nResults saved in:\n{excel_name_path}'


create_genome_summary(peaks_data_path=peaks_bed_path,
                      motifs_transl_path=motifs_translation_path,
                      mapped_motifs=mapped_motifs_set,
                      fimo_path=fimo_output_path,
                      excel_name_path=excel_file_path,
                      peak_encoding_path=peak_encoding_path_sys,
                      skip_rows=skipr)
