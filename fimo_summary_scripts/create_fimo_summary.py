import pandas as pd
import numpy as np
import gc
import sys
import os

if len(sys.argv) != 8:
    print("Usage: create_fimo_summary.py <peaks_bed> <motifs_translation> <fimo_output> "
          "<peak_calling_approach> <species> <output_excel_file> <pwm_ids_txt>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
peaks_bed_path = sys.argv[1]
motifs_translation_path = sys.argv[2]
fimo_output_path = sys.argv[3]
peak_calling_approach = sys.argv[4]
species = sys.argv[5]
excel_file_path = sys.argv[6]
selected_motifs_path = sys.arg[7]

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
    mapped_motifs_set = [line.strip() for line in file]


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
    else:  # Remove the first column
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


# update motif summary for a motif
def motif_summary_update(peaks_bed, updated_df, motif_tsv, motif_name_t):

    merged = pd.merge(peaks_bed, motif_tsv, on='chrom')

    # filter and drop unnecessary columns
    merged = merged.loc[
        (merged['start_m'] >= merged['start']) &
        (merged['stop'] <= merged['end']),
        ['peak_code', 'motif_alt_id']
    ]

    # Group by peak and count overlapping elements
    peak_element_counts = merged.groupby(['peak_code', 'motif_alt_id']).size().reset_index(name=motif_name_t)
    updated_df_n = pd.merge(updated_df, peak_element_counts[['peak_code', motif_name_t]], on='peak_code', how='left')
    updated_df_n.iloc[:, -1] = updated_df_n.iloc[:, -1].fillna(0).astype('int32')

    # Free memory
    del merged, peak_element_counts
    gc.collect()
    return updated_df_n


def create_genome_summary(peaks_data_path, motifs_transl_path, mapped_motifs, fimo_path, skip_rows,
                          excel_name_path):
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

        if first:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=peaks_bed, motif_tsv=motif_tsv,
                                                 motif_name_t=motif_name)
            first = 0

        else:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=motif_summary, motif_tsv=motif_tsv,
                                                 motif_name_t=motif_name)
        del motif_tsv
        gc.collect()
    motif_summary.to_excel(excel_name_path, index=False)

    return f'Done!\nResults saved in:\n{excel_name_path}'


create_genome_summary(peaks_data_path=peaks_bed_path,
                      motifs_transl_path=motifs_translation_path,
                      mapped_motifs=mapped_motifs_set,
                      fimo_path=fimo_output_path,
                      excel_name_path=excel_file_path,
                      skip_rows=skipr)
