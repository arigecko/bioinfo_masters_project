import pandas as pd
import numpy as np
import gc
import sys
import os

if len(sys.argv) != 7:
    print("Usage: create_summary_unassigned_upd.py <peaks_bed> <motifs_translation> <fimo_output> "
          "<peak_calling_approach> <species> <output_excel_file>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
peaks_bed_path = sys.argv[1]
motifs_translation_path = sys.argv[2]
fimo_output_path = sys.argv[3]
peak_calling_approach = sys.argv[4]
species = sys.argv[5]
excel_file_path = sys.argv[6]

if peak_calling_approach == 'CR':
    if species in ['abeo', 'abeoforma']:
        skipr = 636
    elif species in ['capsa', 'capsaspora']:
        skipr = 30
    elif species in ['sub', 'suberites', 'sponge']:
        skipr = 114
elif peak_calling_approach == 'MACS':
    skipr = 0

mapped_motifs_set = ["MA0670.2", "MA0798.3", "MA0754.3", "MA0807.1", "MA1991.2", "MA0664.2", "MA0151.1", "MA0752.2",
                     "MA0052.5", "MA0914.2", "MA0657.2", "MA0596.1", "MA1644.2", "MA0874.2", "MA0691.1", "MA0849.1",
                     "MA0002.3", "MA0863.1", "MA0784.3", "MA0799.3", "MA0852.3", "MA0846.2", "MA0095.4", "MA0838.1",
                     "MA0613.1", "MA0058.4", "MA0831.3", "MA0090.4", "MA1625.2", "MA1990.2", "MA1968.2", "MA1632.2",
                     "MA0639.2", "MA0626.2", "MA0704.2", "MA0839.2", "MA1466.2", "MA2328.1", "MA0876.2", "MA1122.2",
                     "MA0143.5", "MA0147.4", "MA1511.2", "MA0868.3", "MA0469.4", "MA0823.1", "MA0521.3", "MA1153.2",
                     "MA1607.2", "MA0051.2", "MA0861.2", "MA0685.2", "MA0505.3", "MA0770.1", "MA0659.4", "MA0506.3",
                     "MA0037.5", "MA0098.4", "MA0162.5", "MA0063.3", "MA1627.2", "MA0083.3", "MA0502.3", "MA1513.2",
                     "MA0105.4", "MA0593.2", "MA0708.3", "MA0875.2"]


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
        return #  motif_tsv
    else:  # Remove the first column
        motif_tsv = pd.read_csv(fimo_path, sep='\t', usecols=range(1, 5), skipfooter=3,
                                engine='python', dtype=dtype_mapping)

        motif_tsv = motif_tsv.rename(columns={'sequence_name': 'chrom', 'start': 'start_m'})

        motif_tsv['chrom'] = motif_tsv["chrom"].str.replace("sca", "", regex=True).astype('uint16')

        return motif_tsv


def motif_name(motif_id, motifs_transl):
    motif_name = motifs_transl.loc[motifs_transl['matrix_id'] == motif_id, 'name']
    # query_condition = "matrix_id=='" + motif_id + "'"
    # motif_name = motifs_transl.query(query_condition)["name"]
    motif_name = ''.join(motif_name.astype(str))
    return motif_name


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

    del merged, peak_element_counts # Free memory
    gc.collect()
    return updated_df_n


def create_genome_summary(peaks_data_path, motifs_transl_path, mapped_motifs, fimo_path, skip_rows,
                          excel_name_path=False):  # save_excel=False,
    motifs_transl = pd.read_excel(motifs_transl_path)
    peaks_bed = pd.read_csv(peaks_data_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skip_rows)
    peaks_bed['chrom'] = peaks_bed["chrom"].str.replace("sca", "", regex=True).astype('uint16')
    peaks_bed['peak_code'] = np.arange(1, len(peaks_bed) + 1, dtype=np.uint32)
    mapped_ids, mapped_names = motif_id_name(mapped_motifs, motifs_transl, fimo_path)
    # for ii in range(len(mapped_names)):
    #     motif_id = mapped_ids[ii]
    #     motif_name = mapped_names[ii]
    first = 1
    for motif_id, motif_name in zip(mapped_ids, mapped_names):
        motif_tsv = fimo_tsv(fimo_path, motif_id)
        # if ii == 0:
        if first:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=peaks_bed, motif_tsv=motif_tsv, motif_name_t=motif_name)
            # del peaks_bed
            first = 0

        else:
            motif_summary = motif_summary_update(peaks_bed=peaks_bed, updated_df=motif_summary, motif_tsv=motif_tsv, motif_name_t=motif_name)
        del motif_tsv
        gc.collect()
    # if save_excel:
    motif_summary.to_excel(excel_name_path, index=False)

    return f'Done!\nResults saved in:\n{excel_name_path}'


create_genome_summary(peaks_data_path=peaks_bed_path,
                      motifs_transl_path=motifs_translation_path,
                      mapped_motifs=mapped_motifs_set,
                      fimo_path=fimo_output_path,
                      excel_name_path=excel_file_path,
                      # save_excel=True,
                      skip_rows=skipr)
