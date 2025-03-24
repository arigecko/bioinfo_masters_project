import pandas as pd
import sys

if len(sys.argv) != 7:
    print("Usage: create_summary_unassigned_upd.py <peaks_bed> <motifs_translation> <fimo_output> "
          "<peak_calling_approach> <species> <output_excel_file>")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
peaks_bed_path = sys.argv[1]
motifs_translation_path = sys.argv[2]
fimo_output_path = int(sys.argv[3])
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

mapped_motifs_68 = ["MA0670.2", "MA0798.3", "MA0754.3", "MA0807.1", "MA1991.2", "MA0664.2", "MA0151.1", "MA0752.2",
                    "MA0052.5", "MA0914.2", "MA0657.2", "MA0596.1", "MA1644.2", "MA0874.2", "MA0691.1", "MA0849.1",
                    "MA0002.3", "MA0863.1", "MA0784.3", "MA0799.3", "MA0852.3", "MA0846.2", "MA0095.4", "MA0838.1",
                    "MA0613.1", "MA0058.4", "MA0831.3", "MA0090.4", "MA1625.2", "MA1990.2", "MA1968.2", "MA1632.2",
                    "MA0639.2", "MA0626.2", "MA0704.2", "MA0839.2", "MA1466.2", "MA2328.1", "MA0876.2", "MA1122.2",
                    "MA0143.5", "MA0147.4", "MA1511.2", "MA0868.3", "MA0469.4", "MA0823.1", "MA0521.3", "MA1153.2",
                    "MA1607.2", "MA0051.2", "MA0861.2", "MA0685.2", "MA0505.3", "MA0770.1", "MA0659.4", "MA0506.3",
                    "MA0037.5", "MA0098.4", "MA0162.5", "MA0063.3", "MA1627.2", "MA0083.3", "MA0502.3", "MA1513.2",
                    "MA0105.4", "MA0593.2", "MA0708.3", "MA0875.2"]


def expand_df(mapped_motifs, peaks_data_path, motifs_transl_path, skip_rows):
    peaks_bed = pd.read_csv(peaks_data_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skip_rows)
    # mapping peaks to their codes
    peak_codes = [i for i in range(1, len(peaks_bed)+1)]  # 1-based indexing
    peaks_bed.insert(len(peaks_bed.columns), 'peak_code', peak_codes)
    # peaks_bed.to_csv('peaks_bed_code_translation.tsv', sep='\t', index=False)
    motifs_transl = pd.read_excel(motifs_transl_path)

    mapped_motifs_names = []
    # expanding the dataframe
    column_names = list(peaks_bed.columns)
    new_col_dict = {}
    for ii in range(len(mapped_motifs)):
        # get the TF name based on corresponding PFM ID
        query_condition = "matrix_id=='" + mapped_motifs[ii] + "'"
        motif_name = motifs_transl.query(query_condition)["name"]
        motif_name = ''.join(motif_name.astype(str))
        new_col_dict[motif_name] = [0]*len(peaks_bed)
        column_names.append(motif_name)
        mapped_motifs_names.append(motif_name)

    new_col_dict['motif_cont'] = [0]*len(peaks_bed)
    # create new extended df
    df_to_add = pd.DataFrame(new_col_dict)

    peaks_bed_extended = pd.concat([peaks_bed, df_to_add], axis=1)

    key_list = peaks_bed_extended['peak_code'].to_list()
    # create an empty dictionary to store motifs found for each peak
    peaks_dict = {key: [] for key in key_list}
    return peaks_bed_extended, peaks_bed, mapped_motifs_names, peaks_dict


def update_motif(peaks_bed_extended, fimo_path, motif_id, motif_name, peaks_dict):

    fimo_path = fimo_path + '/' + motif_id + "/fimo.tsv"
    motif_tsv = pd.read_csv(fimo_path, sep='\t')
    # remove 3 last rows from the df - don't contain any relevant information
    motif_tsv.drop(motif_tsv.tail(3).index, inplace=True)
    amount = 0
    if motif_tsv.empty:
        return peaks_bed_extended, peaks_dict, amount

    for ind in motif_tsv.index:
        amount += 1
        scaffold = motif_tsv['sequence_name'][ind]
        start, end = motif_tsv['start'][ind] -1, motif_tsv['stop'][ind] - 1
        # filter the summary db to get the index and code of the matched peak
        filtered_df = peaks_bed_extended[(peaks_bed_extended['chrom'] == scaffold) &
                                         (peaks_bed_extended['start'] <= start) & (peaks_bed_extended['end'] >= end)]

        if len(list(filtered_df.index)) != 0:
            ind = list(filtered_df.index)[0]
            code = filtered_df.loc[ind, "peak_code"]
            # update the value
            peaks_bed_extended.loc[peaks_bed_extended.peak_code == code, motif_name] += 1
            if motif_name not in peaks_dict[code]:
                peaks_dict[code].append(motif_name)

    return peaks_bed_extended, peaks_dict, amount


def create_genome_summary(peaks_data_path, motifs_transl_path, mapped_motifs, fimo_path, skip_rows,
                          save_excel=False, excel_name_path=False):
    peaks_bed_extended, peaks_bed, mapped_motifs_names, peaks_dict = expand_df(mapped_motifs, peaks_data_path,
                                                                              motifs_transl_path, skip_rows)
    motif_total_counts = {}
    for ii in range(len(mapped_motifs)):
        motif_id = mapped_motifs[ii]
        motif_name = mapped_motifs_names[ii]
        peaks_bed_extended, peaks_dict, count = update_motif(peaks_bed_extended, fimo_path, motif_id,
                                                             motif_name, peaks_dict)
        motif_total_counts[motif_name] = count
    # calculate number of motifs in peak:
    peaks_bed_extended['present_motifs'] = (peaks_bed_extended[mapped_motifs_names] > 0).sum(axis=1)

    if save_excel:
        peaks_bed_extended.to_excel(excel_name_path, index=False)
        # save the total counts separately
        # total_counts_df = pd.DataFrame(list(data.items()), columns=['TF_name', 'total_count'])
        # total_counts_df.to_csv('TF_motifs_total_count_genome.csv', index=False)
    return peaks_bed_extended, peaks_dict


peaks_bed_extended, peaks_dict = create_genome_summary(peaks_data_path=peaks_bed_path,
                                                       motifs_transl_path=motifs_translation_path,
                                                       mapped_motifs=mapped_motifs_68,
                                                       fimo_path=fimo_output_path,
                                                       excel_name_path=excel_file_path,
                                                       save_excel=False,
                                                       skip_rows=skipr)
