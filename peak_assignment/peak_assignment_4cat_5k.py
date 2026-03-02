import pandas as pd
import itertools
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 9:
    print("Usage: python script.py <gtf_path> <atac_bed_file> <genome_fasta> <plot_title> "
          "<peak calling approach (CR/MACS)> <species> ")
    sys.exit(1)

# paths to the required files
# Accessing command-line arguments
isoforms_data_path_sys = sys.argv[1]
peaks_data_path_sys = sys.argv[2]
genome_f_path_sys = sys.argv[3]
summary_tsv_path_sys = sys.argv[4]
plot_title_sys = sys.argv[5]
plot_file_sys = sys.argv[6]
peak_calling_approach = sys.argv[7]
species = sys.argv[8]

if species in ['abeo', 'abeoforma']:
    scaffold_base = 'AWHISG'
elif species in ['capsa', 'capsaspora']:
    scaffold_base = 'COWC_'
elif species in ['sub', 'suberites', 'sponge', 'Sd']:
    scaffold_base = 'SUB4.g'

if peak_calling_approach == 'CR':
    if species in ['abeo', 'abeoforma']:
        skipr = 636
    elif species in ['capsa', 'capsaspora']:
        skipr = 30
    elif species in ['sub', 'suberites', 'sponge']:
        skipr = 29
# elif peak_calling_approach in ['MACS', 'other']:
else:
    skipr = 0

plot_title_sys = plot_title_sys.replace("\\n", "\n")


def map_introns_exons(isoforms_df, genome_f_path):
    # preparing the dictionary for the future mapping
    # keys are scaffolds values are lists of 0s of length equal to the length of the corresponding sequence
    genome = SeqIO.to_dict(SeqIO.parse(genome_f_path, "fasta"))
    scaffolds = list(genome.keys())

    mapping = {}
    mapping_genes = {}

    for scaffold in scaffolds:
        mapping[scaffold] = [0] * len(genome[scaffold].seq)
        mapping_genes[scaffold] = [[0]] * len(genome[scaffold].seq)

    # based on the isoforms table update the mapping in "mapping" dict
    for scaffold in scaffolds:
        scaffold_df = isoforms_df[isoforms_df['seqname'] == scaffold]
        transcripts_df = scaffold_df[scaffold_df['feature'] == "transcript"]
        transcripts = transcripts_df["transcript_id"].to_list()

        for isoform in transcripts:
            isoform = scaffold_df[(scaffold_df["transcript_id"] == isoform)]
            isoform_end = isoform.iat[0, 4]

            # store names of the genes
            gene_name = isoform.iat[0, 8]
            gene_name = int(gene_name.replace(scaffold_base, ""))

            # store all coordinates for exons and introns
            exons = []
            introns = []

            for row in range(len(isoform)):
                if row == 0:
                    continue
                else:
                    exon_start = isoform.iat[row, 3] - 1  # -1 due to the change in the indexing
                    exon_end = isoform.iat[row, 4] - 1
                    exons.append([exon_start, exon_end])
                    if exon_end + 1 != isoform_end:
                        intron_start = isoform.iat[row, 4]
                        intron_end = isoform.iat[row + 1, 3] - 2
                        introns.append([intron_start, intron_end])
            # map exons
            for exon in exons:
                for ii in range(exon[0], exon[1] + 1, 1):
                    if mapping[scaffold][ii] != 2:  # check if the position wasn't labelled as intron previously
                        mapping[scaffold][ii] = 1
                        if len(mapping_genes[scaffold][ii]) == 1 and mapping_genes[scaffold][ii][0] == 0:
                            mapping_genes[scaffold][ii] = [gene_name]
                        elif gene_name in mapping_genes[scaffold][ii]:
                            continue
                        else:
                            mapping_genes[scaffold][ii].append(gene_name)
            for intron in introns:
                for ii in range(intron[0], intron[1] + 1, 1):
                    # to check if this position was previously assigned as exon and the gene names assigned
                    # to this position need to be reset
                    exon_check = False
                    if mapping[scaffold][ii] == 1:
                        exon_check = True
                    mapping[scaffold][ii] = 2
                    if (exon_check == True) or (len(mapping_genes[scaffold][ii]) == 1 and
                                                mapping_genes[scaffold][ii][0] == 0):
                        mapping_genes[scaffold][ii] = [gene_name]
                    elif gene_name in mapping_genes[scaffold][ii]:
                        continue
                    else:
                        mapping_genes[scaffold][ii].append(gene_name)

    return mapping, mapping_genes


def prepare_peaks_bed(peaks_data_path, skipr):
    # load the file with unassigned peaks. skip 114 rows s they just list the scaffold for wic we have data
    peaks_bed = pd.read_csv(peaks_data_path, sep="\t", names=['chrom', 'start', 'end'], skiprows=skipr)
    # mapping peaks to their codes
    peak_codes = [i for i in range(1, len(peaks_bed) + 1)]  # 1-based indexing
    peaks_bed = peaks_bed.assign(peak_code=peak_codes)
    # translation is also available in the 'peaks_bed_code_translation.tsv' file
    key_list = peaks_bed['peak_code'].to_list()
    # create an empty dictionary to store motifs found for each peak
    peaks_dict_summary = {key: [] for key in key_list}
    return peaks_bed, peaks_dict_summary


def prepare_isoforms_df(isoforms_data_path):
    # Load the isoforms file and assign column names of GTF file
    isof_col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes',
                      'comments']
    isoforms_df = pd.read_csv(isoforms_data_path, sep="\t", header=None, names=isof_col_names)
    attributes = isoforms_df['attributes'].to_list()
    # prepare data for 2 new columns instead of attributes one - gene ids and transcript ids
    attributes = [x.split("; ") for x in attributes]
    gene_ids = [x[0].split(' ') for x in attributes]
    gene_ids = [x[1][1:-1] for x in gene_ids]
    transcript_id = [x[1].split(' ') for x in attributes]
    transcript_id = [x[1][1:-2] for x in transcript_id]
    # remove not-needed columns and add new ones
    isoforms_df = isoforms_df.drop(['attributes', 'comments'], axis=1)
    isoforms_df = isoforms_df.assign(gene_ids=gene_ids)
    isoforms_df = isoforms_df.assign(transcript_id=transcript_id)
    return isoforms_df


def filter_group(group):
    if group['strand'].iloc[0] == '+' or group['strand'].iloc[0] == '.':
        return group.loc[group['start'].idxmin()]  # Keep the row with smallest 'start' for forward orientation
    else:
        return group.loc[group['end'].idxmax()]  # Keep the row with largest 'end' for reverse orientation


def get_closest_dist(scaffold_plus, scaffold_minus, p_start, p_end, midp=False):
    # if peak is not a TSS peak or is not inside a transcript the closest TSS is found by calculating
    # the distance form the middle of the peak to all possible candidates.
    # Later the distances is adjusted and the half of the peak length is subtracted form it
    if midp:
        midpoint = (p_start + p_end) / 2
        peak_half = max(p_start, p_end) - midpoint

    #  orientation +
    indices_plus = list(scaffold_plus.index)
    start_plus = scaffold_plus['start'].to_list()
    #  orientation -
    indices_minus = list(scaffold_minus.index)
    start_minus = scaffold_minus['end'].to_list()

    if midp:
        start_plus = [abs(x - midpoint) for x in start_plus]
        start_minus = [abs(x - midpoint) for x in start_minus]
    else:

        start_plus = [abs(x - p_start) for x in start_plus]
        start_minus = [abs(x - p_end) for x in start_minus]

    # get the smallest distance from all retrieved coordinates to the midpoint of the peak
    if start_plus:
        min_plus_val = min(start_plus)
        min_plus_ind = indices_plus[start_plus.index(min(start_plus))]
    else:
        min_plus_val = float('inf')
    if start_minus:
        min_minus_val = min(start_minus)
        min_minus_ind = indices_minus[start_minus.index(min(start_minus))]
    else:
        min_minus_val = float('inf')

    if min(min_plus_val, min_minus_val) > 10000:
        gene_ids, distance = '-', '-'
    else:
        if midp:
            # compare the distances and output the assigned gene(s) and distance
            if min_plus_val < min_minus_val:
                gene_ids, distance = scaffold_plus.loc[min_plus_ind, "gene_ids"], int(min_plus_val)

            elif min_plus_val > min_minus_val:
                gene_ids, distance = scaffold_minus.loc[min_minus_ind, "gene_ids"], int(min_minus_val)
            elif min_plus_val == min_minus_val:
                gene_ids = [scaffold_plus.loc[min_plus_ind, "gene_ids"], scaffold_minus.loc[min_minus_ind, "gene_ids"]]
                distance = int(min_plus_val)

        else:
            if min_plus_val < min_minus_val:
                gene_ids, distance = scaffold_plus.loc[min_plus_ind, "gene_ids"], int(min_plus_val)

            elif min_plus_val > min_minus_val:
                gene_ids, distance = scaffold_minus.loc[min_minus_ind, "gene_ids"], int(min_minus_val)

            elif min_plus_val == min_minus_val:
                gene_ids = [scaffold_plus.loc[min_plus_ind, "gene_ids"], scaffold_minus.loc[min_minus_ind, "gene_ids"]]
                distance = int(min_plus_val)
    return gene_ids, distance


def intron_exon_peaks(scaffold_df, scaffold, p_start, p_end, mapped_genome, mapped_genes):
    peak_mapping = mapped_genome[scaffold][p_start:p_end + 1]
    mapping_set = set(peak_mapping)
    peak_mapping_genes = mapped_genes[scaffold][p_start:p_end + 1]
    # assign category to the peak:
    category = ""  # category to which the peak will be assigned to
    transcript_intron = True
    if 2 in mapping_set:

        category = "intron"
        if 0 in mapping_set:
            transcript_intron = False
    elif 1 in mapping_set:
        if (len(mapping_set) == 2) and (0 in mapping_set):
            category = "exon_intergenic"
        elif len(mapping_set) == 1:
            category = "pure_exon"
    elif (len(mapping_set) == 1) and (0 in mapping_set):
        category = "intergenic"

    if category == 'intron' and transcript_intron:
        # peak overlaps with an intron. can be fully in the transcript or overlap with and intergenic region too

        # get indices of all intron positions to determine the genes we need to consider
        intron_indices = [i for i in range(len(peak_mapping)) if peak_mapping[i] == 2]
        intron_genes = [peak_mapping_genes[i] for i in intron_indices]

        # flatten the nested list and also de-code the gene names back into the full format
        intron_genes = list(set(list(itertools.chain(*intron_genes))))
        intron_genes = [scaffold_base + str(i) for i in intron_genes]
        # get information on which genes have which orientation
        intron_isoforms_pl = scaffold_df[(scaffold_df['gene_ids'].isin(intron_genes)) &
                                         (scaffold_df['feature'] == 'transcript') & (scaffold_df['strand'] == "+")]
        intron_isoforms_mn = scaffold_df[(scaffold_df['gene_ids'].isin(intron_genes)) &
                                         (scaffold_df['feature'] == 'transcript') & (scaffold_df['strand'] == "-")]

        # get the distance and the assigned genes
        gene_ids, distance = get_closest_dist(intron_isoforms_pl, intron_isoforms_mn, p_start, p_end, midp=False)
        if gene_ids == "-":
            return category, gene_ids, distance
        else:
            return category, intron_genes, distance

    # all other cases - peak is overlapping with the end of the transcript; peak is not in a transcript;
    else:
        # get information on which genes have
        other_isoforms_pl = scaffold_df[(scaffold_df['feature'] == 'transcript') & (scaffold_df['strand'] == "+")]
        other_isoforms_mn = scaffold_df[(scaffold_df['feature'] == 'transcript') & (scaffold_df['strand'] == "-")]

        gene_ids, distance = get_closest_dist(other_isoforms_pl, other_isoforms_mn, p_start, p_end, midp=True)
        if (category == 'intergenic' and distance == '-') or (category == 'intergenic' and distance > 5000):
            category = 'distal intergenic'
        elif category == 'intergenic' and distance < 5000:
            category = 'proximal intergenic'

        return category, gene_ids, distance


def save_results(peaks_dict_summary, summary_tsv_path):
    # create a df from the summary dictionary
    peak_codes = list(peaks_dict_summary.keys())
    genes_list = []
    dist_l = []
    category = []

    for peak in peak_codes:
        genes_list.append(peaks_dict_summary[peak][0])
        dist_l.append(peaks_dict_summary[peak][1])
        category.append(peaks_dict_summary[peak][2])

    genes_list_joined = []
    for ii in range(len(genes_list)):
        if isinstance(genes_list[ii], list):
            genes_list_joined.append("; ".join(genes_list[ii]))
        else:
            genes_list_joined.append(genes_list[ii])

    peaks_assignment_dict = {'peak_code': peak_codes, 'gene_ids': genes_list_joined,
                             'dist': dist_l, 'category': category}
    peaks_assignment_df = pd.DataFrame(peaks_assignment_dict)

    peaks_assignment_df['category'] = peaks_assignment_df['category'].replace(
        ['intron', 'pure_exon', 'exon_intergenic'],
        'intragenic')

    peaks_assignment_df.to_csv(summary_tsv_path, sep='\t', index=False)


def assign_peaks(isoforms_data_path, peaks_data_path, genome_f_path, summary_tsv_path, skiprows):
    # load and modify data
    isoforms_df = prepare_isoforms_df(isoforms_data_path)
    peaks_bed, peaks_dict_summary = prepare_peaks_bed(peaks_data_path, skiprows)
    mapped_genome, mapped_genes = map_introns_exons(isoforms_df, genome_f_path)

    # go through all peaks and check which genes can bbe assigned to it
    for ind in peaks_bed.index:
        # store the data for a peak
        scaffold = peaks_bed['chrom'][ind]
        code = peaks_bed['peak_code'][ind]
        p_start = peaks_bed['start'][ind]
        p_end = peaks_bed['end'][ind]

        scaffold_df = isoforms_df[(isoforms_df['seqname'] == scaffold)]

        # first check if the peak overlaps with TSS of te most 3' isoform.
        # (here TSS is the first nucleotide at the  beginning of the gene)
        scaffold_transcripts = scaffold_df[scaffold_df['feature'] == 'transcript']
        filtered_transcripts_df = scaffold_transcripts.groupby('gene_ids', group_keys=False).apply(
            filter_group)  # , include_groups=False)

        # check for both gene orientations
        dist_0_plus = filtered_transcripts_df[(filtered_transcripts_df['strand'] == "+") &
                                              (filtered_transcripts_df['start'] >= p_start) &
                                              (filtered_transcripts_df['start'] < p_end)]

        dist_0_minus = filtered_transcripts_df[(filtered_transcripts_df['strand'] == "-") &
                                               (filtered_transcripts_df['end'] > p_start) &
                                               (filtered_transcripts_df['end'] <= p_end)]

        if (len(dist_0_minus) > 0) or (len(dist_0_plus) > 0):

            dist_0_minus_ids = dist_0_minus["gene_ids"].to_list()
            dist_0_plus_ids = dist_0_plus["gene_ids"].to_list()
            genes = list(set(dist_0_minus_ids + dist_0_plus_ids))
            category = "TSS"
            distance = 0

            peaks_dict_summary[code] = [genes, distance, category]

        # check for peaks not overlapping with TSS
        else:
            category, genes, distance = intron_exon_peaks(scaffold_df, scaffold, p_start, p_end,
                                                          mapped_genome, mapped_genes)
            peaks_dict_summary[code] = [genes, distance, category]

    save_results(peaks_dict_summary, summary_tsv_path)


# run the de-novo peak assignment
assign_peaks(isoforms_data_path=isoforms_data_path_sys, peaks_data_path=peaks_data_path_sys,
             genome_f_path=genome_f_path_sys, summary_tsv_path=summary_tsv_path_sys, skiprows=skipr)

assigned_peaks_df = pd.read_csv(summary_tsv_path_sys, sep='\t')
total_len = len(assigned_peaks_df)


# Creating dataset
all_categories = assigned_peaks_df['category'].to_list()
# Create a DataFrame from the list
all_categories_df = pd.DataFrame(all_categories, columns=['Category'])

categories = ['TSS', 'intragenic', 'proximal intergenic', 'distal intergenic']
# colors for the bars in the plot
colors = ['#0464a4', '#a41c4c', '#94c454', '#28894f']


# Ensure consistent order by setting a Categorical dtype
all_categories_df['Category'] = pd.Categorical(all_categories_df['Category'], categories=categories, ordered=True)

# Calculate the frequency of each category
category_counts = all_categories_df['Category'].value_counts().reindex(categories).reset_index()
category_counts.columns = ['Category', 'Frequency']

category_counts['Frequency'] = category_counts['Frequency'].fillna(0)

# Calculate the total number of occurrences
total_counts = category_counts['Frequency'].sum()
# Calculate the percentage of each category
category_counts['Percentage'] = (category_counts['Frequency'] / total_counts) * 100

palette = dict(zip(categories, colors))

sns.set_style("white")
# Create the barplot
fig = plt.figure(figsize=(6, 4))
barplot = sns.barplot(x='Category', y='Percentage', data=category_counts, palette=palette)
# setting y-axis limits and ticks
plt.ylim(0, 100)
plt.yticks(range(0, 101, 10))
# Add title and labels
plt.title(plot_title_sys, fontsize=15, pad=20)
# adjusting tick font size
plt.xticks(fontsize=16, rotation=20)
plt.yticks(fontsize=15)
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}%'))
plt.xlabel('Category', fontsize=12)
plt.ylabel('Percentage (%)', fontsize=12)
# Add percentages on top of the bars
for p in barplot.patches:
    barplot.annotate(format(p.get_height(), '.2f') + '%',
                     (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center', xytext=(0, 4),
                     textcoords='offset points')

# Save the plot
fig.savefig(plot_file_sys, dpi=fig.dpi, bbox_inches='tight')
