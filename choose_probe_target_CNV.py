from time import sleep
from typing import Dict
import matplotlib.pyplot as plt

from choose_probe_target import parse_chromosome_position_range
from others import write_output_excel




def divide_CNV_by_gene(
        gene_mutation_type_info_dict: Dict[str, dict]):

    visualize_CNV_on_IGV(gene_mutation_type_info_dict)

    all_tumours_ids = []

    for gene, info in gene_mutation_type_info_dict.items():
        tumour_ids = info['id tumour']
        all_tumours_ids.extend(tumour_ids)

    all_unique_tumour_ids = list(set(all_tumours_ids))

    num_unique_tumour_ids = len(all_unique_tumour_ids)

    input_list_CNV = []
    # already grouped by genes
    for gene, info in gene_mutation_type_info_dict.items():

        tumour_ids = info['id tumour']
        CNV_coordinates = info['CNV coordinates']
        mutation_type = info['mutation type']

        unique_tumours = list(set(tumour_ids))
        num_unique_tumours = len(unique_tumours)

        input_list_CNV.append([
            gene, CNV_coordinates, mutation_type, unique_tumours, num_unique_tumours, num_unique_tumour_ids
        ])

    input_list_CNV.sort(key=lambda x: x[4], reverse=True)

    input_list_CNV.insert(0, ['gene', 'CNV coordinates', 'mutation type', 'unique tumours',
                              'number of unique tumours'])

    return input_list_CNV


def visualize_CNV_on_IGV(gene_mutation_type_info_dict):
    all_CNV_dict = {}
    for gene, info in gene_mutation_type_info_dict.items():

        for i in range(len(info['CNV coordinates'])):
            CNV_location = info['CNV coordinates'][i]
            all_CNV_dict[CNV_location] = ''

    all_CNV_list = []
    for CNV_location, tumour_ids in all_CNV_dict.items():

        position_list, cnv_location_end, cnv_location_start, \
        semicolon_index = parse_chromosome_position_range(CNV_location)

        chromosome = CNV_location[:semicolon_index]

        all_CNV_list.append([
            'chr' + chromosome,
            cnv_location_start,
            cnv_location_end,
            CNV_location,
            '1',
            '+',
            cnv_location_start,
            cnv_location_end,
            '21816532'
        ])

    write_output_excel(all_CNV_list, 'all CNV.xlsx')


def plot_CNV(CNV_genes):

    # find all tumours
    # skip the first row since it is the header
    all_tumour = []
    for row in CNV_genes[1:]:
        gene_tumour = row[3]
        all_tumour.extend(gene_tumour)

    all_unique_tumour = list(set(all_tumour))
    num_all_unique_tumour = len(all_unique_tumour)

    selected_tumour = []
    cover_sizes_percentage = [0]
    for row in CNV_genes[1:]:

        if cover_sizes_percentage[-1] > 0.9:
            break
        gene_tumour = row[3]
        selected_tumour.extend(gene_tumour)
        selected_tumour = list(set(selected_tumour))

        print(len(selected_tumour), num_all_unique_tumour)
        cover_sizes_percentage.append(len(selected_tumour)/num_all_unique_tumour)

    plt.plot(cover_sizes_percentage, linewidth=1)
    plt.ylabel('percentage of tumour covered')
    plt.xlabel('number of CNV target genes selected')
    plt.title('coverage of tumours with increasing ranges in myeloid CNV')
    plt.savefig('coverage.pdf')


def choose_SNP_targets(CNV_genes, reference_genome_filename, needed_minor_allele_frequency,
                       common_snp_filename):

    gene_ranges_dict = {}

    # map gene to transcription
    gene_name_transcription_dict = {}
    with open(reference_genome_filename) as hgnc_genes:
        next(hgnc_genes)
        for gene_transcript in hgnc_genes:

            gene_transcript_info_list = gene_transcript.split('\t')
            gene_name = gene_transcript_info_list[9].upper()
            gene_location = gene_name_transcription_dict.get(gene_name)

            chromosome = gene_transcript_info_list[0]
            # is a scaffold
            if len(chromosome) != 4 and len(chromosome) != 5:
                continue

            # if chromosome == 'chrY' or chromosome == 'chrX':
            #     gene_name = gene_name + ',' + chromosome

            # TODO
            # there are 10 genes with the same name, on X and Y
            # if we need both gene on both chromosome
            # change gene name to gene name + chromosome for X and Y
            # if we cannot find the gene when trying to use
            # gene_name_transcription_dict, then we add chrX or chrY to it
            # and search it again
            # and separate out X or Y before adding it to position_gene_dict

            transcription_start = gene_transcript_info_list[1]
            transcription_end = gene_transcript_info_list[2]

            gene_ranges_dict.setdefault(gene_name, []).append([transcription_start, transcription_end])

            # if gene_location is not None:
            #     earliest_transcription_start = gene_location[1]
            #     earliest_transcription_start = min(earliest_transcription_start, transcription_start)
            # else:
            #     earliest_transcription_start = transcription_start
            #
            # if gene_location is not None:
            #     latest_transcription_end = gene_location[2]
            #     latest_transcription_end = max(latest_transcription_end, transcription_end)
            # else:
            #     latest_transcription_end = transcription_end
            #
            # gene_name_transcription_dict[gene_name] = [chromosome, earliest_transcription_start, latest_transcription_end]
    print("mapped gene to transcription region")

    num_overlap = 0
    for gene, all_transcript_ranges, in gene_ranges_dict.items():
        overlap_merged = merge_overlaps(all_transcript_ranges)
        if len(overlap_merged) > 1:
            print(gene, overlap_merged)
            num_overlap += 1

    print(len(gene_ranges_dict), num_overlap)
    sleep(1000)

    # map each position to gene
    position_gene_dict = {}
    # for each CNV gene, map each position of its transcription region to that gene
    num_gene_not_in_file = 0
    for CNV_gene_section in CNV_genes[1:]:

        # TODO only do it for the top X genes (default 10)

        gene_name = CNV_gene_section[0]
        gene_transcription_region = gene_name_transcription_dict.get(gene_name)
        if gene_transcription_region is None:
            print(gene_name)
            num_gene_not_in_file += 1
            continue
        gene_chromosome = gene_transcription_region[0]
        gene_transcription_start = int(gene_transcription_region[1])
        gene_transcription_end = int(gene_transcription_region[2])

        transcription_region_list = list(range(gene_transcription_start, gene_transcription_end + 1))
        for position in transcription_region_list:
            position_gene_dict[gene_chromosome + ',' + str(position)] = gene_name


    print(len(CNV_genes) - 1, num_gene_not_in_file)
    print("mapped position to gene for all CNV genes")

    # now we have all position that correspond to an CNV gene
    # using position_gene_dict
    # map each snp to a gene (gene to snps dict)
    # i.e. collect each snp under a gene
    gene_snps_dict = {}
    with open(common_snp_filename) as snp_list:
        next(snp_list)

        for snp in snp_list:
            snp_info_list = snp.split('\t')

            snp_chromosome = snp_info_list[1]
            start_position = snp_info_list[2]
            snp_chromosome_start_position = snp_chromosome + ',' + start_position
            print(snp_chromosome_start_position)

            gene_name = position_gene_dict.get(snp_chromosome_start_position)
            if gene_name is not None:
                gene_snps_dict.setdefault(gene_name, []).append(snp)

    print("mapped CNV gene to all its snp")

    # now we have all snp of all CNV genes
    # find all snp that are good (is around MAF)
    for CNV_gene_section in CNV_genes[1:]:

        gene_name = CNV_gene_section[0]
        print(gene_name)
        gene_transcription_region = gene_name_transcription_dict.get(gene_name)

        if gene_transcription_region is None:
            print(gene_name)
        gene_chromosome = gene_transcription_region[0]
        gene_transcription_start = gene_transcription_region[1]
        gene_transcription_end = gene_transcription_region[2]

        gene_all_snps = gene_snps_dict[gene_name]

        gene_snp_list = []
        for snp in gene_all_snps:
            snp_info_list = snp.split('\t')
            num_allele = int(snp_info_list[21])
            if num_allele != 2:
                continue

            snp_chromosome = snp_info_list[1]
            if gene_chromosome == snp_chromosome:
                continue

            start_position = snp_info_list[2]
            end_position = snp_info_list[3]
            if start_position < gene_transcription_start and end_position > gene_transcription_end:
                freqs_allele = snp_info_list[24]
                minor_allele_frequency = min(freqs_allele)

                if needed_minor_allele_frequency - 1 <= float(minor_allele_frequency) <= needed_minor_allele_frequency + 1:
                    snp_location = start_position + '-' + end_position
                    gene_snp_list.append(snp_location)

        CNV_gene_section.append(gene_snp_list)

    print("kept only good snps")

    write_output_excel('CNV_genes', 'CNV_with_snps.xlsx')



def merge_overlaps(transcription_ranges):
    # Sort the unique transcript ranges by their starts
    transcription_ranges.sort()
    overlap_merged = [transcription_ranges[0]]
    # insert first transcript range into overlap_merged
    # for every other transcript
    for exon_range in transcription_ranges[1:]:

        this_exon_start = exon_range[0]
        overlap_merged_last_exon = overlap_merged[-1]
        overlap_merged_last_exon_start = overlap_merged_last_exon[0]
        overlap_merged_last_exon_end = overlap_merged_last_exon[-1]
        # if transcript ranges overlap (in other words, if this exon, starts
        # in the middle of the last transcript in overlap_merged)
        if overlap_merged_last_exon_start <= this_exon_start <= overlap_merged_last_exon_end:
            # then change the last transcript's end to whichever one is bigger (
            # last transcript end or this transcript end
            overlap_merged[-1][-1] = max(overlap_merged[-1][-1], exon_range[-1])
        # if this transcript start after the last transcript
        else:
            overlap_merged.append(exon_range)

    return overlap_merged


