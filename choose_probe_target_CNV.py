
from typing import Dict
import matplotlib.pyplot as plt

from choose_probe_target import parse_chromosome_position_range
from others import write_output_excel, write_output_txt


def divide_CNV_by_gene(
        gene_mutation_type_info_dict: Dict[str, dict]):

    # visualize_CNV_on_IGV(gene_mutation_type_info_dict)

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

        cover_sizes_percentage.append(len(selected_tumour)/num_all_unique_tumour)

    print("number of unique tumours", num_all_unique_tumour)
    plt.plot(cover_sizes_percentage, linewidth=1)
    plt.ylabel('percentage of tumour covered')
    plt.xlabel('number of CNV target genes selected')
    plt.title('coverage of tumours with increasing ranges in myeloid CNV')
    plt.savefig('coverage.pdf')


def choose_SNP_targets(CNV_genes, reference_genome_filename, needed_minor_allele_frequency,
                       common_snp_filename, top_X_CNV_gene_to_be_targeted, targeting_window_size,
                       num_probe_per_gene_individual):

    gene_ranges_dict = {}

    # map gene to transcription
    gene_name_transcription_dict = {}
    map_gene_to_transcription_region(gene_name_transcription_dict, gene_ranges_dict, reference_genome_filename)
    print("mapped gene to transcription region")

    # num_overlap = 0
    # for gene, all_transcript_ranges, in gene_ranges_dict.items():
    #     overlap_merged = merge_overlaps(all_transcript_ranges)
    #     if len(overlap_merged) > 1:
    #         print(gene, overlap_merged)
    #         num_overlap += 1
    #
    # print(len(gene_ranges_dict), num_overlap)

    # map each position to gene
    position_gene_dict = {}
    # for each top X CNV gene, map each position of its transcription region to that gene
    map_transcribed_bases_to_gene(CNV_genes, gene_name_transcription_dict, position_gene_dict,
                                  top_X_CNV_gene_to_be_targeted)

    print("mapped position to gene for all CNV genes")


    # map each snp to a gene (gene to snps dict)
    gene_snps_dict = {}
    map_gene_to_snps(gene_snps_dict, common_snp_filename, position_gene_dict)

    snp_list_igv_format_after = []
    all_snp_per_gene = {}

    # now we have all snp of all CNV genes
    # find all snp that are good (is around MAF)
    find_all_good_snps(CNV_genes, all_snp_per_gene, gene_name_transcription_dict, gene_snps_dict,
                       needed_minor_allele_frequency, snp_list_igv_format_after, top_X_CNV_gene_to_be_targeted)

    print("kept only good snps")

    write_output_txt(snp_list_igv_format_after, 'snp_IGV_after.bed')

    write_output_excel(CNV_genes[1:top_X_CNV_gene_to_be_targeted + 1], 'CNV_with_snps.xlsx')

    # record if each base has a snp (if yes, it is in this dict)
    snp_existence_dict = {}
    for gene_name, gene_snp_start_location_list in all_snp_per_gene.items():
        for gene_snp_start_location in gene_snp_start_location_list:
            snp_existence_dict[gene_snp_start_location] = ''

    gene_snp_range_list = []
    range_placement_snp(CNV_genes, gene_name_transcription_dict, num_probe_per_gene_individual, snp_existence_dict,
                        targeting_window_size, top_X_CNV_gene_to_be_targeted, gene_snp_range_list)

    write_output_excel(gene_snp_range_list, "snp_range.xlsx")


def map_gene_to_transcription_region(gene_name_transcription_dict, gene_ranges_dict, reference_genome_filename):
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

            # there are 10 genes with the same name, on X and Y
            # if we need both gene on both chromosome
            # change gene name to gene name + chromosome for X and Y
            if chromosome == 'chrY' or chromosome == 'chrX':
                gene_name = gene_name + ',' + chromosome

            # later if we cannot find the gene when trying to use
            # gene_name_transcription_dict, then we add chrX or chrY to it
            # and search it again
            # and separate out X or Y before adding it to position_gene_dict

            transcription_start = gene_transcript_info_list[1]
            transcription_end = gene_transcript_info_list[2]

            gene_ranges_dict.setdefault(gene_name, []).append([transcription_start, transcription_end])

            if gene_location is not None:
                earliest_transcription_start = gene_location[1]
                earliest_transcription_start = min(earliest_transcription_start, transcription_start)
            else:
                earliest_transcription_start = transcription_start

            if gene_location is not None:
                latest_transcription_end = gene_location[2]
                latest_transcription_end = max(latest_transcription_end, transcription_end)
            else:
                latest_transcription_end = transcription_end

            gene_name_transcription_dict[gene_name] = [chromosome, earliest_transcription_start,
                                                       latest_transcription_end]


def map_transcribed_bases_to_gene(CNV_genes, gene_name_transcription_dict, position_gene_dict,
                                  top_X_CNV_gene_to_be_targeted):
    num_gene_not_in_file = 0
    for CNV_gene_section in CNV_genes[1:top_X_CNV_gene_to_be_targeted + 1]:

        gene_name = CNV_gene_section[0]
        gene_transcription_region = gene_name_transcription_dict.get(gene_name)
        if gene_name_transcription_dict.get(gene_name) is None:
            gene_name_X = gene_name + ',' + 'chrX'
            gene_transcription_region = gene_name_transcription_dict.get(gene_name_X)
            if gene_name_transcription_dict.get(gene_name_X) is None:
                gene_name_Y = gene_name + ',' + 'chrY'
                gene_transcription_region = gene_name_transcription_dict.get(gene_name_Y)
                if gene_name_transcription_dict.get(gene_name_Y) is None:
                    print(gene_name)
                    num_gene_not_in_file += 1
                    continue
        gene_chromosome = gene_transcription_region[0]
        gene_transcription_start = int(gene_transcription_region[1])
        gene_transcription_end = int(gene_transcription_region[2])

        transcription_region_list = list(range(gene_transcription_start, gene_transcription_end + 1))
        for position in transcription_region_list:
            position_gene_dict[str(gene_chromosome) + ',' + str(position)] = gene_name


def map_gene_to_snps(gene_snps_dict, common_snp_filename, position_gene_dict):
    """
    In the dbsnp file, we know the coordinate of each snp, but we do not know
    which gene is corresponds to. so we need to find the transcription region of
    each gene, map each gene to transcription region, then map each base in the
    transcription region to the gene. Finally, we can map gene to snp, through
    mapping snp to 'coordinate'/'a base in transcription region', a base
    in the transcription to the gene.
    :param common_snp_filename: file of common snps
    :param position_gene_dict: mapping of a base in the transcription to a gene
    :param gene_snps_dict: dict to be populated
    :return:
    """

    # now we have all position that correspond to an CNV gene
    # using position_gene_dict
    # i.e. collect each snp under a gene
    # if the starting position of the snp, matches a gene
    with open(common_snp_filename) as snp_list:
        next(snp_list)

        for snp in snp_list:
            snp_info_list = snp.split('\t')

            snp_chromosome = snp_info_list[1]
            start_position = snp_info_list[2]
            snp_chromosome_start_position = str(snp_chromosome) + ',' + str(start_position)

            gene_name = position_gene_dict.get(snp_chromosome_start_position)
            if gene_name is not None:
                gene_snps_dict.setdefault(gene_name, []).append(snp)
    print("mapped all snp gene to a CNV gene")


def find_all_good_snps(CNV_genes, all_snp_per_gene, gene_name_transcription_dict, gene_snps_dict,
                       needed_minor_allele_frequency, snp_list_igv_format_after, top_X_CNV_gene_to_be_targeted):

    for CNV_gene_section in CNV_genes[1:top_X_CNV_gene_to_be_targeted + 1]:

        # for each CNV gene, get its transcription region
        # to be later compare that the SNP is indeed within it
        # since before I only check that the starting position of the SNP
        # is a position of the transcription region
        # and also check if the chromosome match
        gene_name = CNV_gene_section[0]
        print(gene_name)
        gene_transcription_region = gene_name_transcription_dict.get(gene_name)

        if gene_transcription_region is None:
            print(gene_name)
        gene_chromosome = gene_transcription_region[0]
        gene_transcription_start = gene_transcription_region[1]
        gene_transcription_end = gene_transcription_region[2]

        # find all snps of this gene
        if gene_name in gene_snps_dict:
            gene_all_snps = gene_snps_dict[gene_name]

            gene_snp_list = []
            gene_snp_start_location_list = []
            # for one snp
            for snp in gene_all_snps:
                snp_info_list = snp.split('\t')

                # skip if it is not 2 allele
                num_allele = int(snp_info_list[21])
                if num_allele != 2:
                    continue

                # skip if it is not a single nucleotide polymorphism
                allele_class = snp_info_list[11]
                if allele_class != 'single':
                    continue

                # skip if the gene's chromosome and snp's chromosome are not the same
                # this is just an additional check since, the key to position_gene_dict
                # does involve chromosome
                snp_chromosome = snp_info_list[1]
                if gene_chromosome != snp_chromosome:
                    continue

                # check if snp is in transcription region
                start_position = snp_info_list[2]
                end_position = snp_info_list[3]
                assert start_position <= end_position
                if gene_transcription_start <= start_position and end_position <= gene_transcription_end:

                    # then check if is match our needed MAF
                    freqs_allele = snp_info_list[24]
                    minor_allele_frequency = float(min(freqs_allele.strip(',').split(',')))
                    if needed_minor_allele_frequency - 0.01 <= float(
                            minor_allele_frequency) <= needed_minor_allele_frequency + 0.01:
                        snp_location = gene_chromosome + ':' + start_position + '-' + end_position
                        gene_snp_list.append(snp_location)
                        assert int(start_position) + 1 == int(end_position)
                        gene_snp_start_location_list.append(start_position)

                        snp_list_igv_format_after.append([
                            snp_chromosome,
                            start_position,
                            end_position,
                            minor_allele_frequency,
                            '1',
                            '+',
                            start_position,
                            end_position,
                            '21816532'
                        ])

            # gather all good snp location, then append to the CNV list
            CNV_gene_section.append(gene_snp_list)
            CNV_gene_section.append(len(gene_snp_list))
            all_snp_per_gene[gene_name] = gene_snp_start_location_list


def range_placement_snp(CNV_genes, gene_name_transcription_dict, num_probe_per_gene_individual, snp_existence_dict,
                        targeting_window_size, top_X_CNV_gene_to_be_targeted, gene_snp_range_list):
    print(len(snp_existence_dict))

    # for each CNV, optimize probe placement
    for CNV_gene_section in CNV_genes[1:top_X_CNV_gene_to_be_targeted + 1]:
        gene_name = CNV_gene_section[0]
        gene_transcription_region = gene_name_transcription_dict.get(gene_name)
        if gene_name_transcription_dict.get(gene_name) is None:
            gene_name_X = gene_name + ',' + 'chrX'
            gene_transcription_region = gene_name_transcription_dict.get(gene_name_X)
            if gene_name_transcription_dict.get(gene_name_X) is None:
                gene_name_Y = gene_name + ',' + 'chrY'
                gene_transcription_region = gene_name_transcription_dict.get(gene_name_Y)
                if gene_name_transcription_dict.get(gene_name_Y) is None:
                    print(gene_name)
                    continue

        gene_transcription_start = int(gene_transcription_region[1])
        gene_transcription_end = int(gene_transcription_region[2])

        range_start = gene_transcription_start
        range_end = range_start + targeting_window_size

        # first range

        best_snp_range_list = []
        snp_range_dict = {}
        initial_range_num_snps = 0
        for i in range(range_start, range_end):
            if str(i) in snp_existence_dict:
                initial_range_num_snps += 1
        best_snp_range_list.append([range_start, initial_range_num_snps])
        snp_range_dict[range_start] = initial_range_num_snps

        print(range_start, initial_range_num_snps)

        while range_end != gene_transcription_end:
            previous_range_num_snps = snp_range_dict[range_start]
            range_start += 1
            range_end += 1

            # if the newest addition into the range, has a snp, add one
            if str(range_end) in snp_existence_dict:
                previous_range_num_snps += 1
            # if the latest subtraction has a snp, minus one
            if str(range_start - 1) in snp_existence_dict:
                previous_range_num_snps -= 1

            current_range_num_snps = previous_range_num_snps

            snp_range_dict[range_start] = current_range_num_snps
            best_snp_range_list.append([range_start, current_range_num_snps])


        # sort it
        best_snp_range_list.sort(key=lambda x: x[1])

        print(best_snp_range_list)

        num_probe_used = 0
        selected_ranges = []
        overlapping_ranges = {}
        # while it has not used all the assigned probes
        while num_probe_used < num_probe_per_gene_individual:
            # start with the best range (that cover the most snps)
            for snp_range_start_num_snps in best_snp_range_list:
                snp_range_start = snp_range_start_num_snps[0]
                num_snps = snp_range_start_num_snps[1]
                # if this does not cover any snp, continue
                if num_snps <= 0:
                    break
                # if this range overlap with another is already selected,
                if snp_range_start in overlapping_ranges:
                    continue
                # if not, add it to here
                selected_ranges.append(snp_range_start)

                # add new overlapping ranges
                for i in range(snp_range_start - 80, snp_range_start + 80):
                    overlapping_ranges[snp_range_start] = ''

            num_probe_used += 1

        gene_snp_range_list.append([gene_name, selected_ranges])


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


def visualize_SNP_on_IGV(common_snp_filename):
    snp_list_igv_format = []


    with open(common_snp_filename) as snp_list:
        next(snp_list)

        num = 0
        for snp in snp_list:
            snp_info_list = snp.split('\t')

            snp_chromosome = snp_info_list[1]
            start_position = snp_info_list[2]
            end_position = snp_info_list[3]
            freqs_allele = snp_info_list[24]

            minor_allele_frequency = float(min(freqs_allele.strip(',').split(',')))

            # 50,700,000
            # 50,900,000

            # skip if it is not 2 allele
            num_allele = int(snp_info_list[21])
            if num_allele != 2:
                continue

            # skip if it is not a single nucleotide polymorphism
            allele_class = snp_info_list[11]
            if allele_class != 'single':
                continue

            if '50700000' <= start_position and end_position <= '50900000':

                snp_list_igv_format.append([
                    snp_chromosome,
                    start_position,
                    end_position,
                    minor_allele_frequency,
                    '1',
                    '+',
                    start_position,
                    end_position,
                    '21816532'
                ])

            num += 1
            print(num, num/13347306)

    write_output_txt(snp_list_igv_format, 'snp_IGV_before.bed')
