from typing import Dict

from choose_probe_target import parse_chromosome_position_range
from temp import write_output


def choose_probe_placement_CNV(
        gene_mutation_type_info_dict: Dict[str, dict],
        targeting_window_size: int = 80,
        cumulative_contribution_threshold: int = 90,
        recurrent_definition: int = 2, ):
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

    write_output(input_list_CNV, 'input_gene_CNV_cosmic.xlsx')


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

    write_output(all_CNV_list, 'all CNV.xlsx')
