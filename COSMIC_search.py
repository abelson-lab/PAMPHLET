
import argparse

from choose_probe_target import choose_probe_placement_point
from choose_probe_target_CNV import choose_probe_placement_CNV
from collect_information import define_important_columns, \
    read_file_choose_cancer, \
    read_process_file_CNV_mutation_cbioportal, read_process_file_CNV_mutation_cosmic, read_process_file_point_mutation, \
    user_chose_options
from temp import \
    convert_dict_list, \
    define_output_file_heading, \
    find_num_gene_only_CNV, \
    integrate_CNV_point_mutation_info
from cover_entire_gene import make_probe_these_gene


def main(cosmic_mutation_file_name: str, CNV_source: str, CNV_file_name: str,
         reference_genome_filename: str, use_default: str):

    remove_intronic_mutation = True
    recurrent_definition = 2
    targeting_window_size = 80
    indel_filter_threshold = 30
    cumulative_contribution_threshold = 90
    merge_others = True
    cover_entire_gene = True
    if use_default == 'default':
        print("using default")
        use_default = True
    else:
        remove_intronic_mutation, recurrent_definition, targeting_window_size, indel_filter_threshold,\
        cumulative_contribution_threshold, merge_others, cover_entire_gene = user_chose_options()
        use_default = False

    if cover_entire_gene:
        make_probe_these_gene(reference_genome_filename)


    # # TODO this is where you add more genes, from other papers
    # l_chip_gene_set_1 = read_selected_genes(gene_list_file_name)

    important_column_heading_list, important_column_heading_list_CNV, important_column_number_list, important_column_number_list_CNV = define_important_columns()

    gene_mutation_type_info_dict_CNV = {}
    if CNV_source == 'cosmic':
        chosen_set = read_file_choose_cancer(CNV_file_name, use_default, search_CNV=True)
        read_process_file_CNV_mutation_cosmic(
            CNV_file_name,
            gene_mutation_type_info_dict_CNV,
            important_column_heading_list_CNV,
            important_column_number_list_CNV, chosen_set,
        )
    elif CNV_source == 'cbioportal':
        read_process_file_CNV_mutation_cbioportal(
            CNV_file_name,
            reference_genome_filename
        )


    choose_probe_placement_CNV(
        gene_mutation_type_info_dict_CNV,
        recurrent_definition=recurrent_definition,
        targeting_window_size=targeting_window_size,
        cumulative_contribution_threshold=cumulative_contribution_threshold,
    )

    chosen_set = read_file_choose_cancer(cosmic_mutation_file_name, use_default, search_CNV=False)

    gene_mutation_type_info_dict = {}
    read_process_file_point_mutation(cosmic_mutation_file_name,
                                     gene_mutation_type_info_dict,
                                     important_column_heading_list,
                                     important_column_number_list, chosen_set,
                                     remove_intronic_mutation)

    choose_probe_placement_point(
        gene_mutation_type_info_dict,
        recurrent_definition=recurrent_definition,
        targeting_window_size=targeting_window_size,
        indel_filter_threshold=indel_filter_threshold,
        cumulative_contribution_threshold=cumulative_contribution_threshold,
        merge_others=merge_others)



    find_num_gene_only_CNV(gene_mutation_type_info_dict,
                           gene_mutation_type_info_dict_CNV)

    all_possible_l_chip_dict = {}
    integrate_CNV_point_mutation_info(all_possible_l_chip_dict,
                                      gene_mutation_type_info_dict,
                                      gene_mutation_type_info_dict_CNV)

    # find_genes_in_cosmic_and_paper_in_paper_not_in_cosmic(
    #     all_possible_l_chip_dict, l_chip_gene_set_1)

    # turn dict into list, and sort them based on frequency
    potential_l_chip_list_headings = define_output_file_heading()

    all_possible_l_chip_list = []
    convert_dict_list(all_possible_l_chip_dict,
                      all_possible_l_chip_list)

    # process_healthy_genes_paper_genes(all_possible_l_chip_list,
    #                                   healthy_mutation_file_name,
    #                                   l_chip_gene_set_1)

    sorted_all_possible_l_chip_list = sorted(all_possible_l_chip_list,
                                             key=lambda x: float(x[1]),
                                             reverse=True)

    print('num all genes', len(sorted_all_possible_l_chip_list))

    # write into files
    # write_output(
    #     potential_l_chip_list_headings + sorted_all_possible_l_chip_list,
    #     'potential_l_chip.xlsx')


# def convert_to_tsv():
#     import pandas as pd
#
#     # names of files to read from
#     r_filenameTSV = '../../Data/Chapter01/realEstate_trans.tsv'
#
#     # names of files to write to
#     w_filenameTSV = '../../Data/Chapter01/realEstate_trans.tsv'
#
#     # read the data
#     csv_read = pd.read_csv(r_filenameTSV)
#     tsv_read = pd.read_csv(r_filenameTSV, sep='\t')
#
#     # print the first 10 records
#     print(csv_read.head(10))
#     print(tsv_read.head(10))
#
#     with open(w_filenameTSV, 'w') as write_tsv:
#         write_tsv.write(csv_read.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    # mutation file, CNV file, and gene list file


    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--mutation', help='cosmic mutation file name',
                        type=str, required=True)
    parser.add_argument('-c', '--CNV', help='cosmic or cbioportal CNV file name',
                        type=str, required=True)
    parser.add_argument('-s', '--source', help='CNV from cosmic or cbioportal',
                        type=str, required=True, choices=['cosmic', 'cbioportal'])
    parser.add_argument('-r', '--refgene', help='UCSC gene prediction tracks',
                        type=str, required=True)
    parser.add_argument('-d', '--default', help='for development purposes, do not use',
                        type=str, default='no', required=False)
    args = parser.parse_args()

    # cosmic_mutation_file_name: str, CNV_file_name: str, reference_genome_filename: str,
    #      use_default: str

    main(args.mutation, args.source, args.CNV, args.refgene, args.default)
    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
