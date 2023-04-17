import argparse

from choose_probe_target import choose_probe_placement_point
from choose_probe_target_CNV import choose_SNP_targets, divide_CNV_by_gene, plot_CNV, visualize_SNP_on_IGV
from collect_information import define_important_columns, \
    define_important_columns_CNV, read_file_choose_cancer, \
    read_process_file_CNV_mutation_cbioportal, read_process_file_CNV_mutation_cosmic, read_process_file_point_mutation, \
    user_chose_options, user_chose_options_CNV
# from temp import define_output_file_heading
from cover_entire_gene import make_probe_these_gene
from others import \
    convert_dict_list, \
    find_num_gene_only_CNV, \
    integrate_CNV_point_mutation_info, write_output_excel


def mutation_main(cosmic_mutation_filename: str, reference_genome_filename: str,
                  use_default: str, default_cancer: str):

    if use_default == 'default':
        print("using default")
        remove_intronic_mutation = True
        recurrent_definition = 2
        targeting_window_size = 80
        indel_filter_threshold = 30
        cumulative_contribution_threshold = 90
        merge_others = True
        cover_entire_gene = False
        use_default = True
    else:
        remove_intronic_mutation, recurrent_definition, targeting_window_size, indel_filter_threshold, \
        cumulative_contribution_threshold, merge_others, cover_entire_gene = user_chose_options()
        use_default = False

    if cover_entire_gene:
        make_probe_these_gene(reference_genome_filename)

    # # TODO this is where you add more genes, from other papers
    # l_chip_gene_set_1 = read_selected_genes(gene_list_file_name)

    # all processing relating to point mutations
    important_column_heading_list, important_column_number_list = define_important_columns()

    chosen_set = read_file_choose_cancer(cosmic_mutation_filename, use_default, default_cancer, search_CNV=False)

    gene_mutation_type_info_dict = {}
    read_process_file_point_mutation(cosmic_mutation_filename,
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

    # no longer used processing
    # find_num_gene_only_CNV(gene_mutation_type_info_dict,
    #                        gene_mutation_type_info_dict_CNV)
    #
    # all_possible_l_chip_dict = {}
    # integrate_CNV_point_mutation_info(all_possible_l_chip_dict,
    #                                   gene_mutation_type_info_dict,
    #                                   gene_mutation_type_info_dict_CNV)

    # find_genes_in_cosmic_and_paper_in_paper_not_in_cosmic(
    #     all_possible_l_chip_dict, l_chip_gene_set_1)

    # turn dict into list, and sort them based on frequency
    # potential_l_chip_list_headings = define_output_file_heading()

    # all_possible_l_chip_list = []
    # convert_dict_list(all_possible_l_chip_dict,
    #                   all_possible_l_chip_list)

    # process_healthy_genes_paper_genes(all_possible_l_chip_list,
    #                                   healthy_mutation_file_name,
    #                                   l_chip_gene_set_1)

    # sorted_all_possible_l_chip_list = sorted(all_possible_l_chip_list,
    #                                          key=lambda x: float(x[1]),
    #                                          reverse=True)
    #
    # print('num all genes', len(sorted_all_possible_l_chip_list))

    # write into files
    # write_output_excel(
    #     potential_l_chip_list_headings + sorted_all_possible_l_chip_list,
    #     'potential_l_chip.xlsx')


def CNV_main(CNV_source: str, CNV_filename: str, reference_genome_filename: str,
             common_snp_filename: str, use_default: str, default_cancer: str):

    top_X_CNV_gene_to_be_targeted = 10
    informative_individual_percentage = 5
    num_probe_per_individual = 100

    if use_default == 'default':
        print("using default")
        use_default = True
    else:
        informative_individual_percentage, num_probe_per_individual, \
        top_X_CNV_gene_to_be_targeted = user_chose_options_CNV()

    needed_minor_allele_frequency = informative_individual_percentage / num_probe_per_individual

    # visualize_SNP_on_IGV(common_snp_filename)

    important_column_heading_list_CNV, important_column_number_list_CNV = define_important_columns_CNV()

    gene_mutation_type_info_dict_CNV = {}

    if CNV_source == 'cosmic':
        chosen_set = read_file_choose_cancer(CNV_filename, use_default, default_cancer, search_CNV=True)
        read_process_file_CNV_mutation_cosmic(
            CNV_filename,
            gene_mutation_type_info_dict_CNV,
            important_column_heading_list_CNV,
            important_column_number_list_CNV, chosen_set,
        )
        CNV_genes = divide_CNV_by_gene(gene_mutation_type_info_dict_CNV)
    elif CNV_source == 'cbioportal':
        CNV_genes = read_process_file_CNV_mutation_cbioportal(
            CNV_filename
        )
    else:
        exit('CNV is neither from cosmic nor cbioportal')
        CNV_genes = []

    if CNV_source == 'cosmic':
        plot_CNV(CNV_genes)

    write_output_excel(CNV_genes, 'CNV gene ranked.xlsx')
    choose_SNP_targets(CNV_genes, reference_genome_filename, needed_minor_allele_frequency,
                       common_snp_filename, top_X_CNV_gene_to_be_targeted)


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

    parser.add_argument('-m', '--mutation', help='Path to cosmic mutation file',
                        type=str, required=True)
    parser.add_argument('-c', '--CNV', help='Path to cosmic or cbioportal CNV file name',
                        type=str, required=True)
    parser.add_argument('-s', '--source', help='Is CNV from cosmic or cbioportal',
                        type=str, required=True, choices=['cosmic', 'cbioportal'])
    parser.add_argument('-r', '--refgene', help='Path to UCSC gene tracks file',
                        type=str, required=True)
    parser.add_argument('-p', '--snp', help='Path to snp file name',
                        type=str, required=True)
    parser.add_argument('-d', '--default', help='For development purposes, do not use',
                        type=str, default='no', required=False)
    parser.add_argument('-a', '--defaultcancer', help='For development purposes, do not use',
                        type=str, default='m', required=False, choices=['m', 'l'])
    args = parser.parse_args()

    # cosmic_mutation_filename: str, CNV_filename: str, reference_genome_filename: str,
    #      use_default: str

    main(args.mutation, args.source, args.CNV, args.refgene, args.snp, args.default, args.defaultcancer)
    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
