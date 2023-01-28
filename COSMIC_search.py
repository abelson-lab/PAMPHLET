import csv
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt


def read_file_choose_cancer(cosmic_mutation_file_name):
    all_primary_tissue_set = []
    all_primary_histology_set = []
    all_histology_subtype_one_set = []

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')

        for row in csv_reader:
            all_primary_tissue_set = row[7]

    all_primary_tissue_set = list(set(all_primary_tissue_set))
    print_set_as_numbered_list(all_primary_tissue_set)
    cancer_list_indices_list = read_number_as_list_indices()
    chosen_primary_tissue_set = []
    report_element_at_indices(all_primary_tissue_set, cancer_list_indices_list,
                              chosen_primary_tissue_set)

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')
        for row in csv_reader:
            if row[7] in chosen_primary_tissue_set:
                all_primary_histology_set.append(row[11])

    all_primary_histology_set = list(set(all_primary_histology_set))
    print_set_as_numbered_list(all_primary_histology_set)
    cancer_list_indices_list = read_number_as_list_indices()
    chosen_primary_histology_set = []
    report_element_at_indices(all_primary_histology_set,
                              cancer_list_indices_list,
                              chosen_primary_histology_set)

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')
        for row in csv_reader:
            if row[7] in chosen_primary_tissue_set and row[
                11] in chosen_primary_histology_set:
                all_histology_subtype_one_set.append(row[12])

    all_histology_subtype_one_set = list(set(all_histology_subtype_one_set))
    print_set_as_numbered_list(all_histology_subtype_one_set)
    cancer_list_indices_list = read_number_as_list_indices()
    chosen_histology_subtype_one_set = []
    report_element_at_indices(all_histology_subtype_one_set,
                              cancer_list_indices_list,
                              chosen_histology_subtype_one_set)

    chosen_sets = [chosen_primary_tissue_set, chosen_primary_histology_set,
                   chosen_histology_subtype_one_set]
    return chosen_sets


def print_set_as_numbered_list(a_list):
    i = 1
    for cancer in a_list:
        print(i, '. ', cancer)
        i += 1


def read_number_as_list_indices():
    cancer_list_indices = input(
        "choose primary tissue/primary histology/histology subtype 1"
        " you want to target from the list below, "
        "separate by comma")

    cancer_list_indices_list = cancer_list_indices.split(',')

    return cancer_list_indices_list


def report_element_at_indices(cancer_list, cancer_list_indicies_list,
                              chosen_cancer_list):
    for indices in cancer_list_indicies_list:
        chosen_cancer_list.append(cancer_list[indices])


def user_chose_options():
    remove_intronic_mutation = input(
        "Do you want include intronic mutations, default is no")
    recurrent_definition = input("How do you define recurrent mutation,"
                                 " default is a mutation is only consider"
                                 " recurrent if it appear in more than 2"
                                 " tumours, including 2 tumours")
    targeting_window_size = input("What is your window size, or the number"
                                  " of base pairs that can be effectively"
                                  " targeted, default is 80bp")
    indel_filter_threshold = input("What size of indel is too big to be cover"
                                   " by your window size, default is 30bp or"
                                   "greater will be filtered out")
    cumulative_contribution_threshold = input("what is the coverage of all"
                                              " recurrent tumours you want for"
                                              " the probes, default is reporting"
                                              " until it covers 80% of the"
                                              " probes")
    if remove_intronic_mutation == "yes":
        remove_intronic_mutation = True
    else:
        remove_intronic_mutation = False

    recurrent_definition = int(recurrent_definition)
    targeting_window_size = int(targeting_window_size)
    indel_filter_threshold = int(indel_filter_threshold)
    cumulative_contribution_threshold = int(cumulative_contribution_threshold)

    return remove_intronic_mutation, recurrent_definition, targeting_window_size, indel_filter_threshold, cumulative_contribution_threshold


def main(cosmic_mutation_file_name, cosmic_CNV_file_name,
         gene_list_file, healthy_mutation_file_name):
    remove_intronic_mutation, recurrent_definition, targeting_window_size, indel_filter_threshold, cumulative_contribution_threshold = user_chose_options()

    chosen_set = read_file_choose_cancer(cosmic_mutation_file_name)

    # TODO this is where you add more genes, from other papers
    gene_list = []
    l_chip_gene_set_1 = read_selected_genes(gene_list, gene_list_file)

    important_column_heading_list, important_column_heading_list_CNV, important_column_number_list, important_column_number_list_CNV = define_important_columns()

    gene_cell_mutation_type_info_dict = {}
    read_process_file_point_mutation(cosmic_mutation_file_name,
                                     gene_cell_mutation_type_info_dict,
                                     important_column_heading_list,
                                     important_column_number_list,
                                     remove_intronic_mutation,
                                     chosen_set)

    categorize_then_choose_probe_placement_within(
        gene_cell_mutation_type_info_dict,
        recurrent_definition=recurrent_definition,
        targeting_window_size=targeting_window_size,
        indel_filter_threshold=indel_filter_threshold,
        cumulative_contribution_threshold=cumulative_contribution_threshold)

    gene_cell_mutation_type_info_dict_CNV = {}
    read_process_file_CNV_mutation(cosmic_CNV_file_name,
                                   gene_cell_mutation_type_info_dict_CNV,
                                   important_column_heading_list_CNV,
                                   important_column_number_list_CNV)

    find_num_gene_only_CNV(gene_cell_mutation_type_info_dict,
                           gene_cell_mutation_type_info_dict_CNV)

    all_possible_l_chip_dict = {}
    integrate_CNV_point_mutation_info(all_possible_l_chip_dict,
                                      gene_cell_mutation_type_info_dict,
                                      gene_cell_mutation_type_info_dict_CNV)

    find_genes_in_cosmic_and_paper_in_paper_not_in_cosmic(
        all_possible_l_chip_dict, l_chip_gene_set_1)

    # turn dict into list, and sort them based on frequency
    potential_l_chip_list_headings = define_output_file_heading()

    all_possible_l_chip_list = []
    convert_dict_list(all_possible_l_chip_dict,
                      all_possible_l_chip_list)

    process_healthy_genes_paper_genes(all_possible_l_chip_list,
                                      healthy_mutation_file_name,
                                      l_chip_gene_set_1)

    sorted_all_possible_l_chip_list = sorted(all_possible_l_chip_list,
                                             key=lambda x: float(x[1]),
                                             reverse=True)

    print('num all genes', len(sorted_all_possible_l_chip_list))

    # write into files
    write_output(
        potential_l_chip_list_headings + sorted_all_possible_l_chip_list,
        'potential_l_chip.xlsx')


def define_output_file_heading():
    potential_l_chip_list_headings = [[
        'gene',
        'T Cell Frequency non-CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'T Cell Frequency CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'B Cell Frequency non-CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'B Cell Frequency CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'Other Cell Frequency non-CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'Other Cell Frequency CNV', 'tumour num', 'total tumour num',
        'only non-CNV', 'only CNV', 'both', 'position',
        'CNV or not', 'found in healthy', 'in paper 1', 'in paper 2'
    ]]
    return potential_l_chip_list_headings


def process_healthy_genes_paper_genes(all_possible_l_chip_list,
                                      healthy_mutation_file_name,
                                      l_chip_gene_set_1):
    genes_in_healthy = which_genes_were_mutated_in_healthy(l_chip_gene_set_1,
                                                           healthy_mutation_file_name)
    unhealthy_genes = []
    for i in range(len(all_possible_l_chip_list)):
        gene_info = all_possible_l_chip_list[i]
        gene = gene_info[0]
        if gene in genes_in_healthy:
            all_possible_l_chip_list[i].append('healthy')
        else:
            all_possible_l_chip_list[i].append('___')
            unhealthy_genes.append(gene)
    print('num genes found in healthy individuals', len(genes_in_healthy))
    print('genes found mutated in healthy individuals', genes_in_healthy)
    print('num genes not found in healthy individuals', len(unhealthy_genes))
    for i in range(len(all_possible_l_chip_list)):
        gene_info = all_possible_l_chip_list[i]
        gene = gene_info[0]
        if gene in l_chip_gene_set_1:
            all_possible_l_chip_list[i].append('YES')
        else:
            all_possible_l_chip_list[i].append('___')


def find_genes_in_cosmic_and_paper_in_paper_not_in_cosmic(
        all_possible_l_chip_dict, l_chip_gene_set_1):
    """this section is only for calculating remaining genes"""
    # separate gene_info_dict into 2 outputs files
    # known: in the gene set
    # potential: not in the gene set
    known_l_chip_dict = {}
    potential_l_chip_dict = {}
    for gene_histology, info in all_possible_l_chip_dict.items():
        gene_histology_list = gene_histology.split(';')
        gene = gene_histology_list[0]
        if gene in l_chip_gene_set_1:
            known_l_chip_dict[gene_histology] = info
        else:
            potential_l_chip_dict[gene_histology] = info

    # in the initial 235 ones
    print('known l chip gene cell-type pair', len(known_l_chip_dict))
    # not in the 235
    print('potential l chip gene cell-type pair', len(potential_l_chip_dict))

    remaining_genes = []
    known_l_chip_gene_list = []
    for gene_histology in list(known_l_chip_dict.keys()):
        gene_histology_list = gene_histology.split(';')
        gene = gene_histology_list[0]
        known_l_chip_gene_list.append(gene)
    known_l_chip_dict_set = set(known_l_chip_gene_list)
    for gene in l_chip_gene_set_1:
        if gene not in known_l_chip_dict_set:
            remaining_genes.append(gene)

    # the genes in 235, but not found on cosmic
    print('remaining', len(remaining_genes))
    write_output(remaining_genes, 'remaining_genes.xlsx')


def integrate_CNV_point_mutation_info(all_possible_l_chip_dict,
                                      gene_cell_mutation_type_info_dict,
                                      gene_cell_mutation_type_info_dict_CNV):
    # add CNV info into it
    # TODO can there be CNV that have more than 1 study but still
    #   does not have more than 1 in point mutation?
    for gene_cell_type, info in gene_cell_mutation_type_info_dict.items():
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        cell_type = gene_cell_type_list[1]
        all_possible_l_chip_dict[gene_cell_type] = info
        CNV_key = gene + ';' + cell_type + ';' + 'CNV'
        # if this gene has a CNV
        if CNV_key in gene_cell_mutation_type_info_dict_CNV:
            # take first two part of this gene (gene name + cell type)
            # add 'CNV' to it, to get the CNV version of it
            this_gene_info = gene_cell_mutation_type_info_dict_CNV[
                CNV_key]
            all_possible_l_chip_dict[CNV_key] = this_gene_info


def find_num_gene_only_CNV(gene_cell_mutation_type_info_dict,
                           gene_cell_mutation_type_info_dict_CNV):
    all_mutated_gene = []
    all_CNV_gene = []
    for gene_cell_type in gene_cell_mutation_type_info_dict:
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        all_mutated_gene.append(gene)
    for gene_cell_type in gene_cell_mutation_type_info_dict_CNV:
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        all_CNV_gene.append(gene)
    all_mutated_gene_set = set(all_mutated_gene)
    all_CNV_gene_set = set(all_CNV_gene)
    only_CNV = all_CNV_gene_set.difference(all_mutated_gene_set)
    print('num gene only CNV', len(only_CNV))


def read_process_file_CNV_mutation(cosmic_CNV_file_name,
                                   gene_cell_mutation_type_info_dict_CNV,
                                   important_column_heading_list_CNV,
                                   important_column_number_list_CNV):
    # each gene cell-type refers to the gene ABC in t-cell CNV

    cell_type_num_tumor_dict_CNV = {}
    read_CNV_file(cosmic_CNV_file_name,
                  cell_type_num_tumor_dict_CNV,
                  gene_cell_mutation_type_info_dict_CNV,
                  important_column_heading_list_CNV,
                  important_column_number_list_CNV)
    print('gene cell-type pair CNV', len(gene_cell_mutation_type_info_dict_CNV))
    print('CNV cell type distribution', cell_type_num_tumor_dict_CNV)
    # and store in gene_cell_mutation_type_info_dict_CNV
    calculate_CNV_frequency(cell_type_num_tumor_dict_CNV,
                            gene_cell_mutation_type_info_dict_CNV)


def read_process_file_point_mutation(cosmic_mutation_file_name,
                                     gene_cell_mutation_type_info_dict,
                                     important_column_heading_list,
                                     important_column_number_list,
                                     remove_intronic_mutation,
                                     chosen_set):
    # read mutation, group by gene+cell+mutation type (in gene_info_dict)
    # and count how many sample each cell type has (in cell_type_num_tumour_dict)
    cell_type_num_tumour_dict = {}
    read_mutation_file(cosmic_mutation_file_name,
                       cell_type_num_tumour_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_cell_mutation_type_info_dict,
                       remove_intronic_mutation,
                       chosen_set)
    print('gene cell-type pair non-CNV', len(gene_cell_mutation_type_info_dict))
    print('tumour distribution for each point mutation cell type ',
          cell_type_num_tumour_dict)
    # and store in gene_info_dict
    calculate_mutation_frequency(cell_type_num_tumour_dict,
                                 gene_cell_mutation_type_info_dict)


def define_important_columns():
    # these two should match (i.e. the 6th column should be id tumour)
    important_column_heading_list = ['id tumour', 'histology subtype 2',
                                     'histology subtype 3',
                                     'mutation CDS', 'mutation description',
                                     'GRch', 'mutation genome position',
                                     'study id']
    important_column_number_list = [6, 13, 14, 19, 21, 24, 25, 30]
    important_column_heading_list_CNV = ['id tumour', 'histology subtype 2',
                                         'histology subtype 3',
                                         'study id',
                                         'GRch', 'mutation genome position']
    important_column_number_list_CNV = [4, 11, 12, 17, 18, 19]
    return important_column_heading_list, important_column_heading_list_CNV, important_column_number_list, important_column_number_list_CNV


def read_selected_genes(gene_list, gene_list_file):
    with open(gene_list_file) as gene_list_file:
        lines = gene_list_file.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
        gene_list.extend(stripped_lines)
    l_chip_gene_set_1 = set(gene_list)
    return l_chip_gene_set_1


def read_mutation_file(cosmic_mutation_file_name,
                       cell_type_num_total_tumour_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_cell_mutation_type_dict,
                       remove_intronic_mutation,
                       chosen_set):
    chosen_primary_tissue_set = chosen_set[0]
    chosen_primary_histology_set = chosen_set[1]
    chosen_histology_subtype_one_set = chosen_set[2]

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')

        unique_tumour_before_intronic_mutation = []
        unique_mutation_before_intronic_mutation = []
        for row in csv_reader:
            if row[7] in chosen_primary_tissue_set and row[
                11] in chosen_primary_histology_set and row[
                12] in chosen_histology_subtype_one_set:

                unique_tumour_before_intronic_mutation.append(row[4])
                unique_mutation_before_intronic_mutation.append(row[25])

        unique_tumour_before_intronic_mutation = set(
            list(unique_tumour_before_intronic_mutation))
        unique_mutation_before_intronic_mutation = set(
            list(unique_mutation_before_intronic_mutation))

        # position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(row[25])

        print('total number of mutations before intronic filter',
              len(unique_mutation_before_intronic_mutation))
        print('total number of tumours with mutations before intronic filter',
              len(unique_tumour_before_intronic_mutation))

        for row in csv_reader:
            # 0, 4, 11, 12, 13, 14.  19, 21, 24, 25, 30
            # gene name, tumour name, primary histology and its subtypes (3),
            # mutation CDS, mutation desc, GRch, mutation genome position, study id,
            # 7 is primary site

            # if primary histology is lymphoid_neoplasm
            # if row[11] == "haematopoietic_neoplasm":
            # if row[11] == 'lymphoid_neoplasm':
            if row[7] in chosen_primary_tissue_set and row[
                11] in chosen_primary_histology_set and row[
                12] in chosen_histology_subtype_one_set:
                # group mutations by gene name + histology subtype 1
                gene_ensembl = row[0].split('_')
                gene_name = gene_ensembl[0]

                mutation_CDS = row[19]
                # assert 'c.' in mutation_CDS
                # check_variant_type(row)

                # if it is intron (where the position is known) and not splice site and not
                # ending or start in an exon, and both start and end in two different introns,
                # then skip it
                if remove_intronic_mutation:
                    filter_out_flag = filter_known_pure_intronic_mutation(
                        mutation_CDS)
                    if filter_out_flag:
                        continue
                    # we don't need this, there are no such mutation in this file
                    # filter_out_flag_intron_range = is_intron_range_mutation(
                    #     mutation_CDS)
                    #
                    # if filter_out_flag_intron_range:
                    #     check_intronic_mutation_diff_intron(mutation_CDS)

                if row[12] == 'NS':
                    continue

                if row[21] == 'Substitution - coding silent':
                    continue

                histology = row[12]
                # if 'myelo' in histology:
                #     dict_key_name = gene_name + ';' + 'myelo' + ';' + 'point'
                #     cell_type = 'myelo'
                # else:
                #     continue

                # TODO remove this, since now if we need to only target t-cell
                #   user would input when asked, what histology subtype 1
                #   dict_key_name should just be
                # dict_key_name = gene_name + ';' + 'point'
                # cell_type can be literally anything
                if 'T_cell' in histology or 'anaplastic' in histology \
                        or 'lymphomatoid_papulosis' in histology \
                        or 'post_transplant_lymphoproliferative_disorder' in histology \
                        or 'mycosis_fungoides-Sezary_syndrome' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'point'
                    cell_type = 'T_cell'
                elif 'B_cell' in histology or 'Burkitt' in histology or 'chronic' in histology \
                        or 'follicular_lymphoma' in histology or 'hairy' in histology \
                        or 'hodgkin' in histology or 'lymphoplasmacytic_lymphoma' in histology \
                        or 'mantle' in histology or 'marginal' in histology \
                        or 'plasma_cell_myeloma' in histology or 'plasmacytoma' in histology \
                        or 'effusion' in histology or 'MALT' in histology:
                    dict_key_name = gene_name + ';' + 'B_cell' + ';' + 'point'
                    cell_type = 'B_cell'
                else:
                    dict_key_name = gene_name + ';' + 'Other' + ';' + 'point'
                    cell_type = 'Other'

                specific_gene_cell_type_dict = gene_cell_mutation_type_dict.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list)):
                    column_heading = important_column_heading_list[i]
                    column_num = important_column_number_list[i]
                    specific_gene_cell_type_dict.setdefault(
                        column_heading, []).append(row[column_num])
                    # for example
                    # specific_gene_histology_dict.setdefault('sample name', []).append(row[4])

                cell_type_num_total_tumour_dict.setdefault(cell_type,
                                                           []).append(row[6])

            # counter += 1
            # print(str(counter*100/3544360) + '%')

    for cell_type, tumour_list in cell_type_num_total_tumour_dict.items():
        cell_type_num_total_tumour_dict[cell_type] = len(set(tumour_list))


def check_variant_type(row):
    """
    if all variant were covered by this check then nothing should be printed
    :param row:
    :return:
    """

    # every variants' nomenclature included in this file
    # c.   optional (*-)   any digit one or more times, optional (+ or - any digit plus sign)
    # optional (_ optional(* or -) any digit one or more optional (+ or - any digit plus sign))
    # one of A C G T d i ? >

    # c.? is an unknown mutation
    # c.1_*del is a whole gene deletion (it does mean 1 to stop codon is deleted)

    # there was a couple of them with this format c.*-9A>G
    # i.e. the start of the range being 9 base pair before the stop codon
    # this would be c.* optional (-) any digit one or more and then ACGTdi
    # c.*-6506_*-6505dup

    # there are two variants( where the break point not sequence)
    # c.(5830_5832)insC
    # c.(5830_5832)insC

    # the position in intron, but unknown exactly
    # c.47-?_780+?del

    # break point not sequenced and unknown position in intron
    # c.(801 +?_802 -?)del

    # this case exist, I did not check for it
    # c.493_*del

    if re.search(
            "c\\.[*\\-]?[0-9]+([+\\-][0-9]+)?(_[*\\-]?[0-9]+([+\\-][0-9]+)?)?[ACGTdi>?]",
            row[19]) is None and \
            re.search("c\\.[?]", row[19]) is None and re.search("c\\.1_\\*del",
                                                                row[19]) is None \
            and re.search("c\\.\\*[-]?[0-9]+(_\\*[-]?[0-9]+)?[ACGTdi>?]",
                          row[19]) is None \
            and re.search("c\\.\\([0-9]+_[0-9]+\\)[ACGTdi>?]", row[19]) is None \
            and re.search(
        "c\\.(\\()?[0-9]+[+\\-]\\?(_[0-9]+[+\\-]\\?)?(\\))?[ACGTdi>?]",
        row[19]) is None:
        print(row[19])


def is_intron_range_mutation(mutation_CDS):
    filter_out_flag_intron_range = re.search(
        'c\\.[0-9]+([+\\-][0-9]+)_[0-9]+([+\\-][0-9]+)[ACGTdi>?]',
        mutation_CDS) is not None
    return filter_out_flag_intron_range


def filter_known_pure_intronic_mutation(mutation_CDS):
    # check for splice site (defined as 2 base pair before and after
    #   the nearest nucleotide in the nearest exon)
    # any digit (once)
    # plus sign or minus sign (once)
    # either any digit between 0-2
    # followed by any alphabet or underline character _ or equal
    # character, and bracket
    # it cannot just be [+\\-][0-2] because -2 is a nucleotide upstream of start codon
    # and because 3241+200 will match [0-9][+\\-][0-2]

    # all intron variants in this file, will match
    # either [0-9][+\\-][0-9]

    # not include this, since we don't know if this is a splice site
    # or [0-9][+\\-]\\?

    # ,also if it is a range, and it starts in intron and ends in exon or vice versa

    # is not = does not match = return None
    # is = does match = return not None

    filter_out_flag = (re.search('[0-9][+\\-][0-9]', mutation_CDS) is not None) \
                      and re.search('[0-9][+\\-][0-2][a-zA-Z_=)]',
                                    mutation_CDS) is None \
                      and re.search(
        'c\\.[*\\-]?[0-9]+([+\\-][0-9]+)_[*\\-]?[0-9]+[ACGTdi>?]',
        mutation_CDS) is None \
                      and re.search(
        'c\\.[*\\-]?[0-9]+_[*\\-]?[0-9]+([+\\-][0-9]+)[ACGTdi>?]',
        mutation_CDS) is None
    return filter_out_flag


def check_intronic_mutation_diff_intron(mutation_CDS):
    # this check if (for example c.4072-1234_5155-246del)
    # the deletion is from an intron between exon 29 and exon 30
    # all the way to another intron between exon 35 and exon 36
    # actually if I can check if the "last nucleotide of the directly upstream exon" of the start of the range
    # and the "first nucleotide of the directly downstream exon" of the end of the range
    # differ by 1
    # if they do, then that deletion or duplication is contained in one intron
    # if not, it is not contained in one exon
    # if it involves a UTR (either start or end), then we will include it anyway, so don't need to check
    # since it does not match, we do not need to do anything about it
    first_dot_index = mutation_CDS.find('.')
    first_underscore_index = mutation_CDS.find('_')
    first_posneg_sign_index = mutation_CDS.find('+', first_dot_index,
                                                first_underscore_index)
    if first_posneg_sign_index == -1:
        first_posneg_sign_index = mutation_CDS.find('-', first_dot_index,
                                                    first_underscore_index)
    last_nucleotide_directly_upstream_exon_position = mutation_CDS[
                                                      first_dot_index + 1:first_posneg_sign_index]
    second_posneg_sign_index = mutation_CDS.find('+', first_underscore_index)
    if second_posneg_sign_index == -1:
        second_posneg_sign_index = mutation_CDS.find('-',
                                                     first_underscore_index)
    first_nucleotide_directly_downstream_exon_position = mutation_CDS[
                                                         first_underscore_index + 1:second_posneg_sign_index]
    if last_nucleotide_directly_upstream_exon_position != first_nucleotide_directly_downstream_exon_position and \
            float(
                last_nucleotide_directly_upstream_exon_position) + 1 == first_nucleotide_directly_downstream_exon_position:
        print(mutation_CDS)
        print(last_nucleotide_directly_upstream_exon_position)
        print(first_nucleotide_directly_downstream_exon_position)


def calculate_mutation_frequency(cell_type_num_tumour_dict,
                                 gene_cell_mutation_type_info_dict):
    # calculate mutation frequency
    for gene_cell_mutation_type, info in gene_cell_mutation_type_info_dict.items():
        gene_cell_mutation_type_list = gene_cell_mutation_type.split(';')
        cell_type = gene_cell_mutation_type_list[1]

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_tumour_gene_cell_type = len(set(info['id tumour'])) * 100
        # total number of mutant sample for this histology subtype 1
        num_tumour_total_cell_type = \
            cell_type_num_tumour_dict[cell_type]

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_tumour_gene_cell_type / num_tumour_total_cell_type
        info['frequency percentage'] = frequency_percentage
        info['num_tumour_gene_cell_type'] = len(set(info['id tumour']))
        info['num_tumour_total_cell_type'] = num_tumour_total_cell_type


def read_CNV_file(cosmic_CNV_file_name, cell_type_num_total_tumor_dict_CNV,
                  gene_cell_mutation_type_dict_CNV,
                  important_column_heading_list_CNV,
                  important_column_number_list_CNV):
    with open(cosmic_CNV_file_name) as CNV_file:
        csv_reader = csv.reader(CNV_file, delimiter=',')
        for row in csv_reader:
            # 2, 9, 10, 11, 12, 13.  17, 18, 19
            # gene name, primary histology and its subtypes (3), sample name,
            # study id, GRch, genomic coordinate
            if row[9] == 'lymphoid_neoplasm':
                gene_ensembl = row[2].split('_')
                gene_name = gene_ensembl[0]

                if row[10] == 'NS':
                    continue

                histology = row[10]
                if 'T_cell' in histology or 'anaplastic' in histology \
                        or 'lymphomatoid_papulosis' in histology \
                        or 'post_transplant_lymphoproliferative_disorder' in histology \
                        or 'mycosis_fungoides-Sezary_syndrome' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'CNV'
                    cell_type = 'T_cell'
                elif 'B_cell' in histology or 'Burkitt' in histology or 'chronic' in histology \
                        or 'follicular_lymphoma' in histology or 'hairy' in histology \
                        or 'hodgkin' in histology or 'lymphoplasmacytic_lymphoma' in histology \
                        or 'mantle' in histology or 'marginal' in histology \
                        or 'plasma_cell_myeloma' in histology or 'plasmacytoma' in histology \
                        or 'effusion' in histology or 'MALT' in histology:
                    dict_key_name = gene_name + ';' + 'B_cell' + ';' + 'CNV'
                    cell_type = 'B_cell'
                else:
                    dict_key_name = gene_name + ';' + 'Other' + ';' + 'CNV'
                    cell_type = 'Other'

                specific_gene_histology_dict_CNV = gene_cell_mutation_type_dict_CNV.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list_CNV)):
                    column_heading = important_column_heading_list_CNV[i]
                    column_num = important_column_number_list_CNV[i]
                    specific_gene_histology_dict_CNV.setdefault(
                        column_heading, []).append(row[column_num])

                cell_type_num_total_tumor_dict_CNV.setdefault(cell_type,
                                                              []).append(row[4])

            # counter += 1
            # print(str(counter*100/650643) + '%')

    for cell_type, tumour_list in cell_type_num_total_tumor_dict_CNV.items():
        cell_type_num_total_tumor_dict_CNV[cell_type] = len(set(tumour_list))


def calculate_CNV_frequency(cell_type_num_tumour_dict_CNV,
                            gene_cell_mutation_type_info_dict_CNV):
    # calculate CNV frequency
    for gene_cell_mutation_type, info in gene_cell_mutation_type_info_dict_CNV.items():
        genn_cell_mutation_type_list = gene_cell_mutation_type.split(';')
        cell_type = genn_cell_mutation_type_list[1]

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_tumour_gene_cell_type = len(set(info['id tumour'])) * 100
        # total number of mutant sample for this histology subtype 1
        num_tumour_total_cell_type = \
            cell_type_num_tumour_dict_CNV[cell_type]

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_tumour_gene_cell_type / num_tumour_total_cell_type
        info['frequency percentage'] = frequency_percentage
        info['num_tumour_gene_cell_type'] = len(set(info['id tumour']))
        info['num_tumour_total_cell_type'] = num_tumour_total_cell_type


def convert_dict_list(gene_cell_type_dict,
                      l_chip_list):
    # first group different gene_cell_type into gene
    gene_dict = {}

    for gene_cell_type, info in gene_cell_type_dict.items():
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        cell_type = gene_cell_type_list[1]
        mutation_type = gene_cell_type_list[2]

        specific_gene_dict = gene_dict.setdefault(gene, {})
        specific_gene_dict[cell_type + ' ' + mutation_type] = info

    type_cell_mutation_list = ['T_cell point', 'T_cell CNV', 'B_cell point',
                               'B_cell CNV', 'Other point', 'Other CNV']
    for gene, cell_type_info_dict in gene_dict.items():
        l_chip_sublist = [gene]
        # each cell_type_info is a dict

        for type_cell_mutation in type_cell_mutation_list:
            # if this combination of cell type and mutation type exist
            if type_cell_mutation in cell_type_info_dict:

                # for this cell type, get both non-CNV and CNV tumour ids
                single_type_cell_mutation_list = type_cell_mutation.split(' ')
                cell_type = single_type_cell_mutation_list[0]

                # get cell_type + ' ' + 'point' from cell_type_info_dict
                #
                # if it does not exist, return empty dict
                # getting 'id tumour' from empty dict would not work
                # so empty set is returned
                #
                # if cell_type + ' ' + 'point' does exist, then the info dict
                # is returned, getting 'id tumour' from it would work
                # ,so it is a list of non-unique tumour ids
                # then I make it unique by turn the list into a set
                unique_non_CNV_tumour_ids = set(
                    cell_type_info_dict.get(cell_type + ' ' + 'point', {}).get(
                        'id tumour', set())
                )
                unique_CNV_tumour_ids = set(
                    cell_type_info_dict.get(cell_type + ' ' + 'CNV', {}).get(
                        'id tumour', set())
                )

                # is in non_CNV, but not in CNV (non_CNV minus intersection)
                only_non_CNV = len(unique_non_CNV_tumour_ids.difference(
                    unique_CNV_tumour_ids))
                # is in CNV, but not in non_CNV (CNV minus intersection)
                only_CNV = len(unique_CNV_tumour_ids.difference(
                    unique_non_CNV_tumour_ids))
                both = len(unique_CNV_tumour_ids.intersection(
                    unique_non_CNV_tumour_ids))

                info = cell_type_info_dict[type_cell_mutation]
                l_chip_sublist.append(info['frequency percentage'])
                l_chip_sublist.append(info['num_tumour_gene_cell_type'])
                l_chip_sublist.append(info['num_tumour_total_cell_type'])
                l_chip_sublist.append(only_non_CNV)
                l_chip_sublist.append(only_CNV)
                l_chip_sublist.append(both)
                l_chip_sublist.append(info['mutation genome position'])
            else:
                l_chip_sublist.append(0)
                l_chip_sublist.append(0)
                l_chip_sublist.append('n/a')
                # above is, in reality this is the same as all others
                l_chip_sublist.append(0)
                l_chip_sublist.append(0)
                l_chip_sublist.append(0)
                l_chip_sublist.append('n/a')

        if 'T_cell CNV' not in cell_type_info_dict \
                and 'B_cell CNV' not in cell_type_info_dict \
                and 'Other CNV' not in cell_type_info_dict:
            l_chip_sublist.append('not')
        elif 'T_cell point' not in cell_type_info_dict \
                and 'B_cell point' not in cell_type_info_dict \
                and 'Other point' not in cell_type_info_dict:
            l_chip_sublist.append('CNV')
        else:
            l_chip_sublist.append('both')

            # TODO other stuff we want in the table

        l_chip_list.append(l_chip_sublist)


def write_output(gene_set, file_name):
    df = pd.DataFrame(gene_set)
    writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='gene set', index=False)
    writer.save()


def which_genes_were_mutated_in_healthy(l_chip_gene_set_1_all,
                                        healthy_mutation_file):
    gene_list = []
    with open(healthy_mutation_file) as healthy_mutation_file:
        lines = healthy_mutation_file.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
        gene_list.extend(stripped_lines)
    l_chip_gene_set_1_healthy_mutation = set(gene_list)

    return l_chip_gene_set_1_all.intersection(
        l_chip_gene_set_1_healthy_mutation)


def separate_gene_cell_mutation_type_info_dict(cell_mutation_type_dict,
                                               gene_cell_mutation_type_info_dict):
    for gene_cell_mutation_type, info in gene_cell_mutation_type_info_dict.items():
        gene_cell_mutation_type_list = gene_cell_mutation_type.split(';')
        gene = gene_cell_mutation_type_list[0]
        cell_type = gene_cell_mutation_type_list[1]
        mutation_type = gene_cell_mutation_type_list[2]
        cell_mutation_type = cell_type + ';' + mutation_type
        cell_mutation_type_dict[cell_mutation_type][gene] = info


def categorize_then_choose_probe_placement_within(
        gene_cell_mutation_type_info_dict,
        targeting_window_size=80,
        cumulative_contribution_threshold=90,
        indel_filter_threshold=30,
        recurrent_definition=2):
    cell_mutation_type_dict = {
        'T_cell;point': {}, 'B_cell;point': {}, 'Other;point': {},
    }
    type_probe_dict = {
        'T_cell;point': {}, 'B_cell;point': {}, 'Other;point': {},
    }

    separate_gene_cell_mutation_type_info_dict(cell_mutation_type_dict,
                                               gene_cell_mutation_type_info_dict)

    for cell_mutation_type, gene_info_dict in cell_mutation_type_dict.items():

        rank_all_probe_choose_most_recurrent(
            gene_info_dict, type_probe_dict, cell_mutation_type,
            cumulative_contribution_threshold, targeting_window_size,
            indel_filter_threshold, recurrent_definition)

        chose_most_recurrent_mutation_then_center_probe(
            gene_info_dict, type_probe_dict, cell_mutation_type,
            cumulative_contribution_threshold, targeting_window_size,
            indel_filter_threshold, recurrent_definition)

    for cell_mutation_type, probe_info in type_probe_dict.items():
        # total number of probe to cover at the threshold
        num_probe_percentage_cover_threshold = probe_info[0]
        # a list of list with location and cumulative contribution
        probe_location_cc_list = probe_info[1]
        print("we need", num_probe_percentage_cover_threshold,
              "probes to cover",
              cumulative_contribution_threshold, "% of the tumours in",
              cell_mutation_type, ", their locations are:")
        for probe_location_cc in probe_location_cc_list:
            print('\t', probe_location_cc)

        write_output(probe_location_cc_list, cell_mutation_type + 'probe.xlsx')


def chose_most_recurrent_mutation_then_center_probe(
        gene_info_dict, type_probe_dict, cell_mutation_type,
        cumulative_contribution_threshold, targeting_window_size,
        indel_filter_threshold, recurrent_definition):
    all_tumour_id = []
    gene_tumour_position = {}

    # group tumour id based on genomic position
    get_mapping_position_tumour(all_tumour_id, gene_tumour_position,
                                gene_info_dict,
                                indel_filter_threshold)

    # recurrent_definition,
    #                                        targeting_window_size

    # each gene has one dict holding position mapping to tumour
    all_position_tumour_dict = {}
    for gene_cell_type, position_tumour_dict_list in gene_tumour_position.items():

        position_tumour = position_tumour_dict_list[0]
        position_tumour_deletion = position_tumour_dict_list[2]

        for position, tumour_set in position_tumour:
            all_position_tumour_dict[position_tumour] = tumour_set

        for position, tumour_set in position_tumour_deletion:
            all_position_tumour_dict[position_tumour] = tumour_set

    tumour_set_list = []
    unique_tumour_id = []
    recurrent_unique_tumour_id = count_mutation_and_tumour_after_filters(
        gene_tumour_position, recurrent_definition, tumour_set_list,
        unique_tumour_id)

    print('start finding cover')
    covered_tumour_id = []
    cover_sizes_percentage = [0]
    percentage_cover_threshold = len(
        recurrent_unique_tumour_id) * cumulative_contribution_threshold / 100
    mutation_list = []
    while len(covered_tumour_id) < percentage_cover_threshold:

        # find the position that currently has the most tumour ids
        num_most_tumour = 0
        position_most_tumour = ''
        for position, tumour_set in all_position_tumour_dict.items():
            if len(tumour_set[0]) > num_most_tumour:
                num_most_tumour = len(tumour_set[0])
                position_most_tumour = position

        # get its tumour set
        tumour_set_selected_position = all_position_tumour_dict[
            position_most_tumour]

        # remove that position
        all_position_tumour_dict.pop(position_most_tumour, None)

        # get center of the mutation, then start search within the window size, where the center is
        # that center of mutation
        position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
            position_most_tumour)
        mutation_center = (position_range_start + position_range_end)//2

        probe_start = mutation_center - (targeting_window_size // 2)

        probe_end = mutation_center + (targeting_window_size // 2)

        tumour_set_selected_position_and_probe = []
        tumour_set_selected_position_and_probe.extend(tumour_set_selected_position)
        all_position_in_probe = []
        for i in range(probe_start, probe_end):
            a_position_in_probe = i
            a_position_in_probe_tumour_set = all_position_tumour_dict.get(i)
            if a_position_in_probe_tumour_set is not None and a_position_in_probe_tumour_set >= recurrent_definition:
                all_position_in_probe.append(a_position_in_probe)
                tumour_set_selected_position_and_probe.append(a_position_in_probe_tumour_set)
                all_position_tumour_dict.pop(i, None)

        all_position_in_probe.sort()



        # TODO: get a probe that is center at selected position
        # for window size, if >= recurrent_definition
        # add tumour from those position as well
        # add those tumour ids too, and subtract those tumour ids too
        first_mutation_selected_position_and_probe = all_position_in_probe[0]
        last_mutation_selected_position_and_probe = all_position_in_probe[-1]

        # TODO:
        # print first_mutation_selected_probe
        # and see if you can add one base pair to it
        # or maybe when parsing just add one (where I filter indels)

        # add the newly covered tumour ids, and keep it unique
        covered_tumour_id.extend(list(tumour_set_selected_position_and_probe))
        covered_tumour_id = list(set(covered_tumour_id))
        # calculate new cumulative contribution
        cumulative_contribution = len(covered_tumour_id) / len(
            recurrent_unique_tumour_id)

        # position most tumour is represented by its central base pair
        mutation_list.append([position_most_tumour, cumulative_contribution,
                              first_mutation_selected_position_and_probe,
                              last_mutation_selected_position_and_probe])
        cover_sizes_percentage.append(cumulative_contribution)
        print(cumulative_contribution)

        # for every other set, subtract tumour ids in current, add to new list
        for position, tumour_set in all_position_tumour_dict.items():
            all_position_tumour_dict[position] = tumour_set.difference(
                tumour_set_selected_position_and_probe)

    type_probe_dict[cell_mutation_type] = [percentage_cover_threshold,
                                           mutation_list]

    plt.plot(cover_sizes_percentage, linewidth=1, label=cell_mutation_type)
    plt.ylabel('percentage of tumour covered')
    plt.xlabel('number of ' + targeting_window_size + 'bp probes selected')
    plt.title('coverage of tumours with increasing probe ')
    plt.savefig('stop_at_90.pdf')


def rank_all_probe_choose_most_recurrent(
        gene_info_dict, type_probe_dict, cell_mutation_type,
        cumulative_contribution_threshold, targeting_window_size,
        indel_filter_threshold, recurrent_definition):
    all_tumour_id = []
    gene_tumour_position = {}

    # group tumour id based on genomic position
    get_mapping_position_tumour(all_tumour_id, gene_tumour_position,
                                gene_info_dict,
                                indel_filter_threshold)

    all_possible_probe_tumour_dict = {}
    find_all_potential_probe_placement(gene_tumour_position,
                                       all_possible_probe_tumour_dict,
                                       recurrent_definition,
                                       targeting_window_size)

    tumour_set_list = []
    unique_tumour_id = []
    recurrent_unique_tumour_id = count_mutation_and_tumour_after_filters(
        gene_tumour_position, recurrent_definition, tumour_set_list,
        unique_tumour_id)

    print('total number of possible probes',
          len(all_possible_probe_tumour_dict))

    print('start finding cover')
    covered_tumour_id = []
    cover_sizes_percentage = [0]
    percentage_cover_threshold = len(
        recurrent_unique_tumour_id) * cumulative_contribution_threshold / 100
    probe_list = []
    while len(covered_tumour_id) < percentage_cover_threshold:

        # the tumour ids of the position that current has the most tumour ids
        num_most_tumour = 0
        probe_most_tumour = ''
        for probe, tumour_set in all_possible_probe_tumour_dict.items():
            if len(tumour_set[0]) > num_most_tumour:
                num_most_tumour = len(tumour_set[0])
                probe_most_tumour = probe
        tumour_set_selected_probe = all_possible_probe_tumour_dict[
            probe_most_tumour][0]
        first_mutation_selected_probe = all_possible_probe_tumour_dict[
            probe_most_tumour][1]
        last_mutation_selected_probe = all_possible_probe_tumour_dict[
            probe_most_tumour][2]
        all_possible_probe_tumour_dict.pop(probe_most_tumour, None)

        # TODO:
        # print first_mutation_selected_probe
        # and see if you can add one base pair to it

        # add the newly covered tumour ids, and keep it unique
        covered_tumour_id.extend(list(tumour_set_selected_probe))
        covered_tumour_id = list(set(covered_tumour_id))
        # calculate new cumulative contribution
        cumulative_contribution = len(covered_tumour_id) / len(
            recurrent_unique_tumour_id)

        # probe most tumour is represented by its starting base pair
        probe_list.append([probe_most_tumour, cumulative_contribution,
                           first_mutation_selected_probe,
                           last_mutation_selected_probe])
        cover_sizes_percentage.append(cumulative_contribution)
        print(cumulative_contribution)

        # for every other set, subtract tumour ids in current, add to new list
        for probe, tumour_set in all_possible_probe_tumour_dict.items():
            all_possible_probe_tumour_dict[probe][0] = tumour_set[0].difference(
                tumour_set_selected_probe)

    type_probe_dict[cell_mutation_type] = [percentage_cover_threshold,
                                           probe_list]

    plt.plot(cover_sizes_percentage, linewidth=1, label=cell_mutation_type)
    plt.ylabel('percentage of tumour covered')
    plt.xlabel('number of ' + targeting_window_size + 'bp probes selected')
    plt.title('coverage of tumours with increasing probe ')
    plt.savefig('stop_at_90.pdf')

    "so this is actually a set cover problem, whose naive greedy algorithm produces" \
    "a solution with the upper bound of OPT log_e n where OPT is the optimal" \
    "solution and the n the total number of element (size of the universe)" \
    "our unique n = 35566 , so about the optimal solution is at most 11 times smaller"
    "there is faster greedy algorithms, see wikipedia or CLRS"
    "the number of sets is 9915"
    "there is a 9:5073770 with 17333"

    "since we are trying to cover mutation with a probe of 100 bp, so we " \
    "changing the subset's definition to from 'all unique tumours at a bp' to" \
    "'all unique tumour within a 80bp window' this mean we increase the subset size" \
    "now the approximation ratio does not change since it depend on universe size," \
    "which has not changed, but the optimal solution has become smaller, how much smaller?" \
    "I don't know, not sure if it has to do with the fact that we dont know the" \
    "optimal solution to begin with "


def count_mutation_and_tumour_after_filters(gene_tumour_position,
                                            recurrent_definition,
                                            tumour_set_list, unique_tumour_id):
    # extract all sets of tumour, remove duplicates within each set, remove
    # duplicate set. that is convert to set for all 'subset', and convert to
    # multiset for all the collection of the subset remove duplicates and count
    # number of unique tumours this is to know the number of (recurrent) unique tumours
    for gene_cell_type, position_tumour_dict_list in gene_tumour_position.items():
        position_tumour = position_tumour_dict_list[0]
        position_tumour_deletion = position_tumour_dict_list[2]

        # add tumours ids from deletion mutation
        for deletion, tumour_ids in position_tumour_deletion.items():
            tumour_ids_set = list(set(tumour_ids))
            tumour_set_list.append(tumour_ids_set)
            unique_tumour_id.append(tumour_ids_set)

        # add tumours ids from other mutations
        for position, tumour_ids in position_tumour.items():
            tumour_ids_set = list(set(tumour_ids))
            tumour_set_list.append(tumour_ids_set)
            unique_tumour_id.extend(tumour_ids_set)
    unique_tumour_id = list(set(unique_tumour_id))
    unique_tumour_set_list = list(
        set(frozenset(item) for item in tumour_set_list))
    # we now have
    # a dict with genomic position to a list of tumour ids and its reverse
    # a list of set of tumour ids
    # recurrent mutation here is a base pair that is mutated in multiple tumour
    # since each tumour set is a set a tumour that is mutated at that position
    # keep only recurrent mutations, means only using sets with >=2 tumour
    # which is the default recurrent definition
    recurrent_unique_tumour_set_list = [tumour_set
                                        for tumour_set in unique_tumour_set_list
                                        if
                                        len(tumour_set) >= recurrent_definition]
    # then only keep these tumour ids
    recurrent_unique_tumour_id = list(set([tumour_id
                                           for tumour_set in
                                           recurrent_unique_tumour_set_list
                                           for tumour_id in tumour_set]))
    print('total number of mutations after indel filter',
          len(unique_tumour_set_list))
    print('total number of mutations after indel filter',
          len(unique_tumour_id))
    print('total number of recurrent mutations after indel filter',
          len(recurrent_unique_tumour_set_list))
    print(
        'total number of tumours with only recurrent mutations after indel filter',
        len(recurrent_unique_tumour_id))
    return recurrent_unique_tumour_id


def get_mapping_position_tumour(all_tumour_id, gene_tumour_position,
                                gene_cell_mutation_type_info_dict,
                                indel_filter_threshold):
    all_tumour_id_before_indel_filter = []
    num_mutation_before_indel_filter = 0
    for gene_cell_type, info in gene_cell_mutation_type_info_dict.items():

        for i in range(len(info['mutation genome position'])):
            tumour_id = info['id tumour'][i]
            all_tumour_id_before_indel_filter.append(tumour_id)

        num_mutation_before_indel_filter += len(
            info['mutation genome position'])

    all_tumour_id_before_indel_filter = list(
        set(all_tumour_id_before_indel_filter))

    print("total number of mutation before indel filter",
          num_mutation_before_indel_filter)
    print("total number of tumours before indel filter",
          all_tumour_id_before_indel_filter)

    for gene_cell_type, info in gene_cell_mutation_type_info_dict.items():
        # assert len(info['mutation genome position']) == len(info['id tumour'])

        position_tumour = {}
        tumour_position = {}
        position_tumour_deletion = {}
        gene_tumour_position[gene_cell_type] = [position_tumour,
                                                tumour_position,
                                                position_tumour_deletion]
        for i in range(len(info['mutation genome position'])):
            chromosome_position_range = info['mutation genome position'][i]
            tumour_id = info['id tumour'][i]
            if chromosome_position_range != 'null' and tumour_id != 'null':

                # all genomic position are in the format a:b-c
                # where a is the chromosome number, followed by a semicolon
                # b is the start of the range, followed by a dash
                # c is the end of the range, followed by nothing
                # for single base pair mutation it is denoted as having the
                # same start and end of the range
                # for insertion, it is one more
                # deletion it is a lot more
                position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
                    chromosome_position_range)

                # remove mutation with greater than then indel filter threshold (default is 30)
                if int(position_range_end) > int(
                        position_range_start) + indel_filter_threshold:
                    continue

                if int(position_range_end) > int(position_range_start) + 1:
                    position_tumour_deletion[
                        chromosome_position_range] = tumour_id

                all_tumour_id.append(tumour_id)
                chromosome = chromosome_position_range[:semicolon_index]
                for position in position_list:
                    chromosome_position = chromosome + ':' + str(position)
                    position_tumour.setdefault(chromosome_position, []).append(
                        tumour_id)
                    tumour_position.setdefault(tumour_id, []).append(
                        chromosome_position)


def parse_chromosome_position_range(chromosome_position_range):
    # example 2:3241-3276
    semicolon_index = chromosome_position_range.index(':')
    dash_index = chromosome_position_range.index('-')
    position_range_start = chromosome_position_range[
                           semicolon_index + 1:dash_index]
    position_range_end = chromosome_position_range[dash_index + 1:]
    position_list = list(range(int(position_range_start),
                               int(position_range_end) + 1))
    return position_list, position_range_end, position_range_start, semicolon_index


def find_all_potential_probe_placement(gene_tumour_position,
                                       all_possible_probe_tumour_dict,
                                       recurrent_definition,
                                       targeting_window_size):
    counter = 1
    for gene_cell_type, position_tumour_dict_list in gene_tumour_position.items():

        position_tumour = position_tumour_dict_list[0]
        position_tumour_deletion = position_tumour_dict_list[2]

        # make a list of all position with mutation, sorted to be use in the next section
        all_position_list = []
        for position, tumour in position_tumour.items():
            all_position_list.append(position)
        all_position_list.sort()

        # make a list that denotes the number of unique tumour at all position
        # with or without mutation
        position_tumour_set_list = []
        all_position_list_with_gap = []
        position_num_tumour_list = []
        position_tumour_deletion_dict = {}
        make_position_info_list(all_position_list, all_position_list_with_gap,
                                position_num_tumour_list, position_tumour,
                                position_tumour_set_list,
                                position_tumour_deletion,
                                position_tumour_deletion_dict)

        index_of_mutation = []
        for i in range(len(position_num_tumour_list)):
            # only if the mutation is recurrent, add it in
            if position_num_tumour_list[i] >= recurrent_definition:
                index_of_mutation.append(i)

        # calculate what tumour will be in each probe
        # another list to signify the index where the mutation is
        #   not empty, so that we only include probes that are within 80bp of
        #   each side for every mutation
        probe_tumour_set = []
        make_each_probe(all_position_list_with_gap,
                        all_possible_probe_tumour_dict, index_of_mutation,
                        position_tumour_set_list, probe_tumour_set,
                        targeting_window_size, recurrent_definition,
                        position_tumour_deletion, position_tumour_deletion_dict)

        # the position in all_possible_probe_tumour_dict represent the
        # start of a probe, unlike position_tumour dict

        print('nth gene', counter, len(gene_tumour_position))
        counter += 1


def make_position_info_list(all_position_list, all_position_list_with_gap,
                            position_num_tumour_list, position_tumour,
                            position_tumour_set_list, position_tumour_deletion,
                            position_tumour_deletion_dict):
    # i want a dictionary of all position of any deletion to empty string
    # I need the covered base pairs of the deletion
    for deletion_position, tumour_set in position_tumour_deletion.items():
        position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
            deletion_position)
        for position in position_list:
            position_tumour_deletion_dict[position] = deletion_position

    for i in range(len(all_position_list)):
        last_position = all_position_list[i]
        tumour_set = position_tumour[last_position]
        position_tumour_set_list.append(tumour_set)

        position_num_tumour_list.append(len(tumour_set))

        all_position_list_with_gap.append(last_position)

        if len(all_position_list) == 1:
            break

        # if we have reach the end of the end, do not add gaps, because there
        # are not meant to be any
        if i + 1 == len(all_position_list):
            break

        # if there is a gap of X, append X 0s
        next_position = all_position_list[i + 1]
        gap_length = find_gap_length(next_position, last_position)

        gap_list = ['0' for _ in range(gap_length)]
        position_gap_list = ['' for _ in range(gap_length)]
        num_gap_list = [0] * gap_length

        # a list of tumour set for each position, except where there is
        # no tumour at that position, then it has a '0'
        position_tumour_set_list.extend(gap_list)
        all_position_list_with_gap.extend(position_gap_list)
        position_num_tumour_list.extend(num_gap_list)


def make_each_probe(all_position_list_with_gap, all_possible_probe_tumour_dict,
                    index_of_mutation, position_tumour_set_list,
                    probe_tumour_set, targeting_window_size,
                    recurrent_definition, position_tumour_deletion,
                    position_tumour_deletion_dict):
    # a unique list of list of mutation indices that is covered by each probe
    # list_position_set = []
    for current_position in range(len(position_tumour_set_list)):

        # TODO the function below is the bottle neck (that it ran for so many times)
        # this alone is 1105.556s out of 1653.824s (18.4 min out of 27min)
        # the rest of make each probe is 221
        # it was called 969 million times, but the rest of make each probe is called only 3m times
        # we need pre-filter this to reduce bottleneck
        # check if the current position is at most targeting_window_size before a mutation
        between, index_of_first_mutation_covered = is_current_position_a_probe_length_before_mutation(
            current_position, index_of_mutation, targeting_window_size)

        # skip, if it is not zero and not between
        if current_position != 0 and not between:
            continue

        # I think this was the rest of make each probe
        # # go from the index of the first mutation index covered to the index of the last mutation that
        # # can be covered (if there is one)
        # # a list of mutation indices covered by this probe
        # list_mutation_index_covered_probe = []
        # end_of_probe = current_position + 80
        # for k in range(index_of_first_mutation_covered, len(index_of_mutation)):
        #     if index_of_mutation[k] > end_of_probe:
        #         break
        #     list_mutation_index_covered_probe.append(k)
        #
        # # if this potential probe covers each the same mutation indices as a previous probe, skip
        # # making this probe
        # if list_mutation_index_covered_probe in list_position_set:
        #     continue
        # else:
        #     list_position_set.append(list_mutation_index_covered_probe)

        # no, this does not matter, since this is just probe starting bp
        # if the number of unique tumours at this position is 1 or 0
        # if len(position_tumour_set_list[current_position]) < 2:
        #     continue

        # the starting base pair of this probe
        probe_start_bp = all_position_list_with_gap[current_position]

        all_tumour_probe, first_current_mutation_location, last_recurrent_mutation_location \
            = collect_probe_unique_tumours(current_position,
                                           position_tumour_set_list,
                                           targeting_window_size,
                                           recurrent_definition,
                                           position_tumour_deletion,
                                           position_tumour_deletion_dict,
                                           probe_start_bp)

        # TODO only need to extend them by the index with mutations
        # i need everything from the last all_tumour_probe
        # except the one from the first set

        # first_tumour_set = set(probe_gene_section[i])
        #
        # only_in_first = set(all_tumour_probe) - first_tumour_set

        all_possible_probe_tumour_dict[probe_start_bp] = (all_tumour_probe,
                                                          first_current_mutation_location,
                                                          last_recurrent_mutation_location)
        probe_tumour_set.append(all_tumour_probe)


def is_current_position_a_probe_length_before_mutation(current_position,
                                                       index_of_mutation,
                                                       targeting_window_size):
    between = False
    index_of_first_mutation_covered = 0
    for j in range(len(index_of_mutation)):
        a_mutation_position = index_of_mutation[j]
        if a_mutation_position - targeting_window_size <= current_position <= a_mutation_position:
            between = True
            index_of_first_mutation_covered = j
            break
    return between, index_of_first_mutation_covered


def collect_probe_unique_tumours(current_position, position_tumour_set_list,
                                 targeting_window_size, recurrent_definition,
                                 position_tumour_deletion,
                                 position_tumour_deletion_dict,
                                 probe_start_bp):
    """
    :param probe_start_bp:
    :param position_tumour_deletion_dict:
    :param position_tumour_deletion:
    :param current_position: the starting position of this potential probe
    :param position_tumour_set_list: a list of sets that have of tumours id of that position
    :param targeting_window_size: the number of base pair our probe can target
    :param recurrent_definition: the number of samples for a mutation to be considered a recurrent mutation
    count the unique tumours that is covered by this potential probe
    """

    # the section of this gene that corresponds to this probe
    end_of_probe = current_position + targeting_window_size
    probe_gene_section = position_tumour_set_list[current_position:end_of_probe]
    # all unique tumour covered by this probe
    all_tumour_probe = []

    # if any base of the probe covers the deletion
    probe_end_bp = probe_start_bp + targeting_window_size
    deletion_start_bp_list = []
    deletion_end_bp_list = []
    for i in range(probe_start_bp, probe_end_bp):
        # if this position is covered by one of the deletion
        if i in position_tumour_deletion_dict:
            deletion_position = position_tumour_deletion_dict[i]
            # parse the chromosome position
            position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
                deletion_position)
            # if yes, add its tumour to this probe
            if position_range_end < probe_end_bp:
                all_tumour_probe.extend(
                    position_tumour_deletion[deletion_position])

            deletion_start_bp_list.append(position_range_start)
            deletion_end_bp_list.append(position_range_end)

    deletion_start_bp_list = sorted(list(set(deletion_start_bp_list)))
    deletion_end_bp_list = sorted(list(set(deletion_end_bp_list)))

    position_counter = 0
    subin_location_list = []
    for tumour_set in probe_gene_section:
        position_counter += 1
        if tumour_set == '0':
            continue
        if len(tumour_set) < recurrent_definition:
            continue
        all_tumour_probe.extend(tumour_set)
        subin_location_list.append(position_counter)

    first_recurrent_mutation_location = min(subin_location_list[0],
                                            deletion_start_bp_list[0])
    last_recurrent_mutation_location = max(subin_location_list[-1],
                                           deletion_end_bp_list[-1])
    all_tumour_probe = set(all_tumour_probe)
    return all_tumour_probe, first_recurrent_mutation_location, last_recurrent_mutation_location


def find_gap_length(next_position, last_position):
    next_position_list = next_position.split(':')
    last_position_list = last_position.split(':')
    next_position_wo_chromosome = next_position_list[1]
    last_position_wo_chromosome = last_position_list[1]
    return int(next_position_wo_chromosome) - int(
        last_position_wo_chromosome) - 1


if __name__ == "__main__":
    # mutation file, CNV file, and gene list file
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
