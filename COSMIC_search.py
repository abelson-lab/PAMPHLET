import csv
import sys
import re
import pandas as pd


# how about, I take the cosmic files, calculate the frequency, then remove CNV
# separate into 2 tsv outputs, one of m-chip vs l-chip paper,
# one for every other gene


def rank_cosmic_rows(cosmic_mutation_file_name, cosmic_CNV_file_name,
                     gene_list_file, healthy_mutation_file_name):
    # read known l-chip gene lists
    # TODO this is where you add more genes, from other papers
    gene_list = []
    with open(gene_list_file) as gene_list_file:
        lines = gene_list_file.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
        gene_list.extend(stripped_lines)
    l_chip_gene_set_1 = set(gene_list)

    # these two should match
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

    gene_cell_mutation_type_info_dict = dict()
    cell_type_num_tumour_dict = {}

    # read mutation, group by gene+cell+mutation type (in gene_cell_mutation_type_info_dict)
    # and count how many sample each cell type has (in cell_type_num_tumour_dict)
    read_mutation_file(cosmic_mutation_file_name,
                       cell_type_num_tumour_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_cell_mutation_type_info_dict)

    print('gene cell-type pair non-CNV', len(gene_cell_mutation_type_info_dict))
    print('non-CNV cell type distribution', cell_type_num_tumour_dict)

    # and store in gene_cell_mutation_type_info_dict
    calculate_mutation_frequency(
        cell_type_num_tumour_dict,
        gene_cell_mutation_type_info_dict)

    # filter out those not found in at least 2 studies
    # TODO change this to in more than 1 tumour
    # reproducible = {}
    # filter_non_reproducible(gene_cell_mutation_type_info_dict, reproducible)
    # print('reproducible', len(reproducible))
    reproducible = gene_cell_mutation_type_info_dict

    # read CNV
    gene_cell_mutation_type_info_dict_CNV = {}
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

    # add CNV info into it
    # TODO can there be CNV that have more than 1 study but still
    #   does not have more than 1 in point?
    all_possible_l_chip_dict = {}
    for gene_cell_type, info in reproducible.items():
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

    # TODO extract into functions
    # separate reproducible into 2 outputs files
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

    print('known l chip gene cell-type pair', len(known_l_chip_dict))
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

    print('remaining', len(remaining_genes))

    # turn dict into list, and sort them based on frequency
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

    all_possible_l_chip_list = []
    convert_dict_list(all_possible_l_chip_dict,
                      all_possible_l_chip_list)

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

    sorted_all_possible_l_chip_list = sorted(all_possible_l_chip_list,
                                             key=lambda x: float(x[1]),
                                             reverse=True)

    print('num all genes', len(sorted_all_possible_l_chip_list))

    # write into files
    write_output(
        potential_l_chip_list_headings + sorted_all_possible_l_chip_list,
        'potential_l_chip.xlsx')
    write_output(remaining_genes, 'remaining_genes.xlsx')


def read_mutation_file(cosmic_mutation_file_name,
                       cell_type_num_total_tumour_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_cell_mutation_type_dict):
    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')

        for row in csv_reader:
            # 0, 4, 11, 12, 13, 14.  19, 21, 24, 25, 30
            # gene name, sample name, primary histology and its subtypes (3),
            # mutation CDS, mutation desc, GRch, mutation genome position, study id,
            # if primary histology is lymphoid neoplasm

            if row[11] == 'lymphoid_neoplasm':
                # group mutations by gene name + histology subtype 1
                gene_ensembl = row[0].split('_')
                gene_name = gene_ensembl[0]

                mutation_CDS = row[19]
                # assert 'c.' in mutation_CDS
                # check_variant_type(row)

                # we don't need this, there are no such mutation in this file
                # filter_out_flag_intron_range = is_intron_range_mutation(
                #     mutation_CDS)
                #
                # if filter_out_flag_intron_range:
                #     check_intronic_mutation_diff_intron(mutation_CDS)

                # if it is intron (where the position is known) and not splice site and not
                # ending or start in an exon, and both start and end in two different introns skip it
                filter_out_flag = filter_known_pure_intronic_mutation(
                    mutation_CDS)
                if filter_out_flag:
                    continue

                if row[12] == 'NS':
                    continue

                if row[21] == 'Substitution - coding silent':
                    continue

                histology = row[12]
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


def filter_non_reproducible(lymphoid_neoplasm_gene_histology_dict,
                            reproducible):
    for gene_histology, info in lymphoid_neoplasm_gene_histology_dict.items():
        # for this gene, cell_type if it has more than 2 unique studies
        if len(set(info['study id'])) >= 2:
            reproducible[
                gene_histology] = info


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
                    cell_type_info_dict.get(cell_type + ' ' + 'point', {})
                        .get('id tumour', set())
                )
                unique_CNV_tumour_ids = set(
                    cell_type_info_dict.get(cell_type + ' ' + 'CNV', {})
                        .get('id tumour', set())
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


if __name__ == "__main__":
    # mutation file, CNV file, and gene list file
    rank_cosmic_rows(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
