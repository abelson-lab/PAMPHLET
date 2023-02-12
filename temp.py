import pandas as pd



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


def integrate_CNV_point_mutation_info(all_possible_l_chip_dict,
                                      gene_cell_mutation_type_info_dict,
                                      gene_cell_mutation_type_info_dict_CNV):
    # add CNV info into it
    # FIXME can there be CNV that have more than 1 study but still
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


def convert_dict_list(gene_cell_type_dict,
                      l_chip_list):
    # first group different gene_cell_type into gene
    gene_dict = {}

    for gene_cell_type, info in gene_cell_type_dict.items():
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        mutation_type = gene_cell_type_list[1]

        specific_gene_dict = gene_dict.setdefault(gene, {})
        specific_gene_dict[mutation_type] = info

    type_cell_mutation_list = ['point', 'CNV',]
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
                l_chip_sublist.append(info['num_tumour_gene'])
                l_chip_sublist.append(info['num_tumour_total'])
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


def write_output(gene_set, file_name):
    df = pd.DataFrame(gene_set)
    writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='gene set', index=False)
    writer.save()


"""the rest are 3rd level or beyond"""


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
