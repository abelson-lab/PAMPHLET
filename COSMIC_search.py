import csv
import sys


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

    print('gene-cell type pair', len(gene_cell_mutation_type_info_dict))
    print(cell_type_num_tumour_dict)

    # and store in gene_cell_mutation_type_info_dict
    calculate_mutation_frequency(
        cell_type_num_tumour_dict,
        gene_cell_mutation_type_info_dict)

    # filter out those not found in at least 2 studies
    reproducible = {}
    filter_non_reproducible(gene_cell_mutation_type_info_dict, reproducible)
    print('reproducible', len(reproducible))

    # read CNV
    gene_cell_mutation_type_info_dict_CNV = {}
    CNV_gene_list = []
    cell_type_num_tumor_dict_CNV = {}
    read_CNV_file(cosmic_CNV_file_name,
                  cell_type_num_tumor_dict_CNV,
                  CNV_gene_list,
                  gene_cell_mutation_type_info_dict_CNV,
                  important_column_heading_list_CNV,
                  important_column_number_list_CNV)

    print(cell_type_num_tumor_dict_CNV)

    # and store in gene_cell_mutation_type_info_dict_CNV
    calculate_CNV_frequency(cell_type_num_tumor_dict_CNV,
                            gene_cell_mutation_type_info_dict_CNV)

    lymphoid_neoplasm_CNV_gene_set = set(CNV_gene_list)

    # add CNV info into it
    # TODO can there be CNV that have more than 1 study but still
    #   does not have more than 1 in point?
    reproducible_and_CNV = {}
    for gene_cell_type, info in reproducible.items():
        gene_cell_type_list = gene_cell_type.split(';')
        gene = gene_cell_type_list[0]
        cell_type = gene_cell_type_list[1]
        reproducible_and_CNV[gene_cell_type] = info
        CNV_key = gene + ';' + cell_type + ';' + 'CNV'
        # if this gene has a CNV
        if CNV_key in lymphoid_neoplasm_CNV_gene_set:
            # take first two part of this gene (gene name + cell type)
            # add 'CNV' to it, to get the CNV version of it
            this_gene_info = gene_cell_mutation_type_info_dict_CNV[
                CNV_key]
            reproducible_and_CNV[CNV_key] = this_gene_info

    # TODO extract into functions
    # separate reproducible into 2 outputs files
    # known: in the gene set
    # potential: not in the gene set
    known_l_chip_dict = {}
    potential_l_chip_dict = {}
    for gene_histology, info in reproducible_and_CNV.items():
        gene_histology_list = gene_histology.split(';')
        gene = gene_histology_list[0]
        if gene in l_chip_gene_set_1:
            known_l_chip_dict[gene_histology] = info
        else:
            potential_l_chip_dict[gene_histology] = info

    print('known l chip', len(known_l_chip_dict))
    print('potential l chip', len(potential_l_chip_dict))

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
    known_l_chip_list_headings = [[
        'gene',
        'T Cell Frequency Point', 'tumour num', 'total tumour num',
        'B Cell Frequency Point', 'tumour num', 'total tumour num',
        'Other Cell Frequency Point', 'tumour num', 'total tumour num',
        'T Cell Frequency CNV', 'tumour num', 'total tumour num',
        'B Cell Frequency CNV', 'tumour num', 'total tumour num',
        'Other Cell Frequency CNV', 'tumour num', 'total tumour num',
        'CNV or not', 'found in healthy'
    ]]

    potential_l_chip_list_headings = known_l_chip_list_headings

    known_l_chip_list = []
    convert_dict_list(known_l_chip_dict,
                      known_l_chip_list)
    potential_l_chip_list = []
    convert_dict_list(potential_l_chip_dict,
                      potential_l_chip_list)

    genes_in_healthy = which_genes_were_mutated_in_healthy(l_chip_gene_set_1, healthy_mutation_file_name)
    for i in range(len(known_l_chip_list)):
        gene_info = known_l_chip_list[i]
        gene = gene_info[0]
        if gene in genes_in_healthy:
            known_l_chip_list[i].append('healthy')
        else:
            known_l_chip_list[i].append('___')




    sorted_known_l_chip_list_T = sorted(known_l_chip_list,
                                        key=lambda x: float(x[1]), reverse=True)
    sorted_known_l_chip_list_B = sorted(known_l_chip_list,
                                        key=lambda x: float(x[4]), reverse=True)
    sorted_known_l_chip_list_O = sorted(known_l_chip_list,
                                        key=lambda x: float(x[7]), reverse=True)
    sorted_known_l_chip_list_TC = sorted(known_l_chip_list,
                                        key=lambda x: float(x[10]), reverse=True)
    sorted_known_l_chip_list_BC = sorted(known_l_chip_list,
                                        key=lambda x: float(x[13]), reverse=True)
    sorted_known_l_chip_list_OC = sorted(known_l_chip_list,
                                        key=lambda x: float(x[16]), reverse=True)

    sorted_potential_l_chip_list = sorted(potential_l_chip_list,
                                          key=lambda x: float(x[1]),
                                          reverse=True)


    # write into files
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_T,
                 'known_l_chip T cell.csv')
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_B,
                 'known_l_chip B cell.csv')
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_O,
                 'known_l_chip Other.csv')
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_TC,
                 'known_l_chip T cell CNV.csv')
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_BC,
                 'known_l_chip B cell CNV.csv')
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list_OC,
                 'known_l_chip Other CNV.csv')

    write_output(potential_l_chip_list_headings + sorted_potential_l_chip_list,
                 'potential_l_chip.csv')
    write_output(remaining_genes, 'remaining_genes.csv')


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

                if row[12] == 'NS':
                    continue

                histology = row[12]
                if 'T_cell' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'point'
                    cell_type = 'T_cell'
                elif 'B_cell' in histology:
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
        cell_type_num_total_tumour_dict[cell_type] = len(tumour_list)


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
                  CNV_gene_list, gene_cell_mutation_type_dict_CNV,
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
                if 'T_cell' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'CNV'
                    cell_type = 'T_cell'
                elif 'B_cell' in histology:
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

                CNV_gene_list.append(dict_key_name)

                cell_type_num_total_tumor_dict_CNV.setdefault(cell_type,
                                                              []).append(row[4])

            # counter += 1
            # print(str(counter*100/650643) + '%')

    for cell_type, tumour_list in cell_type_num_total_tumor_dict_CNV.items():
        cell_type_num_total_tumor_dict_CNV[cell_type] = len(tumour_list)


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
        if len(info['study id']) >= 2:
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

    type_cell_mutation_list = ['T_cell point', 'B_cell point', 'Other point',
                               'T_cell CNV', 'B_cell CNV', 'Other CNV']
    for gene, cell_type_info_dict in gene_dict.items():
        l_chip_sublist = [gene]
        # each cell_type_info is a dict

        for type_cell_mutation in type_cell_mutation_list:
            # if this combination of cell type and mutation type exist
            if type_cell_mutation in cell_type_info_dict:
                info = cell_type_info_dict[type_cell_mutation]
                l_chip_sublist.append(info['frequency percentage'])
                l_chip_sublist.append(info['num_tumour_gene_cell_type'])
                l_chip_sublist.append(info['num_tumour_total_cell_type'])
            else:
                l_chip_sublist.append(0)
                l_chip_sublist.append(0)
                l_chip_sublist.append('n/a')

        if 'T_cell CNV' not in cell_type_info_dict \
                and 'B_cell CNV' not in cell_type_info_dict \
                and 'Other CNV' not in cell_type_info_dict:
            l_chip_sublist.append('___')
        else:
            l_chip_sublist.append('CNV')

            # TODO other stuff we want in the table

        l_chip_list.append(l_chip_sublist)


def write_output(gene_set, file_name):
    with open(file_name, mode='w') as genes_file:
        genes_file_writer = csv.writer(genes_file, delimiter='\t',
                                       quotechar='"',
                                       quoting=csv.QUOTE_MINIMAL)
        genes_file_writer.writerows(gene_set)


def which_genes_were_mutated_in_healthy(l_chip_gene_set_1_all, healthy_mutation_file):

    gene_list = []
    with open(healthy_mutation_file) as healthy_mutation_file:
        lines = healthy_mutation_file.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
        gene_list.extend(stripped_lines)
    l_chip_gene_set_1_healthy_mutation = set(gene_list)

    print(l_chip_gene_set_1_healthy_mutation)

    print(len(l_chip_gene_set_1_healthy_mutation))

    return l_chip_gene_set_1_all.intersection(l_chip_gene_set_1_healthy_mutation)



if __name__ == "__main__":
    # mutation file, CNV file, and gene list file
    rank_cosmic_rows(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
