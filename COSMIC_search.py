import csv
import sys


# how about, I take the cosmic files, calculate the frequency, then remove CNV
# separate into 2 tsv outputs, one of m-chip vs l-chip paper,
# one for every other gene


def rank_cosmic_rows(cosmic_mutation_file_name, cosmic_CNV_file_name,
                     gene_list_file):
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

    lymphoid_neoplasm_gene_cell_type_dict = dict()
    histology_subtype_dict = {}

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

    num_mutant_sample_total_histology_dict = {}

    # read mutation, group by gene+histology (in gene_cell_type_dict)
    # and count how many sample each histology has (in histology_subtype_dict)
    read_mutation_file(cosmic_mutation_file_name, histology_subtype_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       lymphoid_neoplasm_gene_cell_type_dict)

    print('gene-cell type pair', len(lymphoid_neoplasm_gene_cell_type_dict))

    # and store in gene_cell_type_dict
    calculate_mutation_frequency_histology_subtype(histology_subtype_dict,
                                                   lymphoid_neoplasm_gene_cell_type_dict,
                                                   num_mutant_sample_total_histology_dict)

    # filter out those not found in at least 2 studies
    reproducible = {}
    filter_non_reproducible(lymphoid_neoplasm_gene_cell_type_dict, reproducible)
    print('reproducible', len(reproducible))

    # read CNV
    lymphoid_neoplasm_gene_cell_type_dict_CNV = {}
    lymphoid_neoplasm_CNV_gene = []
    histology_subtype_dict_CNV = {}
    num_mutant_sample_total_histology_dict_CNV = {}
    read_CNV_file(cosmic_CNV_file_name, histology_subtype_dict_CNV,
                  lymphoid_neoplasm_CNV_gene,
                  lymphoid_neoplasm_gene_cell_type_dict_CNV,
                  important_column_heading_list_CNV,
                  important_column_number_list_CNV)

    # and store in lymphoid_neoplasm_gene_cell_type_dict_CNV
    calculate_CNV_frequency(histology_subtype_dict_CNV,
                            lymphoid_neoplasm_gene_cell_type_dict_CNV,
                            num_mutant_sample_total_histology_dict_CNV)

    lymphoid_neoplasm_CNV_gene_set = set(lymphoid_neoplasm_CNV_gene)

    # add CNV info into it
    # TODO can there be CNV that have more than 1 study but still
    # does not have more than 1 in point?
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
            this_gene_info = lymphoid_neoplasm_gene_cell_type_dict_CNV[
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

    # TODO extract into functions
    # gene is known, but not found on COSMIC
    remaining_genes = []
    # TODO extract gene name from them

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
        'gene', 'T Cell Frequency Point', 'B Cell Frequency Point',
        'Other Cell Frequency Point', 'T Cell Frequency CNV',
        'B Cell Frequency CNV', 'Other Cell Frequency CNV',
    ]]

    potential_l_chip_list_headings = [[
        'gene', 'T Cell Frequency Point', 'B Cell Frequency Point',
        'Other Cell Frequency Point', 'T Cell Frequency CNV',
        'B Cell Frequency CNV', 'Other Cell Frequency CNV',
    ]]

    known_l_chip_list = []
    convert_dict_list(known_l_chip_dict,
                      known_l_chip_list)
    potential_l_chip_list = []
    convert_dict_list(potential_l_chip_dict,
                      potential_l_chip_list)

    sorted_known_l_chip_list = sorted(known_l_chip_list,
                                      key=lambda x: float(x[-1]), reverse=True)
    sorted_potential_l_chip_list = sorted(potential_l_chip_list,
                                          key=lambda x: float(x[-1]),
                                          reverse=True)

    # also need to figure out how to sort CNV, after calculating its frequency

    # write into files
    write_output(known_l_chip_list_headings + sorted_known_l_chip_list,
                 'known_l_chip.csv')
    write_output(potential_l_chip_list_headings + sorted_potential_l_chip_list,
                 'potential_l_chip.csv')
    write_output(remaining_genes, 'remaining_genes.csv')


def calculate_CNV_frequency(histology_subtype_dict_CNV,
                            lymphoid_neoplasm_CNV_dict,
                            num_mutant_sample_total_histology_dict_CNV):
    # histology here is histology subtype 1
    for histology, sample_list in histology_subtype_dict_CNV.items():
        num_mutant_sample_total_histology_dict_CNV[histology] = len(
            set(sample_list))

    num_tumour_total_t_cell_b_cell_dict_CNV = {'T_cell': 0, 'B_cell': 0,
                                               'Other': 0}
    for histology, num_tumor in num_mutant_sample_total_histology_dict_CNV.items():

        if 'T_cell' in histology:
            num_tumour_total_t_cell_b_cell_dict_CNV['T_cell'] += num_tumor
        elif 'B_cell' in histology:
            num_tumour_total_t_cell_b_cell_dict_CNV['B_cell'] += num_tumor
        else:
            num_tumour_total_t_cell_b_cell_dict_CNV['Other'] += num_tumor

        # TODO there may be cancer subtypes that does not have T_cell or B_cell
        # in their name, but are indeed of those cell types

    # calculate CNV frequency
    for gene_cell_type, info in lymphoid_neoplasm_CNV_dict.items():
        genn_cell_type_list = gene_cell_type.split(';')
        cell_type = genn_cell_type_list[1]

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_sample_gene_mutant = len(set(info['id tumour'])) * 100
        # total number of mutant sample for this histology subtype 1
        num_mutant_sample_total_cell_type = \
            num_tumour_total_t_cell_b_cell_dict_CNV[cell_type]

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_sample_gene_mutant / num_mutant_sample_total_cell_type
        info['frequency percentage'] = frequency_percentage


def read_CNV_file(cosmic_CNV_file_name, histology_subtype_dict_CNV,
                  lymphoid_neoplasm_CNV_gene, gene_histology_dict_CNV,
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

                histology = row[10]
                if 'T_cell' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'CNV'
                elif 'B_cell' in histology:
                    dict_key_name = gene_name + ';' + 'B_cell' + ';' + 'CNV'
                else:
                    dict_key_name = gene_name + ';' + 'Other' + ';' + 'CNV'

                specific_gene_histology_dict_CNV = gene_histology_dict_CNV.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list_CNV)):
                    column_heading = important_column_heading_list_CNV[i]
                    column_num = important_column_number_list_CNV[i]
                    specific_gene_histology_dict_CNV.setdefault(
                        column_heading, []).append(row[column_num])

                lymphoid_neoplasm_CNV_gene.append(dict_key_name)
                # TODO could use dict_key_name instead of row 10
                histology_subtype_dict_CNV.setdefault(row[10], []).append(
                    row[4])

            # counter += 1
            # print(str(counter*100/650643) + '%')


def filter_non_reproducible(lymphoid_neoplasm_gene_histology_dict,
                            reproducible):
    for gene_histology, info in lymphoid_neoplasm_gene_histology_dict.items():
        if len(info['study id']) >= 2:
            reproducible[
                gene_histology] = info


def calculate_mutation_frequency_histology_subtype(histology_subtype_dict,
                                                   lymphoid_neoplasm_gene_cell_type_dict,
                                                   num_mutant_sample_total_histology_dict):
    # histology here is histology subtype 1
    # histology_subtype_dict map a histology subtype to a list of sample
    # num_mutant_sample_total_histology_dict counts the number of unique sample
    # for that histology
    for histology, sample_list in histology_subtype_dict.items():
        num_mutant_sample_total_histology_dict[histology] = len(
            set(sample_list))

    num_tumour_total_t_cell_b_cell_dict = {'T_cell': 0, 'B_cell': 0, 'Other': 0}
    for histology, num_tumor in num_mutant_sample_total_histology_dict.items():

        if 'T_cell' in histology:
            num_tumour_total_t_cell_b_cell_dict['T_cell'] += num_tumor
        elif 'B_cell' in histology:
            num_tumour_total_t_cell_b_cell_dict['B_cell'] += num_tumor
        else:
            num_tumour_total_t_cell_b_cell_dict['Other'] += num_tumor

        # TODO there may be cancer subtypes that does not have T_cell or B_cell
        # in their name, but are indeed of those cell types

    # calculate mutation frequency
    for gene_cell_type, info in lymphoid_neoplasm_gene_cell_type_dict.items():
        gene_cell_type_list = gene_cell_type.split(';')
        cell_type = gene_cell_type_list[1]

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_sample_gene_mutant = len(set(info['id tumour'])) * 100
        # total number of mutant sample for this histology subtype 1
        num_mutant_sample_total_cell_type = \
            num_tumour_total_t_cell_b_cell_dict[cell_type]

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_sample_gene_mutant / num_mutant_sample_total_cell_type
        info['frequency percentage'] = frequency_percentage


def read_mutation_file(cosmic_mutation_file_name, histology_subtype_dict,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_cell_type_dict):
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

                histology = row[12]
                if 'T_cell' in histology:
                    dict_key_name = gene_name + ';' + 'T_cell' + ';' + 'point'
                elif 'B_cell' in histology:
                    dict_key_name = gene_name + ';' + 'B_cell' + ';' + 'point'
                else:
                    dict_key_name = gene_name + ';' + 'Other' + ';' + 'point'

                specific_gene_cell_type_dict = gene_cell_type_dict.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list)):
                    column_heading = important_column_heading_list[i]
                    column_num = important_column_number_list[i]
                    specific_gene_cell_type_dict.setdefault(
                        column_heading, []).append(row[column_num])
                    # for example
                    # specific_gene_histology_dict.setdefault('sample name', []).append(row[4])

                # TODO: could change this into just dict_key_name
                histology_subtype_dict.setdefault(row[12], []).append(row[6])

            # counter += 1
            # print(str(counter*100/3544360) + '%')


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
            if type_cell_mutation in cell_type_info_dict:
                info = cell_type_info_dict[type_cell_mutation]
                l_chip_sublist.append(info['frequency percentage'])
            else:
                l_chip_sublist.append(0)

            # TODO other stuff we want in the table

            l_chip_list.append(l_chip_sublist)


def write_output(gene_set, file_name):
    with open(file_name, mode='w') as genes_file:
        genes_file_writer = csv.writer(genes_file, delimiter='\t',
                                       quotechar='"',
                                       quoting=csv.QUOTE_MINIMAL)
        genes_file_writer.writerows(gene_set)


if __name__ == "__main__":
    # mutation file, CNV file, and gene list file
    rank_cosmic_rows(sys.argv[1], sys.argv[2], sys.argv[3])
