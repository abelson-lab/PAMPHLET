import csv
import re
from typing import List, Dict, Tuple


def user_chose_options() -> Tuple[bool, int, int, int, int, bool]:
    """
    Ask the user to choose their options, from do they want to remove intronic mutation,
    what is their definition of recurrent mutation, targeting window size,
    indel filter threshold, cumulative contribution threshold, and merging
    :return: These fives choices
    """

    remove_intronic_mutation = get_bool(
        "Do you want remove intronic mutations, enter yes if you want to, default is yes ")
    recurrent_definition = input("How do you define recurrent mutation,"
                                 " default is a mutation is only consider"
                                 " recurrent if it appear in more than 2"
                                 " tumours, including 2 tumours ")
    targeting_window_size = input("What is your window size, or the number"
                                  " of base pairs that can be effectively"
                                  " targeted, default is 80bp ")
    indel_filter_threshold = input("What size of indel is too big to be cover"
                                   " by your window size, default is 30bp or"
                                   "greater will be filtered out ")
    cumulative_contribution_threshold = input("what is the coverage of all"
                                              " recurrent tumours you want for"
                                              " the probes, default is reporting"
                                              " until it covers 80% of the"
                                              " probes ")
    merge_other = get_bool(
        "Do you want to merge other mutation within the targeting window"
        "or do you want to have one probe target only one mutation,"
        "enter yes if you want to, default is yes ")

    if remove_intronic_mutation == "yes":
        remove_intronic_mutation = True
    else:
        remove_intronic_mutation = False

    if merge_other == "yes":
        merge_other = True
    else:
        merge_other = False

    recurrent_definition = int(recurrent_definition)
    targeting_window_size = int(targeting_window_size)
    indel_filter_threshold = int(indel_filter_threshold)
    cumulative_contribution_threshold = int(cumulative_contribution_threshold)

    return remove_intronic_mutation, recurrent_definition, targeting_window_size, indel_filter_threshold, cumulative_contribution_threshold, merge_other


"""second level function"""


def get_bool(prompt: str) -> bool:
    """a simple helper function to force the user to enter Yes or No"""
    while True:
        try:
            return {"Yes": True, "No": False}[input(prompt).lower()]
        except KeyError:
            print("Invalid input please enter Yes or No")


def read_file_choose_cancer(cosmic_mutation_file_name: str,
                            use_default: bool) -> List[List[str]]:
    """
    Reads the cosmic file, then ask the user to choose the primary tissue,
    the primary histology, and the histology subtype one they want to target
    :param cosmic_mutation_file_name: The file name of the cosmic mutation file,
    it is under COSMIC Mutation on the cosmic website
    :param use_default: chose t-cell cancers
    :return: A list of cancers that the user have chosen
    """
    all_primary_tissue_set = {}
    all_primary_histology_set = {}
    all_histology_subtype_one_set = {}

    print("searching the file to choose cancer")
    if use_default:
        with open(cosmic_mutation_file_name) as mutation_file:
            csv_reader = csv.reader(mutation_file, delimiter=',')
            for row in csv_reader:
                histology = row[12]
                if 'T_cell' in histology or 'anaplastic' in histology \
                        or 'lymphomatoid_papulosis' in histology \
                        or 'post_transplant_lymphoproliferative_disorder' in histology \
                        or 'mycosis_fungoides-Sezary_syndrome' in histology:
                    all_histology_subtype_one_set[histology] = ""
        all_primary_tissue_set["haematopoietic_and_lymphoid_tissue"] = ""
        all_primary_histology_set["lymphoid_neoplasm"] = ""

        chosen_primary_tissue_set = list(set(all_primary_tissue_set))
        chosen_primary_histology_set = list(set(all_primary_histology_set))
        chosen_histology_subtype_one_set = list(
            set(all_histology_subtype_one_set))

        chosen_sets = [chosen_primary_tissue_set, chosen_primary_histology_set,
                       chosen_histology_subtype_one_set]
        return chosen_sets

    # TODO refactor since these 3 have overlap
    # TODO refactor, how can I skip first row

    print("searching all primary tissue")
    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')

        for row in csv_reader:
            if row[7] not in all_primary_tissue_set:
                all_primary_tissue_set[row[7]] = ""

    all_primary_tissue_set = list(set(all_primary_tissue_set))
    print_set_as_numbered_list(all_primary_tissue_set)
    cancer_list_indices_list = read_number_as_list_indices()
    chosen_primary_tissue_set = []
    report_element_at_indices(all_primary_tissue_set, cancer_list_indices_list,
                              chosen_primary_tissue_set)

    print("searching all primary histology")
    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')
        for row in csv_reader:
            if row[7] in chosen_primary_tissue_set \
                    and row[11] not in all_primary_histology_set:
                all_primary_histology_set[row[11]] = ""

    all_primary_histology_set = list(set(all_primary_histology_set))
    print_set_as_numbered_list(all_primary_histology_set)
    cancer_list_indices_list = read_number_as_list_indices()
    chosen_primary_histology_set = []
    report_element_at_indices(all_primary_histology_set,
                              cancer_list_indices_list,
                              chosen_primary_histology_set)

    print("searching all primary histology")
    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')
        for row in csv_reader:
            if row[7] in chosen_primary_tissue_set \
                    and row[11] in chosen_primary_histology_set \
                    and row[12] not in all_histology_subtype_one_set:
                all_histology_subtype_one_set[row[12]] = ""

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


"""3rd level function"""


def print_set_as_numbered_list(a_list: List[str]):
    i = 1
    for cancer in a_list:
        print(i, '. ', cancer)
        i += 1


"""3rd level function"""


def read_number_as_list_indices() -> List[int]:
    cancer_list_indices = input(
        "choose primary tissue/primary histology/histology subtype 1"
        " you want to target from the list above, "
        "separate by comma (type in a number) ")

    cancer_list_indices_list = cancer_list_indices.split(',')

    cancer_list_indices_list_int = []
    for index in cancer_list_indices_list:
        cancer_list_indices_list_int.append(int(index))

    for i in range(len(cancer_list_indices_list_int)):
        cancer_list_indices_list_int[i] -= 1

    return cancer_list_indices_list_int


"""3rd level function"""


def report_element_at_indices(cancer_list: List[str],
                              cancer_list_indices_list: List[int],
                              chosen_cancer_list: List[str]):
    for indices in cancer_list_indices_list:
        chosen_cancer_list.append(cancer_list[indices])
    print(chosen_cancer_list)


def read_selected_genes(gene_list_file_name):
    """
    Read the l-chip genes that was targeted in that l-chip vs m-chip paper,
    this will be written on the output Excel file too
    """
    with open(gene_list_file_name) as gene_list_file_name:
        lines = gene_list_file_name.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
    l_chip_gene_set_1 = set(stripped_lines)
    return l_chip_gene_set_1


# TODO I can remove this, sort of, each columns is named actually,
def define_important_columns() -> Tuple[
    List[str], List[str], List[int], List[int]]:
    """
    Cosmic mutation files are organized as a table, so when reading them,
    which columns are to be read needs to be known beforehand, since the columns
    themselves are not named.
    Do note that if cosmic changes the order of columns, this function will break
    """
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


def read_process_file_point_mutation(cosmic_mutation_file_name: str,
                                     gene_mutation_type_info_dict: Dict[
                                         str, dict],
                                     important_column_heading_list: List[str],
                                     important_column_number_list: List[int],
                                     chosen_set: List[List[str]],
                                     remove_intronic_mutation: bool = True):
    """
    Read the mutation file, then calculate the frequency
    :param cosmic_mutation_file_name: The file name of the cosmic mutation file
    :param gene_mutation_type_info_dict: The dict that the mutation are to be read into
    :param important_column_heading_list: The headings of the columns to be read
    :param important_column_number_list: The index of the columns to be read
    :param remove_intronic_mutation: Whether the user choose to remove intronic mutations
    :param chosen_set: The group of cancers that the user choose to target
    """
    read_mutation_file(cosmic_mutation_file_name,
                       important_column_heading_list,
                       important_column_number_list,
                       gene_mutation_type_info_dict,
                       remove_intronic_mutation,
                       chosen_set)
    # and store in gene_info_dict
    calculate_mutation_frequency(gene_mutation_type_info_dict)


"""3rd level function"""


def read_mutation_file(cosmic_mutation_file_name: str,
                       important_column_heading_list: List[str],
                       important_column_number_list: List[int],
                       gene_mutation_type_dict: Dict[str, dict],
                       remove_intronic_mutation: bool,
                       chosen_set: List[List[str]]):
    """
    Read the mutation file into the dictionary
    gene_mutation_type_dict consist of dictionaries that map gene name to information
    mutation type is just point vs CNV
    :param cosmic_mutation_file_name: The file name of the cosmic mutation file
    :param important_column_heading_list: The headings of the columns to be read
    :param important_column_number_list: The index of the columns to be read
    :param gene_mutation_type_dict: The dict that the mutation are to be read into
    :param remove_intronic_mutation: Whether the user choose to remove intronic mutations
    :param chosen_set: The group of cancers that the user choose to target
    """

    chosen_primary_tissue_set = chosen_set[0]
    chosen_primary_histology_set = chosen_set[1]
    chosen_histology_subtype_one_set = chosen_set[2]

    print(chosen_primary_tissue_set)
    print(chosen_primary_histology_set)
    print(chosen_histology_subtype_one_set)

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')

        count_tumour_mutation_before_intronic(chosen_histology_subtype_one_set,
                                              chosen_primary_histology_set,
                                              chosen_primary_tissue_set,
                                              csv_reader)

    with open(cosmic_mutation_file_name) as mutation_file:
        csv_reader = csv.reader(mutation_file, delimiter=',')
        for row in csv_reader:
            # 0, 4, 6, 11, 12, 13, 14.  19, 21, 24, 25, 30
            # gene name, sample name, tumour id, primary site primary histology and its subtypes (3),
            # mutation CDS, mutation desc, GRch, mutation genome position, study id

            # if this mutation belong to a chosen cancer
            if row[7] in chosen_primary_tissue_set \
                    and row[11] in chosen_primary_histology_set \
                    and row[12] in chosen_histology_subtype_one_set:

                gene_ensembl = row[0].split('_')
                gene_name = gene_ensembl[0]

                filter_out_flag = False
                if remove_intronic_mutation:
                    filter_out_flag = filter_intronic_mutations(row[19])
                if filter_out_flag or row[12] == 'NS' or row[
                    21] == 'Substitution - coding silent':
                    continue

                # if we need to only target t-cell
                # user would input when asked what histology subtype 1 do you want
                #   dict_key_name should just be
                dict_key_name = gene_name + ';' + 'point'

                specific_gene_dict = gene_mutation_type_dict.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list)):
                    column_heading = important_column_heading_list[i]
                    column_num = important_column_number_list[i]
                    specific_gene_dict.setdefault(
                        column_heading, []).append(row[column_num])
                    # for example
                    # specific_gene_dict.setdefault('sample name', []).append(row[4])

            # counter += 1
            # print(str(counter*100/3544360) + '%')


def filter_intronic_mutations(mutation_CDS) -> bool:
    """
    # if it is intron (where the position is known) and not splice site and not
    # ending or start in an exon, and both start and end in two different introns,
    # then skip it (set filter out flag is true)
    :param mutation_CDS: The variant notation, follows HGVS nomenclature
    :return:
    """
    # assert 'c.' in mutation_CDS
    # check_variant_type(row)
    filter_out_flag = filter_known_pure_intronic_mutation(
        mutation_CDS)
    # we don't need this, there are no such mutation in this file
    # filter_out_flag_intron_range = is_intron_range_mutation(
    #     mutation_CDS)
    #
    # if filter_out_flag_intron_range:
    #     check_intronic_mutation_diff_intron(mutation_CDS)
    return filter_out_flag


def count_tumour_mutation_before_intronic(chosen_histology_subtype_one_set,
                                          chosen_primary_histology_set,
                                          chosen_primary_tissue_set,
                                          csv_reader):
    """
    Count the total number of unique tumour and unique mutations
    including intronic mutation
    :param chosen_histology_subtype_one_set: The set of histology subtype one that the user chose
    :param chosen_primary_histology_set: The set of primary histology that the user chose
    :param chosen_primary_tissue_set: The set of primary tissue that the user chose
    :param csv_reader:
    :return:
    """
    unique_tumour_before_intronic_mutation = []
    unique_mutation_before_intronic_mutation = []
    for row in csv_reader:
        if row[7] in chosen_primary_tissue_set \
                and row[11] in chosen_primary_histology_set \
                and row[12] in chosen_histology_subtype_one_set:

            unique_tumour_before_intronic_mutation.append(row[4])
            unique_mutation_before_intronic_mutation.append(row[25])

    unique_tumour_before_intronic_mutation = set(
        list(unique_tumour_before_intronic_mutation))
    unique_mutation_before_intronic_mutation = set(
        list(unique_mutation_before_intronic_mutation))
    print('total number of mutations before intronic filter',
          len(unique_mutation_before_intronic_mutation))
    print('total number of tumours before intronic filter',
          len(unique_tumour_before_intronic_mutation))


"""3rd level functions"""


def check_variant_type(mutation_cds: str):
    """
    if all variant were covered by this check then nothing should be printed
    :param mutation_cds: The variant notation, follows HGVS nomenclature
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
            mutation_cds) is None and \
            re.search("c\\.[?]", mutation_cds) is None and re.search(
        "c\\.1_\\*del",
        mutation_cds[19]) is None \
            and re.search("c\\.\\*[-]?[0-9]+(_\\*[-]?[0-9]+)?[ACGTdi>?]",
                          mutation_cds) is None \
            and re.search("c\\.\\([0-9]+_[0-9]+\\)[ACGTdi>?]",
                          mutation_cds) is None \
            and re.search(
        "c\\.(\\()?[0-9]+[+\\-]\\?(_[0-9]+[+\\-]\\?)?(\\))?[ACGTdi>?]",
        mutation_cds) is None:
        print(mutation_cds)


def calculate_mutation_frequency(gene_mutation_type_info_dict):
    all_tumour = []
    for gene_mutation_type, info in gene_mutation_type_info_dict.items():
        all_tumour.extend(info['id tumour'])
    all_unique_tumour = len(list(set(all_tumour)))

    # calculate mutation frequency
    for gene_mutation_type, info in gene_mutation_type_info_dict.items():

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_tumour_gene = len(set(info['id tumour'])) * 100
        # calculate the total number of tumour
        num_tumour_total = all_unique_tumour

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_tumour_gene / num_tumour_total
        info['frequency percentage'] = frequency_percentage
        info['num_tumour_gene'] = len(set(info['id tumour']))
        info['num_tumour_total'] = num_tumour_total


def is_intron_range_mutation(mutation_CDS: str) -> bool:
    filter_out_flag_intron_range = re.search(
        'c\\.[0-9]+([+\\-][0-9]+)_[0-9]+([+\\-][0-9]+)[ACGTdi>?]',
        mutation_CDS) is not None
    return filter_out_flag_intron_range


def filter_known_pure_intronic_mutation(mutation_CDS: str) -> bool:
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


def check_intronic_mutation_diff_intron(mutation_CDS: str):
    """
    this check if (for example c.4072-1234_5155-246del)
    the deletion is from an intron between exon 29 and exon 30
    # all the way to another intron between exon 35 and exon 36
    # actually if I can check if the "last nucleotide of the directly upstream exon" of the start of the range
    # and the "first nucleotide of the directly downstream exon" of the end of the range
    # differ by 1
    # if they do, then that deletion or duplication is contained in one intron
    # if not, it is not contained in one exon
    # if it involves a UTR (either start or end), then we will include it anyway, so don't need to check
    # since it does not match, we do not need to do anything about it
    :param mutation_CDS: The variant notation, follows HGVS nomenclature
    """

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


# TODO combine these three function with read mutation file and calculate mutation frequency

def read_process_file_CNV_mutation(cosmic_CNV_file_name,
                                   gene_cell_mutation_type_info_dict_CNV,
                                   important_column_heading_list_CNV,
                                   important_column_number_list_CNV):
    # each gene cell-type refers to the gene ABC in t-cell CNV

    cell_type_num_tumor_dict_CNV = {}
    read_CNV_file(cosmic_CNV_file_name,
                  gene_cell_mutation_type_info_dict_CNV,
                  important_column_heading_list_CNV,
                  important_column_number_list_CNV)
    print('gene cell-type pair CNV', len(gene_cell_mutation_type_info_dict_CNV))
    print('CNV cell type distribution', cell_type_num_tumor_dict_CNV)
    # and store in gene_mutation_type_info_dict_CNV
    calculate_CNV_frequency(gene_cell_mutation_type_info_dict_CNV)


def read_CNV_file(cosmic_CNV_file_name,
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
                elif 'B_cell' in histology or 'Burkitt' in histology or 'chronic' in histology \
                        or 'follicular_lymphoma' in histology or 'hairy' in histology \
                        or 'hodgkin' in histology or 'lymphoplasmacytic_lymphoma' in histology \
                        or 'mantle' in histology or 'marginal' in histology \
                        or 'plasma_cell_myeloma' in histology or 'plasmacytoma' in histology \
                        or 'effusion' in histology or 'MALT' in histology:
                    dict_key_name = gene_name + ';' + 'B_cell' + ';' + 'CNV'
                else:
                    dict_key_name = gene_name + ';' + 'Other' + ';' + 'CNV'

                specific_gene_histology_dict_CNV = gene_cell_mutation_type_dict_CNV.setdefault(
                    dict_key_name, {})

                for i in range(len(important_column_heading_list_CNV)):
                    column_heading = important_column_heading_list_CNV[i]
                    column_num = important_column_number_list_CNV[i]
                    specific_gene_histology_dict_CNV.setdefault(
                        column_heading, []).append(row[column_num])

            # counter += 1
            # print(str(counter*100/650643) + '%')


def calculate_CNV_frequency(gene_mutation_type_info_dict_CNV):
    all_tumour = []
    for gene_mutation_type_CNV, info_CNV in gene_mutation_type_info_dict_CNV.items():
        all_tumour.extend(info_CNV['id tumour'])
    all_unique_tumour = len(list(set(all_tumour)))

    # calculate mutation frequency
    for gene_mutation_type, info in gene_mutation_type_info_dict_CNV.items():

        # number of sample this gene is found to have a mutation
        # multiply by 100 to turn into percentage
        num_tumour_gene = len(set(info['id tumour'])) * 100
        # total number of mutant sample for this histology subtype 1
        num_tumour_total = all_unique_tumour

        # this sort of calculates what percentage of people who have this
        # type of cancer has this gene are mutant
        frequency_percentage = num_tumour_gene / num_tumour_total
        info['frequency percentage'] = frequency_percentage
        info['num_tumour_gene_cell_type'] = len(set(info['id tumour']))
        info['num_tumour_total_cell_type'] = num_tumour_total
