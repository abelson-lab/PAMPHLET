import copy

import matplotlib.pyplot as plt
from typing import Any, List, Dict, Set, Tuple, Union

from others import write_output_excel, write_output_txt

"""the only 0th level function in this file"""


def choose_probe_placement_point(
        gene_mutation_type_info_dict: Dict[str, dict],
        targeting_window_size: int = 80,
        cumulative_contribution_threshold: int = 90,
        indel_filter_threshold: int = 30, recurrent_definition: int = 2,
        merge_others: bool = True):
    """
    Choose probe using one of two approach, generate all probe then greedily choose
    the most recurrent one (the set cover greedy algorithm approach), \
    rank all mutations based on recurrent-ness then add all mutation within the
    targeting window centered around that mutation
    :param gene_mutation_type_info_dict: A dictionary mapping gene to info
    :param targeting_window_size: see user_chose_options
    :param cumulative_contribution_threshold: see user_chose_options
    :param indel_filter_threshold: see user_chose_options
    :param recurrent_definition: see user_chose_options
    :param merge_others: see user_chose_options
    """

    # abandoned approach
    # mutation_list = rank_all_probe_choose_most_recurrent(
    #     gene_mutation_type_info_dict, cumulative_contribution_threshold,
    #     targeting_window_size,
    #     indel_filter_threshold, recurrent_definition)

    mutation_list, indel_probe_list = chose_most_recurrent_mutation_then_probe_centered_at_mutation_center(
        gene_mutation_type_info_dict, cumulative_contribution_threshold,
        targeting_window_size,
        indel_filter_threshold, recurrent_definition, merge_others)

    mutation_list_IGV = copy.deepcopy(mutation_list)
    mutation_list_IGV_format = [probe_range[0:9] for probe_range in mutation_list_IGV]
    for probe_range in mutation_list_IGV_format:
        probe_range[1] = str(int(probe_range[1]) - 1)
        probe_range[6] = str(int(probe_range[6]) - 1)
    write_output_txt(mutation_list_IGV_format, 'range_IGV.xlsx')


    mutation_list.insert(0, ['chromosome', 'first_mutation',
                             'last mutation', 'gene name and range number ',
                             'just 1', 'just plus', 'first mutation',
                             'last mutation', 'color code',
                             'cumulative contribution', 'range mutations',
                             'tumours covered', 'number of tumours'])
    indel_probe_list.insert(0, ['gene name and range number', 'chromosome',
                                'first mutation', 'last mutation', 'range mutations'])
    write_output_excel(mutation_list, 'range.xlsx')
    write_output_excel(indel_probe_list, 'indels_ranges.xlsx')


"""1st level function"""
"""the choose most recurrent mutation and centered approach"""


def chose_most_recurrent_mutation_then_probe_centered_at_mutation_center(
        gene_info_dict: Dict[str, dict], cumulative_contribution_threshold: int,
        targeting_window_size: int, indel_filter_threshold: int,
        recurrent_definition: int, merge_others: bool) -> Tuple[
            List[Any], List[Any]]:
    """
    This approach choose the most recurrent mutation based on number of
    tumour that has this mutation (the greater the number of tumour, the more
    recurrent it is), and then these covered tumour are removed from every other
    mutation's list, so that all other mutation can be ranked by the number of
    tumours they have in addition to the ones already covered. Due to that, the already selected
    mutation will have number of tumours of zero. Additionally, for each selected
    mutation, it will attempt to see if any other mutation can be covered
    by the targeting window whose center is the center of the selected mutation.
    Repeat until total number of tumour covered by the selected mutation and
    its targeting window exceeds the cumulative contribution threshold, which
    is a percentage of the total number of unique recurrent tumours

    :param gene_info_dict: A dictionary mapping gene to info
    :param cumulative_contribution_threshold: see user_chose_options
    :param targeting_window_size: see user_chose_options
    :param indel_filter_threshold: see user_chose_options
    :param recurrent_definition: see user_chose_options
    :param merge_others: see user_chose_options
    :return: A list of probe, with two genomic location, indicating the
    5' most mutation and the 3' most mutation of that probe, and the cumulative
    contribution up to and including that probe.
    """
    all_tumour_id = []
    gene_tumour_position = {}

    # group tumour id based on genomic position
    get_mapping_position_tumour(all_tumour_id, gene_tumour_position,
                                gene_info_dict,
                                indel_filter_threshold)

    # each gene has one dict holding position mapping to tumour
    # and another mapping position with deletion to tumour
    # combine them two
    # note that position_tumour already has insertions
    all_position_tumour_dict = {}
    all_position_tumour_dict_substitution = {}
    all_position_tumour_dict_deletion = {}
    all_position_tumour_dict_insertion = {}
    position_gene_dict = {}
    input_list = []
    for gene, position_tumour_dict_list in gene_tumour_position.items():

        # gene_tumour_position[gene_cell_type] = [
        #     position_tumour, position_tumour_substitution, position_tumour_deletion,
        #     position_tumour_insertion, position_mutation]

        position_tumour = position_tumour_dict_list[0]
        position_tumour_substitution = position_tumour_dict_list[1]
        position_tumour_deletion = position_tumour_dict_list[2]
        position_tumour_insertion = position_tumour_dict_list[3]
        position_mutation = position_tumour_dict_list[4]

        for position, tumour_set in position_tumour.items():
            # add every mutation into a these two dict
            position_gene_dict[position] = gene
            all_position_tumour_dict[position] = set(tumour_set)

            if position in position_mutation:
                mutation_notation = set(position_mutation[position])
                unique_tumours = set(tumour_set)
                num_unique_tumours = len(list(set(tumour_set)))
                input_list.append(
                    [position, mutation_notation, unique_tumours,
                     num_unique_tumours, gene])

        for substitution, tumour_set in position_tumour_substitution.items():
            all_position_tumour_dict_substitution[substitution] = set(tumour_set)

        for deletion, tumour_set in position_tumour_deletion.items():
            all_position_tumour_dict_deletion[deletion] = set(tumour_set)

        for insertion, tumour_set in position_tumour_insertion.items():
            all_position_tumour_dict_insertion[insertion] = set(tumour_set)

    # sort by the number of unique tumours
    input_list.sort(key=lambda x: x[3], reverse=True)

    input_list.insert(0, ['mutation_location', 'mutation aa', 'tumour set',
                          'size of tumour set', 'gene name'])

    # show data before choose ranges
    write_output_excel(input_list, 'input.xlsx')

    tumour_set_list = []
    unique_tumour_id = []
    recurrent_unique_tumour_id = count_mutation_and_tumour_after_filters(
        gene_tumour_position, recurrent_definition, tumour_set_list,
        unique_tumour_id)

    # center_probe_dict = all_position_tumour_dict.copy()

    print('start finding cover')
    covered_tumour_id = []
    cover_sizes_percentage = [0]
    percentage_cover_threshold = len(
        recurrent_unique_tumour_id) * cumulative_contribution_threshold / 100
    mutation_list = []
    indel_probe_list = []
    find_probe_cover(all_position_tumour_dict, cover_sizes_percentage,
                     covered_tumour_id, mutation_list,
                     percentage_cover_threshold, recurrent_definition,
                     recurrent_unique_tumour_id, targeting_window_size,
                     all_position_tumour_dict_deletion, merge_others,
                     all_position_tumour_dict_insertion, indel_probe_list,
                     position_gene_dict)

    # print('start centering probes')
    # final_mutation_list = []
    # centralize_each_probe(center_probe_dict, final_mutation_list, mutation_list,
    #                       recurrent_definition, targeting_window_size)

    plt.plot(cover_sizes_percentage, linewidth=1)
    plt.ylabel('percentage of tumour covered')
    plt.xlabel('number of ' + str(targeting_window_size) + 'bp ranges selected')
    plt.title('coverage of tumours with increasing probe in myeloid')
    plt.savefig('stop_at_90.pdf')

    return mutation_list, indel_probe_list


"""2nd level function"""


def get_mapping_position_tumour(all_tumour_id: List[str],
                                gene_tumour_position: Dict[str, List[
                                    Dict[str, List[str]]]],
                                gene_mutation_type_info_dict: Dict[str, dict],
                                indel_filter_threshold: int):
    """
    Turn gene_mutation_type_info_dict's info which for every genes, there are and
    are multiple list and each element corresponds to one mutation, we want to make
    a new dictionary that maps position to tumour set.
    Also, we want to know how many mutation and tumour are there before after
    filtering indels
    :param all_tumour_id: All tumour id present in the file
    :param gene_tumour_position: A dictionary mapping gene to a list of two
     other dictionaries, position tumour (a dict mapping position to tumour) and
     position_range_tumour_deletion (a dict mapping position_range to tumour set
     for deletion mutations)
    :param gene_mutation_type_info_dict: maps gene name to info
    :param indel_filter_threshold: indels greater than this will be filter out
    :return:
    """

    all_tumour_id_before_indel_filter = []
    num_mutation_before_indel_filter = 0
    for gene_cell_type, info in gene_mutation_type_info_dict.items():

        for i in range(len(info['mutation genome position'])):
            tumour_id = info['id tumour'][i]
            all_tumour_id_before_indel_filter.append(tumour_id)

        num_mutation_before_indel_filter += len(
            list(set(info['mutation genome position'])))

    all_tumour_id_before_indel_filter = list(
        set(all_tumour_id_before_indel_filter))

    print("total number of mutation after intronic before indel filter",
          num_mutation_before_indel_filter)
    print("total number of tumours after intronic before indel filter",
          len(all_tumour_id_before_indel_filter))

    for gene_cell_type, info in gene_mutation_type_info_dict.items():
        # assert len(info['mutation genome position']) == len(info['id tumour'])

        position_tumour = {}
        position_tumour_substitution = {}
        position_tumour_deletion = {}
        position_tumour_insertion = {}
        position_mutation = {}
        gene_tumour_position[gene_cell_type] = [
            position_tumour, position_tumour_substitution, position_tumour_deletion,
            position_tumour_insertion, position_mutation]

        for i in range(len(info['mutation genome position'])):
            chromosome_position_range = info['mutation genome position'][i]
            tumour_id = info['id tumour'][i]
            mutation_aa = info['mutation AA'][i]
            if chromosome_position_range != 'null' and tumour_id != 'null':

                # all genomic position are in the format a:b-c
                # where an is the chromosome number, followed by a semicolon
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

                all_tumour_id.append(tumour_id)

                position_tumour.setdefault(chromosome_position_range, []).append(tumour_id)
                position_mutation.setdefault(chromosome_position_range, []).append(mutation_aa)

                # deletion
                if int(position_range_end) > int(position_range_start) + 1:
                    position_tumour_deletion.setdefault(chromosome_position_range, []).append(tumour_id)
                # insertion
                elif int(position_range_end) == int(position_range_start) + 1:
                    position_tumour_insertion.setdefault(chromosome_position_range, []).append(tumour_id)
                # substitution
                else:
                    position_tumour_substitution.setdefault(chromosome_position_range, []).append(tumour_id)


def count_mutation_and_tumour_after_filters(
        gene_tumour_position: Dict[str, List[Dict[str, List[str]]]],
        recurrent_definition: int,
        tumour_set_list: List[List[str]], unique_tumour_id: List[str]):
    """
    # extract all sets of tumour, remove duplicates within each set, remove
    # duplicate set. that is convert to set for all 'subset', and convert to
    # multiset for all the collection of the subset remove duplicates and count
    # number of unique tumours this is to know the number of (recurrent) unique tumours
    :param gene_tumour_position: A dictionary mapping gene to a list of two
     other dictionaries, position tumour (a dict mapping position to tumour) and
     position_range_tumour_deletion (a dict mapping position_range to tumour set
     for deletion mutations)
    :param recurrent_definition: see user_chose_option
    :param tumour_set_list: A list of all tumour sets (so length of list is
     number of mutations)
    :param unique_tumour_id: A list of all unique tumour ids
    :return: All recurrent unique tumour id
    """
    mutations_locations = []
    for gene, position_tumour_dict_list in gene_tumour_position.items():
        position_tumour = position_tumour_dict_list[0]
        position_range_tumour_deletion = position_tumour_dict_list[2]
        position_range_tumour_insertion = position_tumour_dict_list[3]

        # add tumours ids from deletion mutation
        for deletion, tumour_ids in position_range_tumour_deletion.items():
            tumour_ids_set = list(set(tumour_ids))
            tumour_set_list.append(tumour_ids_set)
            unique_tumour_id.extend(tumour_ids_set)
            mutations_locations.append(deletion)

        # add tumours ids from other mutations
        for position, tumour_ids in position_tumour.items():
            tumour_ids_set = list(set(tumour_ids))
            tumour_set_list.append(tumour_ids_set)
            unique_tumour_id.extend(tumour_ids_set)

            position_range_dict = {}
            make_dict_chromosome_position_to_range(
                position_range_tumour_insertion, position_range_dict)

            # if this position does not belong to an insertion, then count it
            if position_range_dict.get(position) is None:
                mutations_locations.append(position)

        for insertion, tumour_ids in position_range_tumour_insertion.items():
            mutations_locations.append(insertion)

    unique_tumour_id = list(set(unique_tumour_id))
    # we now have
    # a dict with genomic position to a list of tumour ids and its reverse
    # a list of set of tumour ids
    # recurrent mutation here is a base pair that is mutated in multiple tumour
    # since each tumour set is a set a tumour that is mutated at that position
    # keep only recurrent mutations, means only using sets with >=2 tumour
    # which is the default recurrent definition
    recurrent_unique_tumour_set_list = [tumour_set
                                        for tumour_set in tumour_set_list
                                        if
                                        len(tumour_set) >= recurrent_definition]
    # then only keep these tumour ids
    recurrent_unique_tumour_id = list(set([tumour_id
                                           for tumour_set in
                                           recurrent_unique_tumour_set_list
                                           for tumour_id in tumour_set]))

    unique_mutations_locations = list(set(mutations_locations))
    print('total number of mutations after indel filter before recurrent',
          len(unique_mutations_locations))
    print('total number of tumours after indel filter before recurrent',
          len(unique_tumour_id))
    print('total number of recurrent mutations after indel filter',
          len(recurrent_unique_tumour_set_list))
    print(
        'total number of tumours with only recurrent mutations after indel filter',
        len(recurrent_unique_tumour_id))
    return recurrent_unique_tumour_id



def find_probe_cover(all_position_tumour_dict: Dict[str, Set[str]],
                     cover_sizes_percentage: List[float],
                     covered_tumour_id: List[str],
                     mutation_list: List[
                         Union[str, float, int, List[Union[str, List[str]]]]],
                     percentage_cover_threshold: float,
                     recurrent_definition: int,
                     recurrent_unique_tumour_id: List[str],
                     targeting_window_size: int,
                     all_position_tumour_dict_deletion: Dict[str, Set[str]],
                     merge_others: bool,
                     all_position_tumour_dict_insertion: Dict[str, Set[str]],
                     indels_probe_list: List[List[str]],
                     position_gene_dict: Dict[str, str]):
    """
    As long as it has not reached the cumulative contribution threshold,
    find position the currently has most tumour ids, then search in a targeting
    window centered around that position,
    :param all_position_tumour_dict: A dict mapping position to tumour set
    :param cover_sizes_percentage: The list of all cumulative contribution of
     each probe
    :param covered_tumour_id: The set of tumour id that is covered by this
     selection of mutation/probes
    :param mutation_list: A list of mutation/probe that contain all information
     about those probes (the position it starts at or central position), the
     cumulative contribution up to that probe, leftmost mutation in the probe,
     rightmost mutation in the probe, and the tumour set that it covers
    :param percentage_cover_threshold: The covered tumour ids of
     the selected probe should only stop after exceed this threshold. It is a
     percentage of recurrent_unique_tumour_id
    :param recurrent_definition: see user_chose_options
    :param recurrent_unique_tumour_id: The tumour ids that are unique and recurrent
    :param targeting_window_size: see user_chose_options
    :param all_position_tumour_dict_deletion: same as all_position_tumour_dict,
     but only with deletions
    :param all_position_tumour_dict_insertion: same as all_position_tumour_dict,
     but only with insertions
    :param merge_others: see user_chose_options
    :param indels_probe_list: a list of probes that captures indels
    :param position_gene_dict: mapping position to gene
    :return:
    """
    num_probe = 0
    while len(covered_tumour_id) < percentage_cover_threshold:

        # find the position that currently has the most tumour ids
        num_most_tumour = 0
        position_most_tumour = ''
        for position, tumour_set in all_position_tumour_dict.items():
            tumour_list = list(tumour_set)
            if len(tumour_list) > num_most_tumour:
                num_most_tumour = len(tumour_list)
                position_most_tumour = position

        # get the 5' most and the 3' most position in the probe, and the tumour
        # set that is covered by this probe
        first_mutation_selected_position_and_probe, \
        last_mutation_selected_position_and_probe, \
        tumour_set_selected_position_and_probe, \
        chromosome, all_mutation_in_probe, capture_indels = get_probe_covered_tumours(
            all_position_tumour_dict, position_most_tumour,
            recurrent_definition, targeting_window_size,
            all_position_tumour_dict_deletion, merge_others,
            all_position_tumour_dict_insertion)

        # add the newly covered tumour ids, and keep it unique
        covered_tumour_id.extend(
            list(set(tumour_set_selected_position_and_probe)))
        covered_tumour_id = list(set(covered_tumour_id))
        # calculate new cumulative contribution
        cumulative_contribution = len(covered_tumour_id) / len(
            recurrent_unique_tumour_id)

        gene_name = position_gene_dict[position_most_tumour].split(';')[0]

        mutation_list.append([
            'chr' + chromosome,
            first_mutation_selected_position_and_probe,
            last_mutation_selected_position_and_probe,
            gene_name + ' ' + str(num_probe),
            '1',
            '+',
            first_mutation_selected_position_and_probe,
            last_mutation_selected_position_and_probe,
            '21816532',
            cumulative_contribution,
            all_mutation_in_probe,
            tumour_set_selected_position_and_probe,
            len(tumour_set_selected_position_and_probe)
        ])
        # position most tumour is represented by its central base pair,
        # and using the dictionary, this will map back to the gene
        # this assumes that no probe covers 2 genes

        # without this it does not work?
        # not sure why, supposedly it was all removed (using pop)
        if position_most_tumour in all_position_tumour_dict:
            all_position_tumour_dict.pop(position_most_tumour)

        cover_sizes_percentage.append(cumulative_contribution)
        print(cumulative_contribution, position_most_tumour)

        all_mutation_in_probe.sort()
        if capture_indels:
            indels_probe_list.append(
                [gene_name + ' ' + str(num_probe),
                 'chr' + chromosome,
                 first_mutation_selected_position_and_probe,
                 last_mutation_selected_position_and_probe,
                 list(set(all_mutation_in_probe))])

        num_probe += 1


"""3rd level function """


# TODO, do docstring for this, I forgot this last time
def get_probe_covered_tumours(all_position_tumour_dict: Dict[str, Set[str]],
                              position_most_tumour: str,
                              recurrent_definition: int,
                              targeting_window_size: int,
                              all_position_tumour_dict_deletion: Dict[
                                  str, Set[str]],
                              merge_others: bool,
                              all_position_tumour_dict_insertion: Dict[
                                  str, Set[str]]) \
        -> Tuple[int, int, List[str], str, List[str], bool]:
    """
    :param all_position_tumour_dict: A dictionary mapping mutation to tumour sets
    :param position_most_tumour: The mutation that currently has the greatest number
     of mutations
    :param recurrent_definition: see user_chose_option
    :param targeting_window_size: see user_chose_option
    :param all_position_tumour_dict_deletion: same as all_position_tumour_dict,
     but only for deletions
    :param merge_others: see user_chose_options
    :param all_position_tumour_dict_insertion: same as all_position_tumour_dict,
     but only for insertions
    :return:
    """
    # get its tumour set
    tumour_set_selected_position = all_position_tumour_dict[
        position_most_tumour]

    # get leftmost and rightmost probe, then start search within the window size,
    # this search where the mutation is still covered
    position_list, most_tumour_range_end, most_tumour_range_start, \
    semicolon_index = parse_chromosome_position_range(position_most_tumour)
    leftmost_probe_start = int(most_tumour_range_end) - targeting_window_size + 1
    rightmost_probe_end = int(most_tumour_range_start) + targeting_window_size

    first_mutation_selected_position_and_probe = int(most_tumour_range_start) - 1
    last_mutation_selected_position_and_probe = int(most_tumour_range_end) + 1
    tumour_set_selected_position_and_probe = []
    tumour_set_selected_position_and_probe.extend(tumour_set_selected_position)
    all_mutation_in_probe = [position_most_tumour]
    chromosome = position_most_tumour[:semicolon_index]
    captures_indels = False

    # TODO if not merging and the only mutation is a deletion or insertion
    #   It does not cover that
    #   if I turn insertion in a range instead (like deletions, then we can)
    #   check right here (without going into merge_others)

    # TODO: this might break if we are not merging
    if merge_others:
        best_probe_tumour_dict = {}
        best_probe_mutation_dict = {}
        best_probe_position_dict = {}
        best_probe_indel_flag_dict = {}
        probe_start = leftmost_probe_start
        probe_end = leftmost_probe_start + targeting_window_size
        print("position most tumour", position_most_tumour)
        while probe_end < rightmost_probe_end:

            tumour_set_selected_position_and_probe = []
            all_position_in_probe = []
            all_mutation_in_probe = []

            captures_indels = merge_other_mutations_in_probe(
                all_position_in_probe, all_position_tumour_dict,
                all_position_tumour_dict_deletion, chromosome, probe_start,
                probe_end, recurrent_definition, tumour_set_selected_position_and_probe,
                all_position_tumour_dict_insertion, all_mutation_in_probe)

            best_probe_indel_flag_dict[probe_start] = captures_indels
            best_probe_tumour_dict[probe_start] = tumour_set_selected_position_and_probe
            best_probe_mutation_dict[probe_start] = all_mutation_in_probe
            best_probe_position_dict[probe_start] = all_position_in_probe

            print('\t', probe_start, probe_end,
                  len(tumour_set_selected_position_and_probe), captures_indels)
            probe_start += 1
            probe_end += 1

        # find the best probe that covers the mutation at position_most_tumour
        best_range_tumour_set_size = 0
        best_range_start = 0
        for probe_start, tumour_set in best_probe_tumour_dict.items():
            if len(tumour_set) > best_range_tumour_set_size:
                best_range_start = probe_start
                best_range_tumour_set_size = len(tumour_set)
        best_probe_tumour_set = best_probe_tumour_dict.get(best_range_start)
        best_probe_mutation = best_probe_mutation_dict.get(best_range_start)
        best_probe_position = best_probe_position_dict.get(best_range_start)
        best_probe_capture_indel = best_probe_indel_flag_dict.get(best_range_start)
        assert best_probe_tumour_set is not None
        assert best_probe_mutation is not None
        assert best_probe_position is not None

        # now remove the mutations of the best range
        for mutation in best_probe_mutation:
            if mutation in all_position_tumour_dict:
                all_position_tumour_dict.pop(mutation)
            if mutation in all_position_tumour_dict_deletion:
                all_position_tumour_dict_deletion.pop(mutation)
            if mutation in all_position_tumour_dict_insertion:
                all_position_tumour_dict_insertion.pop(mutation)

        best_probe_mutation.sort()
        best_probe_position.sort()
        # now one before
        first_mutation_selected_position_and_probe = best_probe_position[0] - 1
        last_mutation_selected_position_and_probe = best_probe_position[-1] + 1
        all_mutation_in_probe = best_probe_mutation
        tumour_set_selected_position_and_probe = best_probe_tumour_set

        captures_indels = best_probe_capture_indel

        print('best outer range', best_range_start,
              best_range_start + targeting_window_size)
        print("best range", first_mutation_selected_position_and_probe,
              last_mutation_selected_position_and_probe)

    return first_mutation_selected_position_and_probe, last_mutation_selected_position_and_probe, \
        tumour_set_selected_position_and_probe, chromosome, all_mutation_in_probe, captures_indels


"""4th level function"""


def merge_other_mutations_in_probe(all_position_in_probe,
                                   all_position_tumour_dict,
                                   all_position_tumour_dict_deletion,
                                   chromosome, probe_start, probe_end,
                                   recurrent_definition,
                                   tumour_set_selected_position_and_probe,
                                   all_position_tumour_dict_insertion,
                                   all_mutation_in_probe: List[str]) -> bool:
    """
    :param all_position_in_probe:
    :param all_position_tumour_dict:
    :param all_position_tumour_dict_deletion:
    :param chromosome:
    :param probe_end:
    :param probe_start:
    :param recurrent_definition:
    :param tumour_set_selected_position_and_probe:
    :param all_position_tumour_dict_insertion:
    :param all_mutation_in_probe: Every recurrent mutation in the probe,
     in range form
    :return: captures_indels: if this probes captures indels
    """
    # TODO docstring and type hits

    captures_indels = False
    # for every position from probe start to probe end, which is centered
    # around the center of the chosen recurrent position
    for i in range(probe_start, probe_end):

        # combine position with chromosome to get genomic location
        a_position_in_probe = i
        position_range_in_genome = chromosome + ':' + str(a_position_in_probe) + '-' + str(a_position_in_probe)
        position_in_genome = chromosome + ':' + str(a_position_in_probe)

        # get it from non-deletions
        a_position_in_probe_tumour_set = all_position_tumour_dict.get(
            position_range_in_genome)

        # if a position is mutated, and is recurrent
        if a_position_in_probe_tumour_set is not None and len(
                a_position_in_probe_tumour_set) >= recurrent_definition:

            # add its tumour set to the tumour set that probe covers
            tumour_set_selected_position_and_probe.extend(
                a_position_in_probe_tumour_set)
            all_position_in_probe.append(a_position_in_probe)
            all_mutation_in_probe.append(position_range_in_genome)

        # TODO: could move making dict outside, I doing the same thing multiple times
        # make a dict mapping from position to a deletion range
        chromosome_position_deletion_dict = {}
        make_dict_chromosome_position_to_range(
            all_position_tumour_dict_deletion,
            chromosome_position_deletion_dict)

        # see which deletion does this position belong to
        the_deletion_range_list = chromosome_position_deletion_dict.get(
            position_in_genome)
        if the_deletion_range_list is not None:
            for a_deletion_range in the_deletion_range_list:
                # get the deletion's tumour set
                a_position_in_probe_tumour_set_deletion = all_position_tumour_dict_deletion.get(
                    a_deletion_range)
                # check again if it is mutated and if it is recurrent
                if a_position_in_probe_tumour_set_deletion is not None and len(
                        a_position_in_probe_tumour_set_deletion) >= recurrent_definition:
                    position_list, deletion_range_end, deletion_range_start, semicolon_index \
                        = parse_chromosome_position_range(a_deletion_range)
                    # if the entire deletion can be covered by the probe
                    if deletion_range_end < str(probe_end) and deletion_range_start > str(probe_start):
                        # if this deletion is not already covered
                        if a_deletion_range not in all_mutation_in_probe:
                            # add its tumour set to the overall tumour set that probe covers
                            tumour_set_selected_position_and_probe.extend(
                                a_position_in_probe_tumour_set_deletion)
                            # then add the start and end of that deletion to all positions
                            all_position_in_probe.append(int(deletion_range_start))
                            all_position_in_probe.append(int(deletion_range_end))
                            all_mutation_in_probe.append(a_deletion_range)
                            captures_indels = True

        # make a dictionary mapping position to insertion range
        # check if the position belongs to an insertion
        chromosome_position_insertion_dict = {}
        make_dict_chromosome_position_to_range(
            all_position_tumour_dict_insertion,
            chromosome_position_insertion_dict)

        # see which insertion does this position belong to
        the_insertion_range_list = chromosome_position_insertion_dict.get(
            position_in_genome)
        if the_insertion_range_list is not None:
            for the_insertion_range in the_insertion_range_list:
                # get the insertion's tumour set
                a_position_in_probe_tumour_set_insertion = all_position_tumour_dict_insertion.get(
                    the_insertion_range)
                # if the insertion does exist and is not recurrent
                if a_position_in_probe_tumour_set_insertion is not None and len(
                        a_position_in_probe_tumour_set_insertion) >= recurrent_definition:
                    position_list, insertion_range_end, insertion_range_start, semicolon_index \
                        = parse_chromosome_position_range(the_insertion_range)
                    # if the entire deletion can be covered by the probe
                    if insertion_range_end < str(probe_end) and insertion_range_start > str(probe_start):
                        # if the insertion range is not been covered
                        if the_insertion_range not in all_mutation_in_probe:
                            tumour_set_selected_position_and_probe.extend(a_position_in_probe_tumour_set_insertion)
                            all_position_in_probe.append(int(insertion_range_start))
                            all_position_in_probe.append(int(insertion_range_end))
                            all_mutation_in_probe.append(the_insertion_range)
                            captures_indels = True

    return captures_indels


"""other utility functions"""


def make_dict_chromosome_position_to_range(
        all_position_tumour_set_specific_mutation,
        chromosome_position_specific_mutation_dict):
    """make a dictionary mapping chromosome position to chromosome position range"""


    # 42-43, 43-44
    for specific_position_range, tumour_set in all_position_tumour_set_specific_mutation.items():
        position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
            specific_position_range)
        chromosome = specific_position_range[:semicolon_index]
        for position in position_list:
            position_in_genome = chromosome + ':' + str(position)
            chromosome_position_specific_mutation_dict.setdefault(position_in_genome, []).append(specific_position_range)


def parse_chromosome_position_range(chromosome_position_range: str):
    # example 2:3241-3276 (for mutation)
    # example 13:18598287..114344403 (for CNV)
    semicolon_index = chromosome_position_range.index(':')

    # this '8:30837618' is possible
    if '..' in chromosome_position_range:
        first_dot_index = chromosome_position_range.index('.')
        position_range_start = chromosome_position_range[
                               semicolon_index + 1:first_dot_index]
        position_range_end = chromosome_position_range[first_dot_index + 2:]
        position_list = []
    elif '-' not in chromosome_position_range:
        position_range_start = chromosome_position_range[
                               semicolon_index + 1:]
        position_range_end = position_range_start
        position_list = list(range(int(position_range_start),
                                   int(position_range_end) + 1))
    else:
        dash_index = chromosome_position_range.index('-')
        position_range_start = chromosome_position_range[
                               semicolon_index + 1:dash_index]
        position_range_end = chromosome_position_range[dash_index + 1:]
        position_list = list(range(int(position_range_start),
                                   int(position_range_end) + 1))


    return position_list, position_range_end, position_range_start, semicolon_index


def make_position_deletion_dict(position_deletion_dict,
                                position_range_tumour_deletion):
    """
    :param position_deletion_dict: A dictionary mapping a position to a position
    range for deletion
    :param position_range_tumour_deletion: A dictionary mapping a position range
    to a list of tumour for deletion mutation
    :return:
    """
    for deletion_position_range, tumour_set in position_range_tumour_deletion.items():
        position_list, position_range_end, position_range_start, semicolon_index = parse_chromosome_position_range(
            deletion_position_range)
        for position in position_list:
            position_deletion_dict[position] = deletion_position_range
