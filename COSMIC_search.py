import csv
import sys
from typing import List


def search_cosmic_lymphoid_mutation_file(file_type: str, cosmic_file_name: str,
                                         gene_list_file: str):

    gene_list = []
    with open(gene_list_file) as gene_list_file:
        line = gene_list_file.readline()
        gene_list.append(line)
    gene_set = set(gene_list)

    important_column_list = []
    important_column_list = determine_important_column(file_type,
                                                       important_column_list)
    wanted_genes_list = []
    match_gene_extract_cell(cosmic_file_name, gene_set, important_column_list,
                            wanted_genes_list)

    with open('wanted_genes.csv', mode='w') as wanted_genes_file:
        wanted_genes_file_writer = csv.writer(wanted_genes_file, delimiter=',',
                                              quotechar='"',
                                              quoting=csv.QUOTE_MINIMAL)
        wanted_genes_file_writer.writerows(wanted_genes_list)


def match_gene_extract_cell(cosmic_file_name, gene_set, important_column_list,
                            wanted_genes_list):
    # turn out the only difference between a tsv and csv is the delimiter
    # in terms of reading it
    with open(cosmic_file_name) as csv_file:
        # csv_reader = csv.reader(csv_file, delimiter=',')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[0] in gene_set:
                # depending on the file, add these columns of this row
                row_sub_list = []
                for column in important_column_list:
                    row_sub_list.append(row[column])
                wanted_genes_list.append(row_sub_list)


def determine_important_column(file_type: str, important_column_list: List[int]):
    if file_type == 'lymphoid mutation':
        # gene name, CDS mutation, AA mutation
        important_column_list = [0, 2, 3]
    elif file_type == 'lymphoid CNV':
        # gene name
        important_column_list = [0]
    elif file_type == 'complete mutation':
        # gene name, histology, subtype 1, subtype 2, subtype 3, mutation CDS
        # mutation AA, GRCh build number, mutation genome position
        important_column_list = [0, 11, 12, 13, 14, 19, 20, 24, 25]
    elif file_type == 'complete CNV':
        # gene name, histology, subtype 1, subtype 2, subtype 3, GRCh build number
        # genomic position
        important_column_list = [3, 9, 10, 11, 12, 18, 19]
    return important_column_list


if __name__ == "__main__":
    # search the genes from sys.argv[3] in sys.argv[2]
    # search in lymphoid mutation
    # search in lymphoid CNV
    # search in complete mutation
    # search in complete CNV
    if not (sys.argv[1] == 'lymphoid mutation' or sys.argv[1] == 'lymphoid CNV'
            or sys.argv[1] == 'complete mutation' or sys.argv[1] == 'complete CNV'):
        sys.exit('not one of the option, lymphoid mutation, lymphoid CNV,'
                 'complete mutation, complete CNV')
    search_cosmic_lymphoid_mutation_file(sys.argv[1], sys.argv[2], sys.argv[3])


