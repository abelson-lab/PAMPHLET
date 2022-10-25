import csv
import sys
from typing import List


def search_cosmic_download_file(cosmic_file_name: str, gene_list_file: str):

    gene_list = []
    with open(gene_list_file) as gene_list_file:
        line = gene_list_file.readline()
        gene_list.append(line)

    gene_set = set(gene_list)
    wanted_genes_list = []

    with open(cosmic_file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if row[0] in gene_set:
                wanted_genes_list.append([row[0], row[2], row[3]])

    with open('wanted_genes.csv', mode='w') as wanted_genes_file:
        wanted_genes_file_writer = csv.writer(wanted_genes_file, delimiter=',',
                                              quotechar='"',
                                              quoting=csv.QUOTE_MINIMAL)

        wanted_genes_file_writer.writerows(wanted_genes_list)


if __name__ == "__main__":
    search_cosmic_download_file(sys.argv[1], sys.argv[2])
