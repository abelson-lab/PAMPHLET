import sys


def which_genes_were_mutated_in_healthy(gene_list_file, healthy_mutation_file):

    gene_list = []
    with open(gene_list_file) as gene_list_file:
        lines = gene_list_file.readlines()
        stripped_lines = []
        for line in lines:
            stripped_line = line.strip()
            stripped_lines.append(stripped_line)
        gene_list.extend(stripped_lines)
    l_chip_gene_set_1_all = set(gene_list)

    print(l_chip_gene_set_1_all)

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

    which_genes = l_chip_gene_set_1_all.intersection(l_chip_gene_set_1_healthy_mutation)

    if which_genes == l_chip_gene_set_1_healthy_mutation:
        print("its a strict subset")

    print(which_genes)

    print(len(which_genes))

if __name__ == "__main__":
    # mutation file, CNV file, and gene list file
    which_genes_were_mutated_in_healthy(sys.argv[1], sys.argv[2])
