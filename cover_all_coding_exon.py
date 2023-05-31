from others import write_output_excel


def make_probe_these_gene(reference_genome_filename):
    # first collect all gene name
    gene_name_dict = {}
    with open(reference_genome_filename) as refseq_genes:
        for gene_transcript in refseq_genes:
            gene_transcript_info_list = gene_transcript.split('\t')
            gene_name = gene_transcript_info_list[12]
            gene_name_dict[gene_name] = ''

    cover_these_gene = input("which gene do you want to cover in its"
                             " entirety, separated by comma ")
    genes_to_cover = cover_these_gene.split(',')

    for gene in genes_to_cover:
        if gene not in gene_name_dict:
            print(gene, "is not found")
            genes_to_cover.remove(gene)

    if len(genes_to_cover) == 0:
        print("no genes to cover")
        return

    # note only coding transcript are included in this file
    # its protein coding (NM_, XM_) vs non-coding (NR_, XR_)
    # see https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
    # next collect all coding exons range for each gene transcript
    with open(reference_genome_filename) as refseq_genes:
        gene_coding_exon_ranges_dict = {}
        # each line correspond to one gene transcript
        # each line is a tab separated string
        for gene_transcript in refseq_genes:

            gene_transcript_info_list = gene_transcript.split('\t')
            chromosome = gene_transcript_info_list[2]
            gene_name = gene_transcript_info_list[12]
            coding_region_start = gene_transcript_info_list[6]
            coding_region_end = gene_transcript_info_list[7]
            exons_starts = gene_transcript_info_list[9]
            exons_ends = gene_transcript_info_list[10]

            # if the user want to cover this gene, and it has coding regions
            # note that if it does not have coding region,
            # coding_region_start == coding_region_end
            if gene_name in genes_to_cover and coding_region_start < coding_region_end:

                exons_starts_list = exons_starts.strip(',').split(',')
                exons_ends_list = exons_ends.strip(',').split(',')

                assert len(exons_starts) == len(exons_ends)
                # for each exon, we need to exclude the UTR
                gene_transcript_coding_exon_ranges = []
                for exon_range in range(len(exons_starts_list)):

                    this_exon_start = exons_starts_list[exon_range]
                    this_exon_end = exons_ends_list[exon_range]

                    # this exon lies entirely outside the coding region
                    if this_exon_end < coding_region_start:
                        continue
                    elif this_exon_start > coding_region_end:
                        continue

                    if this_exon_start < coding_region_start:
                        this_coding_exon_start = coding_region_start
                    elif this_exon_start >= coding_region_start:
                        this_coding_exon_start = this_exon_start
                    else:
                        print('something is wrong')

                    if this_exon_end > coding_region_end:
                        this_coding_exon_end = coding_region_end
                    elif this_exon_end <= coding_region_end:
                        this_coding_exon_end = this_exon_end
                    else:
                        print('something is wrong')

                    gene_transcript_coding_exon_ranges.append([int(this_coding_exon_start), int(this_coding_exon_end)])

                chromosome_num_gene_name = chromosome + ',' + gene_name
                gene_coding_exon_ranges_dict.setdefault(chromosome_num_gene_name, []).extend(
                    gene_transcript_coding_exon_ranges)

    range_list = [['gene name', 'chromosome', 'coding exon range start',
                   'coding exon range end']]
    # for all coding exons of a gene
    for chromosome_num_gene_name, exon_ranges in gene_coding_exon_ranges_dict.items():

        chromosome_num_gene_name_list = chromosome_num_gene_name.split(',')
        chromosome = chromosome_num_gene_name_list[0]
        gene_name = chromosome_num_gene_name_list[1]

        exon_ranges.sort()
        # add 3 bases beyond the actual coding exon to cover splice sites
        for exon_range in exon_ranges:
            exon_range[0] -= 3
            exon_range[1] += 3

        # Sort the unique exon ranges by their starts
        exon_ranges.sort()
        overlap_merged = [exon_ranges[0]]
        # insert first exon range into overlap_merged
        # for every other exon
        for exon_range in exon_ranges[1:]:

            this_exon_start = exon_range[0]
            overlap_merged_last_exon = overlap_merged[-1]
            overlap_merged_last_exon_start = overlap_merged_last_exon[0]
            overlap_merged_last_exon_end = overlap_merged_last_exon[-1]
            # if exon ranges overlap (in other words, if this exon, starts
            # in the middle of the last exon in overlap_merged)
            if overlap_merged_last_exon_start <= this_exon_start <= overlap_merged_last_exon_end:
                # then change the last exon's end to whichever one is bigger (
                # last exon end or this exon end
                overlap_merged[-1][-1] = max(overlap_merged[-1][-1], exon_range[-1])
            # if this exon start after the last exon
            else:
                overlap_merged.append(exon_range)

        for exon_range in overlap_merged:
            range_list.append([
                gene_name, chromosome, exon_range[0], exon_range[1]
            ])

    write_output_excel(range_list, 'range_entire_gene.xlsx')
