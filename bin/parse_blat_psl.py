'''

   Calculation of [Total # of Contigs], [Total Length], [Total # of trimmed Contigs], [Trimmed Length], [GC content],
   [Min Contig Size [bp]], [Median Contig Size [bp]], [Mean Contig Size [bp]], [Max Contig Size [bp]],
   [N50[bp] [# of Contigs]], [NG50[bp] [# of Contigs]], [N90 [bp] [# of Contigs]], [NG90 [bp] [# of Contigs]],
   [Total # of Contigs > Average Gene Size]

   This code creates an output.txt file with all of the statistics

   Usage: python assembly_stats.py <assembly.fasta> <sample_ID> <minimum contig size> <estimated genome size> <average gene size>

   Author: Nicolas Schmelling
   Edited by: Enrique Doster

'''

#python parse_blat_psl.py SRR4425656_blat_results.psl SRR4425656 megares_modified_annotations_v2.0.csv MEGmobile_annotations_clean.csv

from __future__ import division
from Bio import SeqIO
from Bio import SearchIO
import sys



def parse_psl(blat_psl,sample_id):
    blat_qresult = SearchIO.parse(blat_psl, "blat-psl")
    search_dict = SearchIO.to_dict(blat_qresult)
    num_contig_multiple_hits = 0
    mean_hit_counts = 0
    out_gene_alignments = open(sample_id + '_gene_alignments.csv', 'w')
    out_gene_alignment_stats = open(sample_id + '_gene_alignment_stats.csv', 'w')
    # Go through list of results with the contig name as the index
    for contig in search_dict:
        contig_hit_object = search_dict[contig]
        #print(contig)
        #print(len(contig_hit_object))
        if len(contig_hit_object) > 1:
            contig_alignments = set()
            num_contig_multiple_hits += 1
            # make set of all gene alignments to a each contig
            for result in contig_hit_object:
                contig_alignments.add(result.id)
            # check if any headers in megares and megmobile annotations are both present in the set of alignments
            if any([x in megares_annot.keys() for x in contig_alignments]) and any([x in megmobile_annot.keys() for x in contig_alignments]):
                #print(contig)
                #print(contig_alignments)    
                for gene in contig_alignments:
                    write_alignments = ('{},{},{}\n'.format(sample_id,contig, gene))
                    out_gene_alignments.write(write_alignments)   
            contig_alignments.clear()
    print(len(search_dict)) # total number of contigs
    print(num_contig_multiple_hits) # number of contigs with multiple hits
    write_stats_len =  ('{},num_total_contigs,{}\n'.format(sample_id, len(search_dict)))
    write_stats_mapped =  ('{},num_multiple_mapped_contigs,{}\n'.format(sample_id,num_contig_multiple_hits ))
    out_gene_alignment_stats.write(write_stats_len)
    out_gene_alignment_stats.write(write_stats_mapped)
    out_gene_alignments.close()
    out_gene_alignment_stats.close()
    
    
def load_megares_annotations(file):
    ret = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')
        for line in data[1:]:
            if not line:
                continue
            entry = line.split(',')
            ret.setdefault(entry[0], entry[1:])
    return ret
    # returns dictionary object with headers as the indices and the values containing all other annotation columns
    
def load_megmobile_annotations(file):
    ret = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')
        for line in data[1:]:
            if not line:
                continue
            entry = line.split(',')
            ret.setdefault(entry[0], entry[1:])
    return ret

    

if __name__ == "__main__":
    blat_psl = sys.argv[1]
    sample_id = sys.argv[2]
    megares_annot = load_megares_annotations(sys.argv[3])
    megmobile_annot = load_megmobile_annotations(sys.argv[4])
    parse_psl(blat_psl, sample_id)
