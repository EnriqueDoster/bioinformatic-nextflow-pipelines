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

from __future__ import division
from Bio import SeqIO
import sys

def assembly_stats(contigsMultifasta, sample_id, mini, est_genome_size, average_gene_size):
    
    contigsLength = []
    trimmedLength = []
    total = 0
    sum = 0
    thres = mini - 1
    GC_count = 0
    seqs = open(contigsMultifasta, 'r')

    # Create lists for Total Length and Trimmed Length
    for seq_record in SeqIO.parse(open(contigsMultifasta), 'fasta'):
        contigsLength.append(len(seq_record.seq))
        total += len(seq_record.seq)
        # Min Contig Length Threshold
        if len(seq_record.seq) > thres: 
            sum += len(seq_record.seq)
            trimmedLength.append(len(seq_record.seq))
        
    # Calculating GC content
    for seq_record in seqs.read():
        if seq_record.startswith('>'):
            continue
        else:
            if 'G' in seq_record:
                GC_count += 1
            if 'C' in seq_record:
                GC_count += 1
                
    GC_cont = float((GC_count/total)*100)

    # Sorting the Trimmed Contigs from Large to Small
    trimmedLength.sort()
    trimmedLength.reverse()
    
    # Theoretic NXX and NGXX Sizes
    teoN50 = sum / 2.0 
    teoNG50 = est_genome_size / 2.0
    teoN90 = sum * 0.9
    teoNG90 = est_genome_size * 0.9

    # Calculating Mean Contig Size
    meancontig = int(sum/len(trimmedLength))

    # Calculating Median Contig Size
    median = []
    medcon = []

    for con in trimmedLength:
        medcon.append(con)
        if len(medcon) > len(trimmedLength)/2:
            median.append(con)
            break

    # Checking N50 [bp] [# of Contigs]
    testSum = 0
    N50 = 0
    N50con = 0
    for con in trimmedLength:
        testSum += con
        N50con += 1
        if teoN50 < testSum:
            N50 = con
            break

    # Checking NG50 [bp] [# of Contigs]
    testSum = 0
    NG50 = 0
    NG50con = 0
    for con in trimmedLength:
        if sum < (est_genome_size/2):
            break
        testSum += con
        NG50con += 1
        if teoNG50 < testSum:
            NG50 = con
            break
        
    # Checking N90 [bp] [# of Contigs]
    testSum = 0
    N90 = 0
    N90con = 0
    for con in trimmedLength:
        testSum += con
        N90con += 1
        if teoN90 < testSum:
            N90 = con
            break

    # Checking NG90 [bp] [# of Contigs]
    testSum = 0
    NG90 = 0
    NG90con = 0
    for con in trimmedLength:
        if sum < est_genome_size/2:
            break
        testSum += con
        NG90con += 1
        if teoNG90 < testSum:
            NG90 = con
            break
            
    Xkb = 0
    for con in trimmedLength:
        if con > average_gene_size:
            Xkb += 1
    
    
    out = open(sample_id + '_assembly_stats.txt', 'w')
    out.write(str(sample_id) + ',num_contigs,' + str(len(contigsLength)) + '\n')
    out.write(str(sample_id) + 'num_of_contigs_longer_than' + str(average_gene_size) + 'bp,' + str(Xkb) + '\n')
    out.write(str(sample_id) + ',total_bp_length,' + str(total) + '\n')
    #out.write(str(sample_id) + '# trimmed contigs: ' + str(len(trimmedLength)) + '\n')
    #out.write(str(sample_id) + 'trimmed length [bp]: ' + str(sum) + '\n')
    out.write(str(sample_id) + ',GC_content,' + str(GC_cont) + '\n')
    out.write(str(sample_id) + ',min_contig_length,' + str(min(trimmedLength)) + '\n')
    out.write(str(sample_id) + ',median_contig_length,' + str(median[0]) + '\n')
    out.write(str(sample_id) + ',mean_contig_length,' + str(meancontig) + '\n')
    out.write(str(sample_id) + ',max_contig_length,' + str(max(trimmedLength)) + '\n')
    out.write(str(sample_id) + ',N50_length,' + str(N50) + '\n')
    out.write(str(sample_id) + ',num_N50_contig,' + str(N50con) + '\n')
    out.write(str(sample_id) + ',NG50_length,' + str(NG50) + '\n')
    out.write(str(sample_id) + ',num_NG50_contig,' + str(NG50con) + '\n')
    out.write(str(sample_id) + ',N90_length,' + str(N90) + '\n')
    out.write(str(sample_id) + ',num_N90_contig,' + str(N90con) + '\n')
    out.write(str(sample_id) + ',NG90_length,' + str(NG90) + '\n')
    out.write(str(sample_id) + ',num_NG90_contig,' + str(NG90con) + '\n')
    out.close()

if __name__ == "__main__":
    contigsMultifasta = sys.argv[1]
    sample_id = sys.argv[2]
    mini = int(sys.argv[3])
    est_genome_size = int(sys.argv[4])
    average_gene_size = int(sys.argv[5])
    
    assembly_stats(contigsMultifasta, sample_id, mini, est_genome_size, average_gene_size)

