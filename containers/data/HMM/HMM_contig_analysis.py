#!/usr/bin/env python3

# cat AT_C2_S21_L002.non.host.contig.alignment.sam | ~/Dropbox/github/HMM_amrplusplus/test/HMM_contig_analysis.py mmarc_model_annotations.tsv /home/enrique/Dropbox/github/HMM_amrplusplus/test/test_counts mmarc_snpsearch_metadata.csv test_out

import re
import os.path
import sys

ltrans = {1: 'Class', 2: 'Mechanism', 3: 'Group'}
counts = {1: {}, 2: {}, 3: {}}

reads_mapped = 0
total_reads = 0


def load_mmarc_metadata(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')[1:]
        for line in data:
            if not line:
                continue
            entry = line.split('\t')
            ret[entry[1]] = entry[-3:]
    return ret


def load_mmarc_snp_metadata(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')[1:]
        for line in data:
            if not line:
                continue
            entry = line.split(',')
            ret[entry[0]] = entry
    return ret


def load_tblout(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        for line in data:
            if line.startswith('#') or not line:
                continue
            entry = line.split()
            contig = entry[0]
            if int(entry[6]) < int(entry[7]):
                ali = (int(entry[6]), int(entry[7]))
            else:
                ali = (int(entry[7]), int(entry[6]))
            ret.setdefault(contig, []).append((ali[0], ali[1], entry[2])) ## contig name, ali is start, stop position, annotation (checks if forward or reverse orientation)
    return ret


def parse_cigar(s):
    length = 0
    ret = re.findall(r'(\d+)([A-Z=]{1})', s)
    universe = {'X', 'P', 'I', 'N', 'D', '=', 'M'}
    for occ, op in ret:
        if op in universe:
            length += int(occ)
    return length


class SamParser:
    """This object takes as input a SAM file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input raw SAM file.
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.sam_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a SAM file")
        self.current_line = None
        self.reads_mapping = 0
        self.reads_total = 0
        self.header_lens = {}

    def __iter__(self):
        return self

    @property
    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                sam_line = sys.stdin.readline()  # read from stdin
            else:
                sam_line = self.sam_file.readline()  # read from file
            if not sam_line:
                return  # End of file
            if sam_line[0] != '@':  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if (int(temp[1]) & 4) == 0:
                    self.reads_mapping += 1
                    return temp[2], int(temp[3]), temp[5], temp[0], temp[8], temp[9], temp[10]  # RefName, 1-start, CIGAR |ED| ,readname, seq, dir, qual
        self.sam_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def __next__(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not fileIO
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            sys.stdout.write(
                "\n{:d} reads mapped out of {:d} total reads\n".format(self.reads_mapping, self.reads_total))
            sys.stdout.flush()
            global reads_mapped
            global total_reads
            reads_mapped = self.reads_mapping
            total_reads = self.reads_total
            raise StopIteration()
        else:
            return value


if __name__ == '__main__':
    R = load_mmarc_metadata(sys.argv[1]) ## Sam file is in stdin, model metadata
    D = load_tblout(sys.argv[2]) # Scan file
    SNP = load_mmarc_snp_metadata(sys.argv[3]) # SNP annotation file, the key is the mmarc model name
    SNP_contigs = 0
    non_SNP_contigs = 0
    list_SNP_contig = set()
    list_nonSNP_contig = set()
    for_reads = set()
    rev_reads = set()
    total_mapped_reads = set()
    SNP_mapped_reads = {}
    for refname, start, cigar, readname, dir, seq, qual in SamParser('-'):
        # RefName (contig), 1-start, CIGAR, readname, seq, qual
        stop = start + parse_cigar(cigar) - 1 # calculates stop position, based on start position of alignment in cigar string (to account for insertions, deletions)
        if refname not in D: ## Hammer table out file, checks if contig was not classified as a hmm model
            continue
        model_hits = set()
        #print(D[refname])
        for triplets in D[refname]: ## pulls reference name, looks it up in hammer file (has all annotations) ED. Contigs are dictionary keys
            if max(start, triplets[0]) <= min(stop, triplets[1]): #for each annotation, is it overlapping?
                if R[triplets[2]] != 'NA':
                    model_hits.add(triplets[2]) ## If sensical, add to annotations
        ## CHECK IF SNP
        if any(model in SNP for model in model_hits):
            SNP_contigs += 1
            list_SNP_contig.add(refname)
            if int(dir) > 0:
                for_reads.add(readname)
                SNP_readname = (str(readname) + "/1")
            else:
                rev_reads.add(readname)
                SNP_readname = (str(readname) + "/2")
            SNP_mapped_reads.setdefault(SNP_readname, (seq, qual))
            #print(model_hits)
        else:
            non_SNP_contigs += 1
            list_nonSNP_contig.add(refname)
            if model_hits: # list of models
                for x in model_hits: # x is model name
                    for e, a in enumerate(R[x]): ## Goes over each AMR taxa level, numbers from 0-2 starting with Class
                        annots = a.split('|') ## Split and divide count into each group ( can be more than two groups)
                        correct = len(model_hits) * len(annots) # to not double count reads
                        for annot in annots:
                            try:
                                counts[e + 1][annot] += float(1) / correct
                            except KeyError:
                                counts[e + 1][annot] = float(1) / correct
        total_mapped_reads.add(readname)

    samplename = sys.argv[4].split('/')[-1].split('.hmm')[0]
    with open(sys.argv[4] + ".SNP.fastq", 'w') as out:
        for SNP_read, seq in SNP_mapped_reads.items():
            out.write('@{}\n{}\n+\n{}\n'.format(SNP_read,seq[0],seq[1]))

    print(non_SNP_contigs,"Non-SNP contigs")
    print(SNP_contigs, "SNP contigs")
    print(SNP_contigs/non_SNP_contigs * 100 , "Percent SNP contigs")
    print(len(for_reads),"Forward reads")
    print(len(rev_reads),"Reverse reads")
    print(len(total_mapped_reads),"Total mapped reads")
    print(len(SNP_mapped_reads),"SNP mapped reads")


    with open(sys.argv[4] + ".SNP.stats.tsv", 'w') as out:
        out.write('{}\t{}\t{}\n'.format(samplename,"mapped_reads",total_mapped_reads))
        out.write('{}\t{}\t{}\n'.format(samplename,"SNP_reads",SNP_mapped_reads))
        out.write('{}\t{}\t{}\n'.format(samplename,"SNP_contigs",SNP_contigs))
        out.write('{}\t{}\t{}\n'.format(samplename,"nonSNP_contigs",non_SNP_contigs))

    with open(sys.argv[4] + ".counts", 'a') as out:
        #samplename = sys.argv[4].split('/')[-1].split('.hmm')[0]
        for level, ldict in counts.items():
            for k, count in ldict.items():
                out.write(
                    '{},{},{},{}\n'.format(samplename,ltrans[level], k,
                                                       str(count)))
    with open(sys.argv[4] + ".group_counts.tsv", 'a') as cout:
        #samplename = sys.argv[4].split('/')[-1].split('.hmm')[0]
        cout.write('Sample\tGene\tHits\tGene Fraction\n')
        for level, ldict in counts.items():
            if level == 3:
                for k, count in ldict.items():
                    if not k:
                        continue
                    cout.write(
                        '{}\t{}\t{}\t80\n'.format(samplename, k,
                                                       str(count)))
