#!/usr/bin/env python
#Magnify - A sam parser for insertion site detection

import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--reads_source', type=str, required=True)
parser.add_argument('--genome_source', type=str, required=True)
parser.add_argument('--marker_sources', type=str, required=True, nargs='+')
parser.add_argument('--granularity',type=int, required=False)
parser.add_argument('--threads',type=int,required=False)
parser.add_argument('--supplement_reads',type=int, required=False)
parser.add_argument('--min_matches', type=int, required=False)
args = parser.parse_args()

read_source = args.reads_source
read_source_path = os.getcwd()+"/"+args.reads_source
genome_source = args.genome_source
genome_source_path = os.getcwd()+"/"+args.genome_source
marker_sources = [os.getcwd()+"/"+i for i in args.marker_sources]

granularity_threshold = args.granularity
min_matches = args.min_matches

if args.threads is not None:
    threads = args.threads
else:
    threads = 1

if args.supplement_reads is not None:
    supplement_reads_tag = args.supplement_reads
    if supplement_reads_tag == 0:
        supp_tags = ['']
    if supplement_reads_tag == 1:
        supp_tags = ['SA']
    if supplement_reads_tag == 2:
        supp_tags = ['XA']
    if supplement_reads_tag == 3:
        supp_tags = ['SA','XA']
else:
    supp_tags = ['SA','XA']

#['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
#transgenes = ['pX330-U6-Cas9','pHDR2_backbone','pCAGGSS-Dre','CK557_LHX8','CK580_FIGLA','CK583_SOHLH1','CK596_ZNF281','T2A-tdTomato','T2A-mGreenLantern']

#Classes/Functions to deal with generating the final SAM file

for file in marker_sources:
    os.system(f'echo "$(cat {file})" >> magnify_{genome_source}_markers_combined.txt')

os.system(f'cp {genome_source_path} magnify_{genome_source}_collated.fa')

os.system(f'echo "$(cat magnify_{genome_source}_markers_combined.txt)" >> magnify_{genome_source}_collated.fa')


os.system(f'bwa index magnify_{genome_source}_collated.fa')

os.system(f'samtools fastq -@ {threads} {read_source_path} -1 magnify_{read_source}_R1.fastq -2 magnify_{read_source}_R2.fastq')

os.system(f'repair.sh in1=magnify_{read_source}_R1.fastq in2=magnify_{read_source}_R2.fastq out1=magnify_{read_source}_R1_fixed.fq out2=magnify_{read_source}_R2_fixed.fq outs=magnify_{read_source}_04_singletons.fq repair')

os.system(f'bwa mem -M -t {threads} magnify_{genome_source}_collated.fa magnify_{read_source}_R1_fixed.fastq magnify_{read_source}_R2_fixed.fastq > magnify_{read_source}_markers.sam')

with open(f'magnify_{genome_source}_markers_combined.txt') as f:
    marker_list = [i[1:] for i in f if i.startswith('>')]

markers = "|".join(marker_list)
#markers = 'pX330-U6-Cas9|pHDR2_backbone|pCAGGSS-Dre|CK557_LHX8|CK580_FIGLA|CK583_SOHLH1|CK596_ZNF281|T2A-tdTomato|T2A-mGreenLantern'
os.system(f"samtools view magnify_{read_source}_markers.sam |gawk '{{if (($3 ~ /^({markers})/)&&($7 ~ /chr([0-9]+|[XYM])/)) print $0}};' > magnify_{read_source}_markers_diff_chr.sam")
#--------------------------------------------------

#Classes/Functions that deal with final SAM file

class SamAlignment:
    def __init__(self, alignment): #defines all required parameters as specified in the SAM Manual to avoid dealing with slicing in future, as well as ALIGNMENT, a parameter containing the entire string, and OPTIONAL, a list containing all the optional fields

        self.ALIGNMENT = alignment

        fields = alignment.split('\t')

        self.QNAME = fields[0]
        self.FLAG = fields[1]
        self.RNAME = fields[2]
        self.POS = fields[3]
        self.MAPQ = fields[4]
        self.CIGAR = fields[5]
        self.RNEXT = fields[6]
        self.PNEXT = fields[7]
        self.TLEN = fields[8]
        self.SEQ = fields[9]
        self.QUAL = fields[10]
        try:
            self.OPTIONAL = fields[11:]
        except IndexError:
            self.OPTIONAL = None #if optional fields do not exist, self.OPTIONAL returns None

    def __str__(self): #making an alignment into a string returns the original line from the sam file
        return self.ALIGNMENT

    def optionalTag(self,tag_code): #finds optional tag (if it exists, otherwise returns None) and separates it into [TAG,TYPE,VALUE] as determined in the SAM manual
        if self.OPTIONAL == None:
            return None
        for tag in self.OPTIONAL:
            if tag.startswith(tag_code):
                return tag.split(':')

    def supplementaryAlignments(self,tag_list=['SA','XA']): #finds supplementary alignments corresponding to the tags input (it expects a list). tags that follow this convention are SA (chimeric reads) and XA (split reads)
        supplementary_alignment_list = [] 
        for i in tag_list:
            optional_tag = self.optionalTag(i)
            if optional_tag == None:
                continue
            supplementary_alignment_string = optional_tag[2] #finds the SA tag (indicating 'other canonical alignments') and sets its value 
            supplementary_alignments = supplementary_alignment_string.split(';')
            supplementary_alignment_list = supplementary_alignment_list + [a.split(',') for a in supplementary_alignments if len(a) >= 2] #returns a list of lists of type [rname, pos, strand, CIGAR, mapQ, NM]
        return supplementary_alignment_list            

    def position(self):
        return [self.RNAME,int(self.POS)]

    def matePosition(self):
        return [self.RNEXT,int(self.PNEXT)]

    def supplementaryPosition(self,tag_list=['SA','XA']):
        positions = []
        supplementaries = self.supplementaryAlignments(tag_list)
        for match in supplementaries:
            position = match[1]
            try:
                position = int(position)
            except:
                position = int(position[1:])

            positions.append([match[0],position])
        return positions

    def allPositions(self,tag_list=['SA','XA']): #returns position of the alignment, of its mate, and any supplementary alignments as a three-item list [[position],[mateposition],[supplementalaligmentposition]]
        required_positions = [self.position(),self.matePosition()]
        required_positions.extend(self.supplementaryPosition(tag_list))
        return required_positions

def group(samfile): #returns a dictionary containing all reference transgene chromosomes. each of those then contains another dictionary containing all the chromosomes that got mapped to, and each of those contains a dictionary containing the position of the mappings and how many times they were mapped to it.
    alignments = {}
    with open(samfile) as sam:
        for alignment in sam:
            alignment = SamAlignment(alignment)
            matches = alignment.allPositions(tag_list=supp_tags)
            ref_chromosome = matches[0][0] #takes name of ref chromosome
            mates_and_supplementaries = matches[1:] #takes positions of all matches and excludes the original read (we don't care where the read maps to on the original chromosome)
            try:
                chromosome_dict = alignments[ref_chromosome]
                for position in mates_and_supplementaries:
                    aligned_chrm_name = position[0]
                    try:
                        try:
                            chromosome_dict[aligned_chrm_name][abs(position[1])] += 1 
                        except:
                            chromosome_dict[aligned_chrm_name][abs(position[1])] = 1 
                    except:
                        chromosome_dict[aligned_chrm_name] = {abs(position[1]):1}  
            except KeyError:
                alignments[ref_chromosome] = {}
                for position in mates_and_supplementaries: #if there is no dict for the chromosome the sequence aligns to, it creates one. this for loop is necessary as to not discard the information in that alignment
                    alignments[ref_chromosome][position[0]] = {abs(position[1]):1} 
    return alignments

def compress(alignment_dict, granularity=500):
    readout_dict = {}
    for reference_chromosome, alignments in alignment_dict.items():
        readout_dict[reference_chromosome] = {}
        for aligned_chromosome, alignment_data in alignments.items():
            readout_dict[reference_chromosome][aligned_chromosome] = {}
            for location, repetitions in alignment_data.items():
                saved_locations = readout_dict[reference_chromosome][aligned_chromosome].keys()
                distances = [abs(location - saved_location) for saved_location in saved_locations]
                above_granularity = [distance > granularity for distance in distances]
                if all(above_granularity):
                    try:
                        readout_dict[reference_chromosome][aligned_chromosome][location] += repetitions
                    except KeyError:
                        readout_dict[reference_chromosome][aligned_chromosome][location] = 1
                else:
                    for possible_location, repetition in readout_dict[reference_chromosome][aligned_chromosome].items():
                        if abs(possible_location - location) < granularity:
                            readout_dict[reference_chromosome][aligned_chromosome][possible_location] += repetitions
    return readout_dict

def readout(insertion_dict, chr_filter, min_matches=1):
    print("Insertion Sites Found:")
    for read_chromosome, alignments in insertion_dict.items():
        print(read_chromosome+":")
        for align_chr, sites in alignments.items():
            if align_chr in chr_filter:
                continue
            else:
                print('\t'+align_chr+":")
                for site, repetitions in sites.items():
                    if repetitions >= min_matches:
                        print('\t'+'\t'+str(site)+": "+str(repetitions)+" matched")
                    else:
                        continue
    return True

#---------------------------------------------------



data = group(f'magnify_{read_source}_markers_diff_chr.sam')
insertions = compress(data,granularity_threshold)
readout(insertions,marker_list,min_matches)



