#! /usr/bin/env python3
#Authors: Diego Franco, Junchen Deng
#Affiliation: Symbio Group, Institute of Enviromental Sciences, Jagiellonian University, Poland
#Version 1.0


import os, argparse, itertools, sys, multiprocessing, csv


#parser.add_argument("--blast_evalue", help="setting for e-value cutoff for blast, must be in form 1e-X", type=str, default="1e-10")

parser = argparse.ArgumentParser(
			description='This script is intended to produce a table with both contig information (e.g. average coverage, GC content) and the corresponding taxonomy for output contigs from long-reads and short-reads assemblers (e.g. Canu, Flye for Nanopore, Megahit for illumina)'	
)
parser.add_argument("fasta", metavar='<FASTA>', help="the path to the contig fasta file/folder")
parser.add_argument("-ONT_fastq", "--ONT_fastq", metavar='<FASTQ>', help="the path to the Nanopore reads fastq file/folder")
parser.add_argument("-r1", "--r1", metavar='<FASTQ>', help="the path to the pair-end reads file/folder")
parser.add_argument("-r2", "--r2", metavar='<FASTQ>', help="the path to the pair-end reads file/folder")
parser.add_argument("db_nucl", metavar='<nuclotide database>', type=str, help="the path to the nuclotide database")
parser.add_argument("-db_prot", "--db_prot", metavar='<path to protein database>', type=str, help="the path to the protein database")
parser.add_argument("-o", "--output_dir", metavar='<output dir>', type=str, help="output directory name (default: Assigned_Taxonomy)", default="Assigned_Taxonomy")
parser.add_argument('-c', '--cores', metavar='<Number of Cores>', type=int,
			help='The number of CPU cores the script will use (default: half of cores of the current server)', default=multiprocessing.cpu_count()/2)

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

fasta = args.fasta
ONT_fastq = args.ONT_fastq
r1_fastq = args.r1
r2_fastq = args.r2
db_nucl = args.db_nucl
db_prot = args.db_prot
prefix = fasta.split("/")[-1].split(".fa")[0]
cores = args.cores
output_dir = args.output_dir

#os.system("conda activate base")

os.system("mkdir %s" % output_dir)

print("\n\n Performing blastn... ")
os.system("blastn -task blastn -query %s -db %s -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads %d > %s/%s_blastn.txt" % (fasta, db_nucl, cores, output_dir, prefix))
print("blastn done\n")

if db_prot != None:
    print("\n\n Performing blastx...")
    os.system("blastx -query %s -db %s -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads %d > %s/%s_blastx.txt" % (fasta, db_prot, cores, output_dir, prefix))
    print("blastx done\n")

print("\n\n Mapping... ")

if ONT_fastq != None: 
    print("Mapping long reads with Minimap2...")
    os.system ("minimap2 -ax map-ont %s %s -t %d > %s/%s.sam" % (fasta, ONT_fastq, cores, output_dir, prefix))
    print("Minimap2 done\n")
elif r1_fastq != None and r2_fastq != None:
    print("Mapping short reads with bowtie2...")
    os.system("bowtie2-build {} {}".format(fasta, fasta))
    os.system ("bowtie2 --threads {} -x {} -1 {} -2 {} --no-unal -S {}/{}.sam".format(cores, fasta, r1_fastq, r2_fastq, output_dir, prefix))
    
os.system ("pileup.sh in=%s/%s.sam out=%s/%s_coverage.txt" % (output_dir, prefix, output_dir, prefix))
print("Coverage done\n")
os.system ("rm %s/*sam" % output_dir)


OUTPUT = open("%s/%s_taxonomy_cov.txt" % (output_dir, prefix), "w")


## import taxonomy table as a dictionary
Table = open("%s/%s_blastn.txt" % (output_dir, prefix), "r")
Table_tax = {}
for line in Table: 
    LINE = line.strip().split()
    Table_tax[LINE[0]] = LINE[1]
Table.close()

## import coverage table as a list of lists
Table = open("%s/%s_coverage.txt" % (output_dir, prefix), "r")
Table_cov = []
for line in Table: 
    LINE = line.strip().split()
    Table_cov.append(LINE)
Table.close()

row_count = len(Table_cov)

## assign the taxonomy of blastn
Table_cov[0].append("blastn")
for row in range(1, row_count):
    if Table_cov[row][0] in Table_tax.keys(): 
        contig = Table_cov[row][0]
        Table_cov[row].append(Table_tax[contig])
    else: 
        Table_cov[row].append("NA")

## assign the taxonomy of blastx
if db_prot != None:
    Table = open("%s/%s_blastx.txt" % (output_dir, prefix), "r")
    Table_tax = {}
    for line in Table: 
        LINE = line.strip().split()
        Table_tax[LINE[0]] = LINE[1]
    Table.close()
    
    Table_cov[0].append("blastx")
    for row in range(1, row_count):
        if Table_cov[row][0] in Table_tax.keys(): 
            contig = Table_cov[row][0]
            Table_cov[row].append(Table_tax[contig])
        else: 
            Table_cov[row].append("NA")
            
    Table_cov[0].append("Taxonomy")
    for row in range(1, row_count):
        if Table_cov[row][-2] == "NA" and Table_cov[row][-1] == "NA": 
            Table_cov[row].append("NA")
        elif Table_cov[row][-2] == "NA":
            tax = Table_cov[row][-1].split("_")[0]
            Table_cov[row].append(tax)
        else:
            tax = Table_cov[row][-2].split("_")[0]
            Table_cov[row].append(tax) 
        
        if "mito" in Table_cov[row][-3] or "mito" in Table_cov[row][-2]:
            Table_cov[row][-1] = "mitochondrion"
        if Table_cov[row][-1] == "Acyrthosiphon": 
            Table_cov[row][-1] = "Host"
            
           
## output table
col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn"]
if db_prot != None: 
    col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn", "blastx", "Taxonomy"]
    
col_pos = []
for  name in col_keep: 
    col_pos.append(Table_cov[0].index(name))

for row in Table_cov: 
    col_position = 0
    for col in row[:-1]:
        if col_position in col_pos:
            print(col, file = OUTPUT, end = "\t")
        col_position += 1
    print(row[-1], file = OUTPUT)
    

print("Na zdrowie! Salud! 乾杯! gānbēi(干杯^-^)! Saude! Cheers! Skål!")




