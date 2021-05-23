#! /usr/bin/env python3
# Authors: Diego Franco, Junchen Deng
# Affiliation: Symbio Group, Institute of Enviromental Sciences, Jagiellonian University, Poland
# Version 1.2


import os, argparse, itertools, sys, multiprocessing, csv, time

parser = argparse.ArgumentParser(
			description='This script is intended to produce a table with both contig information (e.g. average coverage, GC content) and the corresponding taxonomy for output contigs from long-reads and short-reads assemblers (e.g. Canu, Flye for Nanopore, Megahit for illumina)'	
)

# positional arguments
parser.add_argument("fasta", metavar='<FASTA>', help="the path to the contig fasta file/folder")
parser.add_argument("db_nucl", metavar='<nuclotide database>', type=str, help="the path to the nuclotide database")

# optional arguments
parser.add_argument("-ONT_fastq", "--ONT_fastq", metavar='<FASTQ>', help="the path to the Nanopore reads fastq file/folder")
parser.add_argument("-r1", "--r1", metavar='<FASTQ>', help="the path to the pair-end reads file/folder")
parser.add_argument("-r2", "--r2", metavar='<FASTQ>', help="the path to the pair-end reads file/folder")
parser.add_argument("-db_prot", "--db_prot", metavar='<path to protein database>', type=str, help="the path to the protein database")
parser.add_argument("-o", "--output_dir", metavar='<output dir>', type=str, help="output directory name (default: Assigned_Taxonomy)", default="Assigned_Taxonomy")
parser.add_argument('-c', '--cores', metavar='<Number of Cores>', type=int,
			help='The number of CPU cores the script will use (default: half of cores of the current server)', default=multiprocessing.cpu_count()/2)

# if no arguments were given, printing the help message (args = "--help")
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

fasta = args.fasta
db_nucl = args.db_nucl
ONT_fastq = args.ONT_fastq
r1_fastq = args.r1
r2_fastq = args.r2
db_prot = args.db_prot
prefix = fasta.split("/")[-1].split(".fa")[0]
cores = args.cores
output_dir = args.output_dir

start_time = time.time()
# creating output directory according to the input (default="Assigned_Taxonomy")
os.system("mkdir %s" % output_dir)

# retrieving taxonomy with blastn
print("\n\n Performing blastn... ")
os.system("blastn -task blastn -query {} -db {} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {} > {}/{}_blastn.txt".format(fasta, db_nucl, cores, output_dir, prefix))
print("blastn done\n")

# retrieving taxonomy with blastx (protein database)
if db_prot != None:
    print("\n\n Performing blastx...")
    os.system("blastx -query {} -db {} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {} > {}/{}_blastx.txt".format(fasta, db_prot, cores, output_dir, prefix))
    print("blastx done\n")

# import taxonomy table as a dictionary
with open("{}/{}_blastn.txt".format(output_dir, prefix), "r") as Table:
    Table_tax = {}
    for line in Table: 
        LINE = line.strip().split()
        header = LINE[0]
        taxonomy = LINE[1]
        Table_tax[header] = taxonomy

if db_prot != None:
    with open("{}/{}_blastx.txt".format(output_dir, prefix), "r") as Table:
        Table_tax_prot = {}
        for line in Table: 
            LINE = line.strip().split()
            header = LINE[0]
            taxonomy = LINE[1]
            Table_tax_prot[header] = taxonomy
        
## output a fasta file that contains contigs getting hits from blast
# make sure that each sequence only occupy one line
with open(fasta, "r") as FASTA:
    contig_dict ={}
    contig = ""
    for line in FASTA:
        if line.startswith(">"):
            if contig != "":    # save the existing header and sequence before overwritting 
                contig_dict[header] = contig
                contig = ""
            header = line.strip("\n").split()[0][1:]    # remove the ">" from the header for the following steps
        else:
            contig += line.strip("\n").upper()
    contig_dict[header] = contig

# output
with open("{}/{}.TrueContigs.fasta".format(output_dir, prefix), "w") as OUTPUT:
    for header in contig_dict.keys():
        if db_prot != None:
            if header in Table_tax_prot.keys() or header in Table_tax.keys():
                print(">"+ header, "\n", contig_dict[header], "\n", end = "", sep = "", file = OUTPUT)
        elif header in Table_tax.keys():
            print(">"+ header, "\n", contig_dict[header], "\n", end = "", sep = "", file = OUTPUT)

# mapping reads onto true contigs
print("\n\n Mapping... ")
if ONT_fastq != None: 
    print("Mapping long reads with Minimap2...")
    os.system ("minimap2 -ax map-ont {2}/{3}.TrueContigs.fasta {0} -t {1} > {2}/{3}.sam".format(ONT_fastq, cores, output_dir, prefix))
    print("Minimap2 done\n")
elif r1_fastq != None and r2_fastq != None:
    print("Mapping short reads with bowtie2...")
    os.system("bowtie2-build {0}/{1}.TrueContigs.fasta {0}/{1}.TrueContigs.fasta".format(output_dir, prefix))
    os.system ("bowtie2 --threads {2} -x {3}/{4}.TrueContigs.fasta -1 {0} -2 {1} --no-unal -S {3}/{4}.sam".format(r1_fastq, r2_fastq, cores, output_dir, prefix))
    
os.system ("pileup.sh in={0}/{1}.sam out={0}/{1}_coverage.txt".format(output_dir, prefix))
print("Coverage done\n")
os.system ("rm {0}/*.sam {0}/*.bt2".format(output_dir))

## combining information in taxonomy and coverage tables
# import coverage table as a list of lists
# the first line will be categories, i.e. "#ID, Avg_fold, Length, Ref_GC, Covered_percent, Covered_bases, Plus_reads, Minus_reads, Read_GC, Median_fold, Std_Dev"
with open("{}/{}_coverage.txt".format(output_dir, prefix), "r") as Table:
    Table_cov = []
    for line in Table: 
        LINE = line.strip().split()
        Table_cov.append(LINE)

row_count = len(Table_cov)

# add the taxonomy from blastn table
Table_cov[0].append("blastn")
for row in range(1, row_count):
    if Table_cov[row][0] in Table_tax.keys(): 
        header = Table_cov[row][0]
        taxonomy = Table_tax[header]
        Table_cov[row].append(taxonomy)
    else: 
        Table_cov[row].append("NA")

# add the taxonomy from blastx table
if db_prot != None: 
    Table_cov[0].append("blastx")
    for row in range(1, row_count):
        if Table_cov[row][0] in Table_tax_prot.keys(): 
            header = Table_cov[row][0]
            taxonomy = Table_tax_prot[header]
            Table_cov[row].append(taxonomy)
        else: 
            Table_cov[row].append("NA")
    
    # provide a summerized taxonomy only when protein database was given        
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
                      
# output table
col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn"]
if db_prot != None: 
    col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn", "blastx", "Taxonomy"]
    
col_pos = []
for name in col_keep: 
    col_pos.append(Table_cov[0].index(name))

with open("{}/{}_taxonomy_cov.txt".format(output_dir, prefix), "w") as OUTPUT:
    for row in Table_cov: 
        col_position = 0
        for col in row[:-1]:
            if col_position in col_pos:
                print(col, file = OUTPUT, end = "\t")
            col_position += 1
        print(row[-1], file = OUTPUT)

end_time = time.time()
print("running time: {} seconds".format(end_time-start_time))     
print("Na zdrowie! Salud! 乾杯! gānbēi(干杯^-^)! Saude! Cheers! Skål!")



