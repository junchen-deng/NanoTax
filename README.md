# NanoTax
NanoTax is intended to produce a table with both **contig information** and the corresponding **taxonomy** for output contigs from assembliers, such as [Canu](https://github.com/marbl/canu) and [Flye](https://github.com/fenderglass/Flye) for Nanopore long reads, and [Megahit](https://github.com/voutcn/megahit) for illumina short reads. 

**contig information** includes contig ID (**#ID**), average coverage (**Avg_fold**), contig length (**Length**) and GC content (**Read_GC**). **taxonomy** includes results from blastn (**blastn**), and/or results from blastx (**blastx**), and/or final taxonomy (**Taxonomy**) based on blastn and blastx. 

NanoTax takes at least three files as input: contigs in FASTA format, sequencing reads in FASTQ format, and customized nucleotide database in FASTA format. Options "**-ONT_fastq**" and "**-r1 -r2**" allow you to choose nanopore reads and illumina paired-end reads as input. If customized protein datase (in FASTA format) was provided, blastx (**blastx**) and final taxonomy (**Taxonomy**) will be added to the final table.

NanoTax runs in four steps: 
* Retrieving taxonomy of contigs by blastn (and blastx if protein database is provided)
* Producing a **TrueContigs.fasta** file that contains contigs with certain hits from blastn and/or blastx 
* Mapping raw reads onto true contigs and calculating basic contig information (GC %, coverage, length)
* Joining contig information and taxonomy to produce the final table

## Software requirement
* Python 3
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Minimap2](https://github.com/lh3/minimap2)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* BBMap [pileup.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh)
* [Diamond] (https://github.com/bbuchfink/diamond)

## Usage
```usage: NanoTax_v2.1.py [-h] [-db_2 <.FASTA>] [-bt <blast/diamond>] [-blastx_task <str>] [-blastn_task <str>] [-diamond_sen <str>] [-BAM] [-SAM] [-cov <from_pileup.sh>] [-bnf <outfmt6>] [-bxf <outfmt6>]
                       [--truecontigs] [-ONT_fastq <FASTQ>] [-r1 <FASTQ>] [-r2 <FASTQ>] [-o] [-prefix] [-c]
                       <contigs> <database.FASTA> <database type>

This script is intended to produce a table with both contig information (e.g. average coverage, GC content) and the corresponding taxonomy for output contigs from long-reads and short-reads assemblers
(e.g. Canu, Flye for Nanopore, Megahit for illumina)

positional arguments:
  <contigs>             the path to the contig fasta file/folder
  <database.FASTA>      the path to the database
  <database type>       database type (option: 'protein' 'dna')

optional arguments:
  -h, --help            show this help message and exit
  -db_2 <.FASTA>, --database_2 <.FASTA>
                        the path to the 2nd database. The db_type of the 2nd database has to be different from the 1st database
  -bt <blast/diamond>, --blast_tool <blast/diamond>
                        the type of tools for blast (default: blast) (option: 'blast', 'diamond')
  -blastx_task <str>, --blastx_task <str>
                        the task type of blastx (option: 'blastx' 'blastx-fast') (default: blastx)
  -blastn_task <str>, --blastn_task <str>
                        the task type of blastn (option: 'blastn' 'megablast') (default: megablast)
  -diamond_sen <str>, --diamond_sensitivity <str>
                        the sensitivity of diamond (option: 'fast' 'mid-sensitive' 'sensitive' 'more-sensitive' 'very-sensitive' 'ultra-sensitive') (default: sensitive)
  -BAM , --BAM          path to BAM file; enable skipping the mapping step
  -SAM , --SAM          path to SAM file; enable skipping the mapping step
  -cov <from_pileup.sh>, --coverage <from_pileup.sh>
                        path to coverage file from pileup.sh; enable to skip both mapping and pileup.sh
  -bnf <outfmt6>, --blastn_file <outfmt6>
                        path to blastn output; enable skipping the blastn step
  -bxf <outfmt6>, --blastx_file <outfmt6>
                        path to blastx output; enable skipping the blastx step
  --truecontigs         enable analysing only contigs with blast hits
  -ONT_fastq <FASTQ>, --ONT_fastq <FASTQ>
                        the path to the Nanopore reads fastq file/folder
  -r1 <FASTQ>, --r1 <FASTQ>
                        the path to the r1 pair-end reads file/folder
  -r2 <FASTQ>, --r2 <FASTQ>
                        the path to the r2 pair-end reads file/folder
  -o , --output_dir     output directory name (default: Assigned_Taxonomy)
  -prefix , --prefix    prefix of each output file name (default: output)
  -c , --cores          The number of CPU cores the script will use (default: 8)


```

## Contact information
[Diego Franco](https://github.com/diecasfranco) (diego.franco@uj.edu.pl)
[Junchen Deng](https://github.com/junchen-deng) (junchen.deng@doctoral.uj.edu.pl) 

## Upcoming features
* checking the required environment before running the script
