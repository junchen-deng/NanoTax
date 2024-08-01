# NanoTax
NanoTax is intended to produce a table with both **contig information** and the corresponding **taxonomy** for output contigs from assembliers, such as [Canu](https://github.com/marbl/canu) and [Flye](https://github.com/fenderglass/Flye) for Nanopore and PacBio long reads, and [SPAdes](https://github.com/ablab/spades) and [Megahit](https://github.com/voutcn/megahit) for illumina short reads. 

**contig information** includes contig ID (**#ID**), average coverage (**Avg_fold**), contig length (**Length**) and GC content (**Read_GC**). **taxonomy** includes results from blastn (**blastn**), and/or results from blastx (**blastx**), and final taxonomy (**Taxonomy**) based on blastn and blastx. 

NanoTax takes at least one file as input: contigs in FASTA format. Options "**-ONT_fastq**", "**-PB_fastq**" and "**-r1 -r2**" allow you to choose nanopore/PacBio long reads and illumina paired-end reads as input. 

In general, NanoTax runs in four steps: 
* Retrieving taxonomy of input contigs by blastn or blastx/diamond 
* Mapping raw reads, either long reads with **minimap2** or illumina with **bowtie2**, onto contigs 
* Calculating basic contig information (GC %, coverage, length) with **pileup.sh** using **bam** output from the mapper
* Joining contig information and taxonomy to produce the final taxonomy table

However, you could skip some of the steps by providing blastn/blastx/diamond output (options **-bnf -bxf**), BAM/SAM file (option **-BAM -SAM**), or coverage file (option **-cov**). 

If you are not interested in contigs without any blast hits, you could set **--truecontigs** and the script will keep only contigs with blast hits (output as **prefix.TrueContigs.fasta**). This way, the script will run faster as the mapping only works on filtered contigs. These removed contigs will not appear in the final taxonomy table.     

## Software requirement
* Python 3
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Minimap2](https://github.com/lh3/minimap2)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* BBMap [pileup.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh)
* [Diamond](https://github.com/bbuchfink/diamond) (optional; required if "**-bt diamond**" is set)

## database requrement
* You need to provide a custom database for blastn/blastx/diamond to function
* Only one type of database (DNA or Protein) is enough for the script to run
* Providing both databases will allow the script to summarize results from different databases in the final taxonomy table.   

## Help message
```
usage: NanoTax_v2.2.py [-h] [-db_nucl] [-db_prot] [-bt <blast/diamond>]
                       [-blastx_task <str>] [-blastn_task <str>]
                       [-diamond_sen <str>] [-BAM] [-SAM]
                       [-cov <from_pileup.sh>] [-bnf <outfmt6>]
                       [-bxf <outfmt6>] [--truecontigs] [-ONT_fastq <FASTQ>]
                       [-PB_fastq <FASTQ>] [-r1 <FASTQ>] [-r2 <FASTQ>] [-o]
                       [-prefix] [-c]
                       <contigs>

This script is intended to produce a table with both contig information (e.g.
average coverage, GC content) and the corresponding taxonomy for output
contigs from long-reads and short-reads assemblers (e.g. Canu, Flye for
Nanopore, Megahit for illumina)

positional arguments:
  <contigs>             the path to the contig fasta file/folder

optional arguments:
  -h, --help            show this help message and exit
  -db_nucl , --nucleotide_database 
                        the path to nucleotide database
  -db_prot , --protein_database 
                        the path to protein database
  -bt <blast/diamond>, --blast_tool <blast/diamond>
                        the type of tools for blast (default: blast) (option:
                        'blast', 'diamond')
  -blastx_task <str>, --blastx_task <str>
                        the task type of blastx (option: 'blastx' 'blastx-
                        fast') (default: blastx)
  -blastn_task <str>, --blastn_task <str>
                        the task type of blastn (option: 'blastn' 'megablast')
                        (default: megablast)
  -diamond_sen <str>, --diamond_sensitivity <str>
                        the sensitivity of diamond (option: 'fast' 'mid-
                        sensitive' 'sensitive' 'more-sensitive' 'very-
                        sensitive' 'ultra-sensitive') (default: sensitive)
  -BAM , --BAM          path to BAM file; enable skipping the mapping step
  -SAM , --SAM          path to SAM file; enable skipping the mapping step
  -cov <from_pileup.sh>, --coverage <from_pileup.sh>
                        path to coverage fiel from pileup.sh; enable to skip
                        both mapping and pileup.sh
  -bnf <outfmt6>, --blastn_file <outfmt6>
                        path to blastn output; enable skipping the blastn step
  -bxf <outfmt6>, --blastx_file <outfmt6>
                        path to blastx output; enable skipping the blastx step
  --truecontigs         enable analysing only contigs with blast hits
  -ONT_fastq <FASTQ>, --ONT_fastq <FASTQ>
                        the path to the Nanopore reads fastq file/folder
  -PB_fastq <FASTQ>, --PB_fastq <FASTQ>
                        the path to the PacBio reads fastq file/folder
  -r1 <FASTQ>, --r1 <FASTQ>
                        the path to the r1 pair-end reads file/folder
  -r2 <FASTQ>, --r2 <FASTQ>
                        the path to the r2 pair-end reads file/folder
  -o , --output_dir     output directory name (default: Assigned_Taxonomy)
  -prefix , --prefix    prefix of each output file name (default: output)
  -c , --cores          The number of CPU cores the script will use (default:
                        8)

```
## Examples: 
NanoTax_v2.2.py -db_prot proteins_references_wo_host_2.fasta -db_nucl F5_refs_bacteria.fasta -bt diamond -PB_fastq hifi_len500.fastq.gz --truecontigs -prefix STEMIN -o taxonomy -c 16 assembly.fasta

## Contact information
[Junchen Deng](https://github.com/junchen-deng) (junchen.deng@doctoral.uj.edu.pl) 
[Diego Franco](https://github.com/diecasfranco) (diego.franco@uj.edu.pl)

## Upcoming features
* checking the required environment before running the script
* fix issues that cannot create the nested repository
* fix issues that '--truecontigs' doesn't work if coverage file is provided  
