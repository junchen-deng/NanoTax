# NanoTax
NanoTax is intended to produce a table with both **contig information** and the corresponding **taxonomy** for output contigs from assembliers, such as [Canu](https://github.com/marbl/canu) and [Flye](https://github.com/fenderglass/Flye) for Nanopore long reads, and [Megahit](https://github.com/voutcn/megahit) for illumina short reads. 

**contig information** includes contig ID (**#ID**), average coverage (**Avg_fold**), contig length (**Length**) and GC content (**Read_GC**). **taxonomy** includes results from blastn (**blastn**), and/or results from blastx (**blastx**), and/or final taxonomy (**Taxonomy**) based on blastn and blastx. 

NanoTax takes three files as input: contigs in FASTA format, nanopore reads in FASTQ format, and customized nucleotide database in FASTA format. if customized protein datase (in FASTA format) was also provided, blastx (**blastx**) and final taxonomy (**Taxonomy**) will be added to the final table.

NanoTax runs in three steps: 
* Retrieving taxonomy of contigs by blastn (and blastx if protein database is provided)
* Mapping raw reads onto contigs and calculating basic contig information (GC %, coverage, length)
* Joining contig information and taxonomy to produce the final table

## Software requirement
* Python 3
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Minimap2](https://github.com/lh3/minimap2)
* BBMap [pileup.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh)

## Usage
```
usage: NanoTax_v1.2.py [-h] [-ONT_fastq <FASTQ>] [-r1 <FASTQ>] [-r2 <FASTQ>]
                       [-db_prot <path to protein database>] [-o <output dir>]
                       [-c <Number of Cores>]
                       <FASTA> <nuclotide database>

This script is intended to produce a table with both contig information (e.g.
average coverage, GC content) and the corresponding taxonomy for output
contigs from long-reads and short-reads assemblers (e.g. Canu, Flye for
Nanopore, Megahit for illumina)

positional arguments:
  <FASTA>               the path to the contig fasta file/folder
  <nuclotide database>  the path to the nuclotide database

optional arguments:
  -h, --help            show this help message and exit
  -ONT_fastq <FASTQ>, --ONT_fastq <FASTQ>
                        the path to the Nanopore reads fastq file/folder
  -r1 <FASTQ>, --r1 <FASTQ>
                        the path to the pair-end reads file/folder
  -r2 <FASTQ>, --r2 <FASTQ>
                        the path to the pair-end reads file/folder
  -db_prot <path to protein database>, --db_prot <path to protein database>
                        the path to the protein database
  -o <output dir>, --output_dir <output dir>
                        output directory name (default: Assigned_Taxonomy)
  -c <Number of Cores>, --cores <Number of Cores>
                        The number of CPU cores the script will use (default:
                        half of cores of the current server)
```

## Contact information
[Junchen Deng](https://github.com/junchen-deng) 
**Email**: junchen.deng@doctoral.uj.edu.pl 

## Upcoming features
* checking the required environment before running the script
