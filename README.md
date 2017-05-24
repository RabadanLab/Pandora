Pandora
=======

Identification and Discovery of Tumor Associated Microbes via RNAseq

**Introduction**

Pandora is a multi-step pipeline to find pathogen sequences in RNAseq data. 
It includes modules for host separation, assembly, blasting contigs, and orf discovery.
As input, Pandora takes paired fastq files; as output, it produces a report.

**Dependencies**

The following programs must be in your `PATH`:

- python 2.7.x
- [Samtools](http://www.htslib.org/) 1.4 (*note*: use older versions of samtools at your own risk - samtools is often not backwards compatible)
- [STAR](https://github.com/alexdobin/STAR)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK279671/) 2.3.x
- [featureCounts](http://subread.sourceforge.net/) (Subread)

Pandora depends on the following Python modules:

- [Biopython](http://biopython.org/wiki/Main_Page)
- [Pandas](http://pandas.pydata.org/)
- scipy

The exact list, with versions, is provided in the `requirements.txt` file.
And the best way to install these is the usual best practice of starting a [virtualenv](https://virtualenv.pypa.io/en/stable/) and running:

```
pip install -r requirements.txt
```

**Workflow**

To accomplish diverse tasks, Pandora has various subcommands (like, say, the program git).
The primary subcommand is `scan`, which is a pipeline comprising the following steps:

1. Subtraction of reads mapping to host genome
2. De-Novo assembly of remaining reads
3. BLAST of assembled contigs
4. ORF search in contigs of unknown origin
5. Filter and parse blast results into tidy human-readable report

The `aggregate` subcommand [...].

**Additional Files**

Pandora requires various references and annotation files.

For `scan` step 1, please provide:
- a host genome indexed for STAR
- a host genome indexed for bowtie2
- (optional) a gtf describing the genes of the host

For `scan` step 3, please provide:
- the BLAST nucleotide collection nt database at ftp://ftp.ncbi.nlm.nih.gov/blast/db/

For `scan` step 4, you can optionally provide:
- the BLAST protein collection nr database at ftp://ftp.ncbi.nlm.nih.gov/blast/db/

For `scan` step 5, you can optionally provide:
- a text file of "blacklist" non-pathogen taxids for filtering. If you do not provide one, the script will use `resources/blacklist.txt` by default. This list contains any taxid children of the nodes chordata (Taxonomy ID: 7711) or "other sequences" (Taxonomy ID: 28384)
- the `names.dmp` file mapping taxID to names from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

Because there are a considerable number of files involved, you can specify their paths with a configuration file instead of command line flags.
See `pandora.config.txt` for example formatting.
Note that options specified as flags take precedence over options specified via the configuration file.

**Usage Examples**

```
pandora.py scan -id patient1 -r1 mate_1.fastq.gz -r2 mate_2.fastq.gz --gzip --refstar /path/ref/STAR --refbowtie /path/ref/bowtie/hg19 -db /path/ref/blastdb/nt --taxid2names /path/names.dmp
```

Here is an example command using a configuration file:

```
pandora.py scan -id patient1 -r1 mate_1.fastq.gz -r2 mate_2.fastq.gz --gzip --verbose -c pandora.config.txt
```

**Notes**

Currently, Pandora makes use of the [Oracle Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) by default.
The reason for this is that blast is computationally intensive, embarrassingly parallelizable, and lends itself very nicely to cluster computing.
You can turn this off with the `--noSGE` flag, but blast will be very slow.

Note that RNA-seq enriched for poly-A transcripts will miss prokaryotic pathogens.

Status: Active Development
