# Bam filter

## Functionality

This is a utility for filtering a bam file by a given string (regex pattern). It's got two interfaces: 1) lower level bash script and 2) R functions.

## Installation

### Docker installation (recommended)

The easiest way to run this filter is using an rstudio/rocker docker interface.

- Firstly copy a ucsc.hg19 reference into `reference` directory in the root directory.
- Then build the docker using: `docker build -t bam_filter_eleniel .`
- After it is built create a new container using _create\_docker.sh_ script: `./create_docker.sh <port> <name> <mount_directory> <password>`

Now you can connect to this docker via rocker/rstudio interface on your chosen port (Most often <http://localhost>:\<port>). In case of any difficulties during the setup see [rocker/rstudio website](https://hub.docker.com/r/rocker/rstudio) for more information.

### Installation into an existing unix enviroment

If for any reason running this utility in a docker is not possible, it is possible to import this functionality into an existing system. This might require more advanced steps and unforeseen difficulties.

- Firstly install samtools if needed using `sudo apt-get install samtools` or `sudo yum install samtools`
- The shell interface should now work as intended (see section Usage/Shell interface)
- For usage in Rstudio you need to source the `bamfilter.R` file into you R enviroment. This is done by `source("<path to bamfilter.R>")` in Rstudio.
- Now the Rstudio interface should work as intented, only `reference_path` argument will not work on basi settings and you need to enter your path to the reference every time.

## Usage

### Rstudio interface

Rstudio interface offers two functions, `filter_fq(...)` and `filter_bam(...)`:

#### filter_fq(...)

This function takes as input two fastQ files and a filter. It runs_bwa mem_ algorithm on the fastQ files and filters the resulting BAM file by given filter. The outputs are then saved to folder `<NAME>_<regex_filter>`.
It takes the following arguments:

- `regex_filter`      -> regex of searched sequence (e.g.: "ATTGA[GC]AG") [NULL]
- `reads1_path`       -> path to the first reads file  (_.fastq_) [NULL]
- `reads2_path`       -> path to the second reads file (_.fastq_) [NULL]
- `reference_path`    -> path to reference file (e.g.: "/reference/ucsc.hg19.fasta") [/reference/ucsc.hg19.fasta]
- `unfiltered_bam`    -> Is ignored, if reads1\_path & reads2_path are filled
                      else path to bam to be filtered. [NULL]
- `NAME`              -> Folder name for outputs of this filtering is 'NAME\_<regex\_filter>'
- `remove_unfiltered` -> Remove unfiltered bam from output folder? [FALSE]

#### filter_bam(...)

This function takes as input a BAM file and a filter. It filters the BAM by the given filter. The outputs are then saved to folder `<NAME>_<regex_filter>`.
It takes the following arguments:

- `regex_filter` -> regex of searched sequence (e.g.: "ATTGA[GC]AG")
- `unfiltered_bam` -> path to bam that should be filtered.
- `NAME`              -> Folder name for outputs of this filtering is 'NAME_<regex_filter>'
- `reference_path` -> path to reference file (e.g.: "/reference/ucsc.hg19.fasta")
- `remove_unfiltered` -> Remove unfiltered bam from output folder? [FALSE]
