FROM rocker/rstudio

LABEL NAME="Bam filter" Version="1.0"
LABEL author="Samuel Soukup"
LABEL contact="soukup.sam(at)gmail.com"

RUN apt-get update
RUN apt-get install -y samtools \
    bwa \
    && apt-get clean
RUN apt-get clean

RUN mkdir /bam_filter
COPY bam_reduce.sh /bam_filter
COPY filter.awk /bam_filter
COPY bamfilter.R /bam_filter
COPY reference /reference

CMD ["/init"]

RUN echo 'source("/bam_filter/bamfilter.R")' >> /usr/local/lib/R/etc/Rprofile.site
