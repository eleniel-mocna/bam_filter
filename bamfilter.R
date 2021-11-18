#' regex_filter -> regex of searched sequence (e.g.: "ATTGA[GC]AG")
#' unfiltered_bam -> path to bam that should be filtered.
#' NAME              -> Folder name for outputs of this filtering is 'NAME_<regex_filter>'  [~reads1_path]
#' reference_path -> path to reference file (e.g.: "/reference/ucsc.hg19.fasta") [/reference/ucsc.hg19.fasta]
#' remove_unfiltered -> Remove unfiltered bam from output folder? [FALSE]
filter_bam <- function(regex_filter,
                       unfiltered_bam,
                       NAME=NULL,
                       reference_path=NULL,
                       remove_unfiltered=FALSE)
{
    filter_fq( regex_filter=regex_filter,
               unfiltered_bam=unfiltered_bam,
               NAME=NAME,
               reference_path, reference_path
               remove_unfiltered=remove_unfiltered)
}

#' regex_filter      -> regex of searched sequence (e.g.: "ATTGA[GC]AG") [NULL]
#' reads1_path       -> path to the first reads file (e.g.:"2843_Bc_F1_S46_L001_R1_001.fastq") [NULL]
#' reads2_path       -> path to the second reads file (e.g.:"2843_Bc_F1_S46_L001_R2_001.fastq") [NULL]
#' reference_path    -> path to reference file (e.g.: "/reference/ucsc.hg19.fasta") [/reference/ucsc.hg19.fasta]
#' unfiltered_bam    -> Is ignored, if reads1_path & reads2_path are filled
#'                     else path to bam to be filtered. [NULL]
#' NAME              -> Folder name for outputs of this filtering is 'NAME_<regex_filter>' [~reads1_path]
#' remove_unfiltered -> Remove unfiltered bam from output folder? [FALSE]
filter_fq <- function( regex_filter=NULL,
                        reads1_path=NULL,
                        reads2_path=NULL,
                        NAME=NULL,
                        unfiltered_bam=NULL,
                        reference_path=NULL,
                        remove_unfiltered = FALSE)
{
    if (is.null(NAME))
    {
        NAME=
            paste(
                strsplit(
                        tail(
                            strsplit("../../M554-10-4_S28_L001_R1_001.fastq", "/")[[1]],
                            1),
                        ".", fixed=TRUE)[[1]][1],
                 "_",
                 regex_filter,
                 sep="")            
    }
    skip_filter <- is.null(regex_filter)
    if (skip_filter)
    {
        print("Filter not given, exiting.")
        return()
    }
    system(paste("mkdir", NAME))
    skip_alignment <- is.null(reads1_path) || is.null(reads2_path)

    old_unfiltered_bam = unfiltered_bam
    unfiltered_bam = paste(NAME, "/unfiltered.bam", sep="")
    if(!skip_alignment)
    {
        if(is.null(reference_path))
        {
            reference_path = "/reference/ucsc.hg19.fasta"
            print(paste("Using default reference:", reference_path))
        }
        system(paste("bwa mem -t 20", reference_path, reads1_path, reads2_path, "| samtools view -h | samtools sort | samtools view -h -b >", unfiltered_bam), intern = TRUE)
        system(paste("samtools index", unfiltered_bam))
    }
    else 
    {
        print("SKIPPING ALIGNMENT") # Get old bam into the resuls folder
        system(paste("samtools view -h ", old_unfiltered_bam, " | samtools sort | samtools view -h -b > ", unfiltered_bam, sep=""))        
        system(paste("samtools index", unfiltered_bam))
    }

    filtered_bam=paste(NAME, "/filtered.bam", sep="")
    system(paste("/bam_filter/bam_reduce.sh", regex_filter, unfiltered_bam, filtered_bam))
    system(paste("samtools index", filtered_bam))

    if (remove_unfiltered)
    {
        system(paste("rm",  unfiltered_bam))
        system(paste("rm ", unfiltered_bam, ".bai", sep=""))
    }
}