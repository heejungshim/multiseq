#' @useDynLib multiseq
#' @title A list containing counts, covariates and metadata from a ChIP-Seq experiment
#' with a pattern of differential Transcription Factor Binding.
#' @description This dataset consists of count data from a HudsonAlpha Institute (HaibTf)
#' ChIP-Seq experiment on Gm12878 and H1hesc cell lines with factor Yy1 that is part of
#' the ENCODE project. Mapped reads were downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/,
#' peaks were downloaded from http://ebi.edu.au/ftp/software/software/ensembl/encode/integration_data_jan2011/byFreeze/june2012/peaks/spp/optimal/.
#' Read counts in region chr1:11740409-11756792 were extracted from the bam files
#' with function \code{\link{get.counts}}.
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 4 by 2^14 matrix of counts (4 is the
#' number of samples and 2^14 is the number of adjacent bases where reads have been
#' counted),\code{read.depth} a 4-dimensional vector specifiying the total number
#' of counts per sample, \code{g} a 4-dimensional vector specifying a
#' binary covariate for each sample, \code{region} a string specifiying the region
#' reads were mapped to, \code{assembly} a string specifying the genome that reads
#' were mapped to, \code{samples} containing metadata about the samples.
#' @name dat
NULL

#' @title A list containing counts, covariates and metadata from a RNA-Seq experiment
#' with a pattern of differential expression.
#' 
#' @description This is a subset of the data collected in an RNA-Seq experiment on
#' lymphoblastoid cell lines derived from unrelated Nigerian individuals (Pickrell et al.,
#' Nature 464, 768-772, 2010). Specifically, it consists of a matrix of counts for 66
#' samples over a region of length 2^12 (chr12:113354417-113358512) where SNP rs10774671
#' was found to be associated with a differential expression pattern. Reads were
#' downloaded from http://eqtl.uchicago.edu/Home.html amd mapped to \code{hg19} with bwa;
#' read counts were extracted from the bam files with function \code{\link{get.counts}}.
#' Covariates (genotypes) are available from the International HapMap Project.
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 66 by 2^12 matrix of counts, \code{read.depth}
#' a 66-dimensional vector specifiying the total number of counts per sample, a string \code{SNP}
#' specifying the SNP name, \code{g} a 66-dimensional vector specifying the genotype (covariate)
#' for each sample, \code{assembly} specifying the genome that reads were mapped to,
#' \code{region} a string specifying the region of length 2^12 that reads were mapped to.
#'
#' @name OAS1
NULL


#' Multiseq: analysis of sequence data from multiple samples
#' @docType package
#' @name multiseqr
#NULL
