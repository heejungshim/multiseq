#' @useDynLib multiseq
#' @title OAS1 gene data set 
#' @description This is a subset of the data collected in an RNA-Seq experiment on
#' lymphoblastoid cell lines derived from unrelated Nigerian individuals (Pickrell et al.,
#' Nature 464, 768-772, 2010). Specifically, it consists of a matrix of counts for 66
#' samples over a region of length 2^12 (gene OAS1; chr12:113354417-113358512) where SNP rs10774671
#' was found to be associated with splicing (see Figure 3 in Pickrell et al.). Reads were
#' downloaded from http://eqtl.uchicago.edu/Home.html amd mapped to \code{hg19} with bwa.
#' Genotypes are available from the International HapMap Project.
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 66 by 2^12 matrix of counts, \code{read.depth}
#' a 66-dimensional vector specifiying the total number of counts per sample, a string \code{SNP}
#' specifying the SNP name, \code{g} a 66-dimensional vector specifying the genotype for each sample, 
#' \code{assembly} specifying the genome that reads were mapped to,
#' \code{region} a string specifying the region of length 2^12 that reads were mapped to.
#' @name OAS1
NULL


#' @title DNase seq data over chr17:10160989-10162012  
#' @description This is a subset of the data collected in a DNase-seq experiment on
#' lymphoblastoid cell lines derived from unrelated Nigerian individuals (Degner et al., Nature
#' 482, 390–4, 2012). Specifically, it consists of a matrix of counts for 70
#' samples over a region of length 2^10 (chr17:10160989-10162012) where SNP chr17.10161485 
#' was found to be associated with chromatin accessibility (see Figure 2 in Shim and Stephens, 
#' Ann. Appl. Stat. 9 665–686, 2015). DNase-seq reads and genotype data were downloaded from 
#' http://eqtl.uchicago.edu/dsQTL_data/.
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 70 by 2^10 matrix of counts, \code{read.depth}
#' a 70-dimensional vector specifiying the total number of counts per sample, a string \code{SNP}
#' specifying the SNP name, \code{g} a 70-dimensional vector specifying the genotype at the SNP for each sample, 
#' \code{region} a string specifying the region of length 2^10 that reads were mapped to.
#' @name chr17.10160989.10162012.DNase.seq
NULL

#' @title ATAC-seq data over chr1:111764939-111765962
#' @description This is a subset of the data collected in an ATAC-seq experiment on 3 copper treated 
#' samples and 3 control samples (Shim et al. 2021). Specifically, it consists of a matrix of counts 
#' for 6 samples over a region of length 2^10 (chr1:111764939-111765962) where chromatin accessibility 
#' is differentially expressed between two conditions (see Figure 3 in Shim et al. 2021). 
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 6 by 2^10 matrix of counts, \code{read.depth}
#' a 6-dimensional vector specifiying the total number of counts per sample, \code{g} a 
#' 6-dimensional vector specifying the group of each sample (copper treated samples vs control samples),
#' \code{region} a string specifying the region of length 2^10 that reads were mapped to, \code{overall.logLR} 
#' a number indicating log likelihood ratio for overall expression obtained from DESeq2 output, 
#' \code{overall.effect} a number indicating the effect size for overall expression obtained from DESeq2 output, 
#' \code{overall.effect.var} a number indicating the variance of the effect size obtained from DESeq2 output.
#' @name chr1.111764939.111765962.ATACseq
NULL

#' @title Empirical distribution of the multiseq test statistic from ATAC-seq analysis in Shim et al, 2021. 
#' 
#' @description This is the output from our ATAC-seq analysis in Shim et al, 2021. We applied multiseq to two controls 
#' (media vs. ethanol) for the 242,714 regions and this data set contains the resulting 242,714 test statistics. 
#' The two controls are expected to have no differences in chromatin accessibility. Thus, we can construct the empirical 
#' null distribution of the multiseq test statistic using the 242,714 test statistics in this data set. 
#' 
#' @docType data
#' @keywords datasets
#' @format A vector of length 242714. 
#' @name ATACseq.multiseq.stat.null
NULL


#' multiseq: multi-scale Poisson process approaches for differential expression analysis or association analysis of high-throughput sequencing data
#' @docType package
#' @name multiseqr
#NULL
