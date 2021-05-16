#' @title Split a region string into sequence name, region start position, and region end position.
#'
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position
#' @examples
#' region=split_region("chr1:11740409-11756792")
#' print(region)
#' print(region$chr)
#' print(region$start)
#' @export
#' @keywords internal
#' @return a list with elements \code{chr}, \code{start}, \code{end}. 
split_region <- function(region){
    split <- unlist(strsplit(region, "\\:|\\-"))
    if (length(split) != 3)
        stop("invalid region: example of a valid region is chr1:2345-234567")
    chr = split[1]
    start=as.numeric(split[2])
    end=as.numeric(split[3])

    locus.length=end-start+1
    
    if (start%%1 | end%%1 | locus.length<1) #check that start and end are integers
        stop("Incorrect parameters start and/or end")

    return(list(chr=chr, start=start, end=end))
}

#' @title Prepare input for multiseq function extracting counts in a genomic region from \code{bam}, \code{hdf5}, or \code{bigWig} files.
#'
#' @description This functions reads in a `samplesheet` and extracts read counts in a genomic region specified by parameter `region`, preparing input for \code{\link{multiseq}}. If in file `samplesheet` \code{samples$bamPath} is specified then this function extracts reads from the bam files in \code{samples$bamPath} using \code{samtools} (which must be in the USER's path) (no filter applied but option \code{onlyonend} can be used to extract reads from only one end if reads are paired end). Else if \code{samples$hdf5Path} is specified this function extracts reads from the \code{hdf5} files in \code{samples$hdf5Path} using the R package \code{rhdf5}. If \code{samples$bigWigPath} is specified this function extracts reads from \code{bigWig} files using the executable \code{bigWigToWig} (which must be in the USER's path).
#'
#' @param samplesheet: a string specifying the path to a samplesheet; the samplesheet should contain a column with header \code{SampleID} and a column with header either \code{bamPath} or \code{hdf5Path} or \code{bigWigPath}.
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position (see example below).
#' @param onlyonend: a bool, defaults to FALSE; use TRUE if input is in bam format and only first end of the paired end read should be used.
#' @export 
#' @return a matrix with \code{N} rows and \code{end-start+1} columns containing the number of reads that start at each base in the specified region in each sample. Rownames correspond to \code{samples$SampleID}
#' @examples
#'\dontrun{
#' setwd(file.path(path.package("multiseq"),"extdata"))
#' samplesheet="samplesheetEncode.txt"
#' region="chr1:11740409-11756792"
#' x=get.counts(samplesheet, region)
#' }
get.counts <- function(samplesheet=NULL, region=NULL, onlyoneend=FALSE){
    if (is.null(region) | is.null(samplesheet))
        stop("Invalid argument")
    #load samplesheet
    region       <- split_region(region)
    locus.length <- region$end-region$start+1
    #load counts
    samples <- read.table(samplesheet, stringsAsFactors=F, header=T)
    
    #load data
    M <- NULL
    if ("bamPath" %in% colnames(samples) & ! noExecutable("samtools")){
        #read bam files into matrix
        cmd <- paste0(region$chr,
                      ":",
                      region$start,
                      "-",
                      region$end)
        if (onlyoneend==TRUE)
            cmd <- paste0(cmd,"| awk '{if (!($7==\"=\" && $4>$8)) print}'")
        cmd <- paste0(cmd,                           
                      " | awk -v s='",
                      region$start,
                      "' 'BEGIN{start=0; count=0}",
                      "{st=$4; if ($4>=s){if (start==st) count+=1; ",
                      "else {if (start>0) print start, count; start=st; count=1} }}' ")
        
        for (bamfile in samples$bamPath){
            command <- paste0("samtools view ",
                              bamfile,
                              " ",
                              cmd)
            print(command)
            con <- pipe(command, open = "r")
            v   <- rep(0, locus.length)
            while(length(oneLine <- readLines(con, n = 1)) > 0){
                oneLine <- as.numeric(unlist(strsplit(oneLine," ")))
                v[oneLine[1]-region$start+1] <- oneLine[2]
            }
            close(con)
            M <- rbind(M, v)
        }
    }else if ("hdf5Path" %in%  colnames(samples)) {
                                        #read hdf5 files into R matrix
        for (h5file in samples$hdf5Path){
            print(paste0("Loading ", h5file))
            print(paste0("h5read(", h5file, ",", region$chr, ", index=list(", region$start, ":", region$end, ")"))
            M <- rbind(M, h5read(h5file, region$chr, index=list(region$start:region$end)))
        }
    }else if ("bigWigPath" %in%  colnames(samples) & !noExecutable("wigToBigWig")){  
        #read bigWig files into R matrix
        Mcol <- region$end-region$start+1
        for (bigWigfile in samples$bigWigPath){
            print(paste0("Loading ", bigWigfile))
            wigfile <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".wig")
            command=paste0("bigWigToWig ", bigWigfile,
                " -chrom=", region$chr,
                " -start=", region$start-1,
                " -end=", region$end,
                " stdout | grep -v fixed > " ,
                wigfile)
            print(command)
            system(command)
            M <- rbind(M, as.numeric(as.matrix(read.table(wigfile, stringsAsFactors=F, header=FALSE))))
            if(ncol(M)!=Mcol)
                stop(paste0("Error: input bigWig file does not (completely) cover region ", region$chr, ":", region$start, "-", region$end))
            file.remove(wigfile)
        }
    }else
    stop(paste("no input file provided and/or no executables available to read input:",
               "provide paths to input files (in hdf5, bigWig or bam format) in the samplesheet",
               "file and/or add required executables to user's PATH."))
    
    row.names(M) <- samples$SampleID
    
    return(M)
}

#' @title Extract transcript annotation (in a specified genomic region) from a file in \code{GenePred}.
#' format (a standard gene annotation format).
#' @param GenePredIn: a file in \code{GenePred} format containing gene annotation
#' @param region: a string specifying a genomic region
#' @export
#' @examples
#' GenePredIn <- file.path(path.package("multiseq"),"extdata","hg19.OAS1.refGene.part.gp")
#' region     <- "chr12:113354417-113358512"
#' get.transcripts(GenePredIn, region)
get.transcripts <- function(GenePredIn, region){
    genePred    <- data.frame(lapply(read.table(GenePredIn,
                                                fill=1,
                                                comment.char="",
                                                header=FALSE),
                                     as.character),
                              stringsAsFactors=FALSE)
    region      <- split_region(region)
    Transcripts <-  genePred[genePred[,2]==region$chr &((genePred[,4]<=region$start & genePred[,5]>=region$start)|genePred[,4]<=region$end & genePred[,5]>=region$start),]
    t           <- list(t=Transcripts)
    return(structure(t,class="transcripts"))
}

get.exons.start.end <- function(transcript){	
    exst <- as.numeric(strsplit(as.character(transcript[9]), ",")[[1]]) + 1
    exen <- as.numeric(strsplit(as.character(transcript[10]), ",")[[1]])
    return(list(exst=exst, exen=exen))
}


#' @title Plot transcripts; if \code{region} is specified then plot only the specified region.
#' @description The only required field for this function is Transcripts.
#' @param Transcripts:  output of \code{\link{get.transcripts}}
#' @param region: region to be plotted; defaults to NULL.
#' @param is.xaxis: bool, if TRUE plot \code{x} axis otherwise don't plot \code{x} axis; defaults to TRUE.
#' @export
#' @examples
#' region     <- "chr12:113354417-113358512" 
#' GenePredIn  <- file.path(path.package("multiseq"),"extdata","hg19.OAS1.refGene.part.gp") 
#' Transcripts <- get.transcripts(GenePredIn, region)
#' plot(Transcripts, region)
plot.transcripts <- function(Transcripts, region=NULL, is.xaxis=TRUE, axes=F, xlim=NULL, type=NULL, cex=NULL, font.main=1, ...){
    if (!is.null(region)){
        region <- split_region(region)
        chr    <- region$chr
        if (is.null(xlim))
            xlim <- c(region$start, region$end)
    }else{
        chr="Unknown Chromosome"
    }

    par(mgp=c(0,1,0))
    #draw x axis
    if (is.null(Transcripts$t)){
        y=1
        type="n"
    }else{    
        nr   <- nrow(Transcripts$t)
        tchr <- Transcripts$t[1,2]
        if (!is.null(region)){
            if (tchr!=region$chr)
                stop("Chromosomes are not consistent.")
        }else{
            chr=tchr
        }
        if (is.null(xlim))
            xlim=c(min(as.numeric(Transcripts$t[,4])), max(as.numeric(Transcripts$t[,5])))
        y="NA"
    }
    if (is.xaxis)
        xlab=paste("Position (Mb) on", chr)
    else
        xlab=""
    plot(y,
                 type=type,
                 xlim=xlim,
                 axes=axes,
                 ylim=c(0,1),
                 yaxt="n",
                 xlab=xlab,
                 ylab="",
                 font.main=1,
                 cex.main=cex,
                 cex.lab=cex,
                 ...)
   
    if (is.xaxis){
        tck=axTicks(1)
        tcklab=format(tck/1000000, nsmall=3)
        axis(1,
             at=tck,
             labels=tcklab,
             cex.lab=cex,
             cex.axis=cex)
    }
    #draw transcripts
    if (!is.null(Transcripts$t)){
        for (i in 1:nr){#set colors
            if (Transcripts$t[i,3] == "-") fill.col="red" else fill.col="blue"
            border.col=fill.col;
            
            if (nr==1){y=0.5} else y=(nr+1-i)*(1/(nr+1))
            
            rect(Transcripts$t[i,4],
                 y-0.001,
                 Transcripts$t[i,5],
                 y+0.001,
                 col="black",
                 border="black")
            trans <- get.exons.start.end(Transcripts$t[i,])
            
            ll=length(trans$exst)
            for (j in 1:ll)
                rect(trans$exst[j],
                     y-0.02,
                     trans$exen[j],
                     y+0.02,
                     col=fill.col,
                     border=border.col )
        }    
    }
}

#' @title Plot the output of \code{\link{multiseq}} (either the effect or the baseline).
#'
#' @description If `what=="effect"` this function will plot the effect and intervals where \code{\link{multiseq}} found strong effect, i.e., when zero is outside of +/- \code{z.threshold} * posterior standard deviation). If  `what=="baseline"` or `what=="log_baseline"` this function will plot either the baseline \code{exp(res$baseline.mean)} or the log of the baseline \code{res$baseline.mean}, respectively, and  intervals where \code{\link{multiseq}} found strong peaks at a specified threshold, i.e., when log(p.threshold) is below +/- \code{z.threshold} * posterior standard deviation of the log baseline. If x$region is defined then the \code{x} axis will use genomic coordinates.
#' 
#' @param x: multiseq output; if x$region is defined then the \code{x} axis will use genomic coordinates.
#' @param is.xaxis: bool, if TRUE plot \code{x} axis otherwise don't plot \code{x} axis.
#' @param z.threshold: a multiplier of the standard deviation.
#' @param p.threshold: this argument is only used when \code{what=="baseline"} or \code{what=="log_baseline"} to specify a threshold for the detection of peaks: if \code{res$baseline.mean-z.threshold*sqrt(baseline.var)>log(p.threshold)} then a peak is called; defaults to 1e-09.
#' @param what: a string, it can be either "baseline" or "log_baseline" or "effect".
#' @param highlight: a bool, if TRUE highlight intervals with strong peaks or strong signal; defaults to TRUE.shell
#' @export
#' @examples
#'\dontrun{
#' data(dat, package="multiseq")
#' res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
#' res$region <- dat$region
#' plot(res)
#' plot(res, what="baseline")
#' }
plot.multiseq <- function(x, is.xaxis=TRUE, z.threshold=2, p.threshold=1e-09, what="effect", highlight=TRUE, axes=F, type="l", col="green", main=NULL, ylim=NULL, xlab=NULL, ylab=NULL, cex=NULL, ...){
    if ((is.null(x$baseline.mean) | is.null(x$baseline.var)) & what=="baseline")
        stop("Error: no baseline in multiseq output")
    if ((is.null(x$effect.mean) | is.null(x$effect.var)) & what=="effect")
        stop("Error: no effect in multiseq output x")
    if (!(what=="effect" | what=="log_baseline" | what=="baseline"))
        stop("Error: wrong parameter 'what'")

    if (what=="effect"){
        ybottom     <- x$effect.mean - z.threshold*sqrt(x$effect.var)
        ytop        <- x$effect.mean + z.threshold*sqrt(x$effect.var)
        N           <- length(x$effect.mean)
        y           <- x$effect.mean
        wh.bottom   <- which(ybottom > 0)
        wh.top      <- which(ytop < 0)
        high.wh     <- sort(unique(union(wh.bottom, wh.top)))
        k           <- 0
    }else{
        ybottom     <- x$baseline.mean - z.threshold*sqrt(x$baseline.var)
        ytop        <- x$baseline.mean + z.threshold*sqrt(x$baseline.var)
        N           <- length(x$baseline.mean)
        k           <- log(p.threshold)
        high.wh     <- which(ybottom > k)
        if (what=="baseline"){
                y           <- exp(x$baseline.mean)
                main        <- paste(what, main)
                ymax        <- max(y)
                ymin        <- min(y)
                k           <- p.threshold
        }else{
            y           <- x$baseline.mean
        }
    }
    if (what!="baseline"){
        main        <- paste0(what, " +/- ", z.threshold, " s.d. ", main)
        ymax        <- max(ytop)
        ymin        <- min(ybottom)
    }
    
    if (is.null(ylim))
        ylim=c(ymin,ymax)
    if (is.xaxis){
        if (is.null(xlab)){
            if (!is.null(x$region)){
                region <- split_region(x$region)
                xlab    <- paste("Position (Mb) on", region$chr)
            }else{
                xlab    <- "Position"
            }
        }
    }else{
        xlab <- ""
    }
    if (is.null(ylab))
        ylab=""
    par(mgp=c(0,1,0))
    plot(y,
         type=type,
         ylim=ylim,
         main=main,
         col=paste("dark",col),
         xlab=xlab,
         ylab=ylab,
         axes=axes)
    #draw axes
    axis(2)
    if (is.xaxis){
        if (!is.null(x$region)){
            axis(1,
                 at=seq(1,N,ceiling(N/5)),
                 labels=format(seq(region$start,region$end,ceiling(N/5))/1000000, nsmall=3),
                 cex.lab=cex,
                 cex.axis=cex)
        }else
            axis(1)
    }
    if (what=="effect" | what=="log_baseline"){
        points(ytop, type=type, col=col)
        points(ybottom, type=type, col=col)
    }
    
    #draw intervals with effect or with peaks
    if (highlight){
        abline(h=k, col="red")  
        N.polygons  <- length(high.wh)
        if (N.polygons > 0)
            for(j in 1:N.polygons)
                rect(high.wh[j]-0.5, ymin-abs(ymin/5), high.wh[j]+0.5, ymax+abs(ymax/5),
                     col=rgb(1, 0, 0,0.5), border=NA, lty=NULL)
    }
}


#' @title helper function to get.intervals
#'
#' @description Output interval is in \code{bed} format (\code{start} is 0-based, \code{end} is 1-based).
#' @param mean: an estimated vector of means.
#' @param var: an estimated vector of variances.
#' @param what: a string, it can be either "baseline" or "log_baseline" or "effect".
#' @param z.threshold: a multiplier of the standard deviation.
#' @param p.threshold: this argument is only used when \code{what=="baseline"} or \code{what=="log_baseline"} to specify a threshold for the detection of peaks: if \code{mean-z.threshold*sqrt(var)>p.threshold} then a peak is called.
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position; defaults to NULL; if provided, the function will output the interval in genomic coordinates.
#' @export
#' @keywords internal 
get.intervals.utils <- function(mean, var, what, z.threshold, p.threshold, region){
    if (is.null(mean) | is.null(var))
        stop(paste("ERROR: no", what, "mean or", what, "var in multiseq output"))
    start    <- NULL
    end      <- NULL
    toreturn <- NULL
    sign     <- NULL
    if (what=="log_baseline" | what=="baseline"){
        toreturn$p.threshold <- p.threshold
    }else if (what=="effect"){
        x=(mean + z.threshold * sqrt(var) < p.threshold)
        xs=sum(x)
        if (xs>0){
            rlex=rle(x)
            boundaries=c(0,cumsum(rlex$lengths))
            for (i in 1:(length(boundaries)-1)){
                if (rlex$values[i]){ #is TRUE
                    start=c(start,boundaries[i]) #0-based
                    end=c(end,boundaries[i+1]) #1-based
                    sign=c(sign,"-")
                }
            }
        }
    }else{
        stop("Error: wrong parameter 'what'")
    }

    y=(mean - z.threshold * sqrt(var) > p.threshold)
    ys=sum(y)
    if (ys>0){
        rley=rle(y)
        boundaries=c(0,cumsum(rley$lengths))
        for (i in 1:(length(boundaries)-1)){
            if (rley$values[i]){ #is TRUE
                start=c(start,boundaries[i]) #0-based
                end=c(end,boundaries[i+1]) #1-based
                sign=c(sign,"+")
            }
        }
    }
    
    type="local"
    if (!is.null(start)){
        if (!is.null(region)){
            region       <- split_region(region)
            if (!is.null(region$start)&!is.null(region$end)&!is.null(region$chr)){
                start       <- region$start+start-1
                end         <- region$start+end
                type               <- "sequence"
                toreturn$chr       <- region$chr
            }else
                warning("WARNING: missing region start or region end; effect start and end are local and not relative to the sequence")
        }
    }
    toreturn$start       <- start
    toreturn$end         <- end
    toreturn$sign        <- sign
    toreturn$z.threshold <- z.threshold
    toreturn$type        <- type
    
    return(toreturn)
}
   


#' @title Print intervals where \code{\link{multiseq}} found strong effect or strong peaks. 
#'
#' @description If \code{what=="effect"} this function will print intervals where \code{\link{multiseq}} found strong effect, i.e., when zero is outside of +/- \code{z.threshold} * posterior standard deviation). If  \code{what=="baseline"} or \code{what=="log_baseline"}  this function will print intervals where \code{\link{multiseq}} found strong peaks at a specified threshold, i.e.,  when p.threshold is below +/- \code{z.threshold} * posterior standard deviation of the log baseline.
#'
#' @return This function returns a list with elements \cr
#' \item{chr}{a string specifying the sequence}
#' \item{start}{a vector specifying the start of each interval}
#' \item{end}{a vector specifying the end of each interval}
#' \item{sign}{can be "+" or "-" and indicates the sign of the effect or is always positive when \code{what="baseline"}} 
#' \item{z.threshold}{multiplier of the standard deviation}
#' \item{p.threshold}{threshold for peak detection} 
#' \item{type}{where \code{type} is either "local" - if \code{res$region} or \code{region} are not specified - or "sequence" otherwise.}
#' Output interval is in \code{bed} format (\code{start} is 0-based, \code{end} is 1-based).
#' @param res: \code{\link{multiseq}} output.
#' @param z.threshold: a multiplier of the standard deviation; default is 2
#' @param p.threshold: this argument is only used when \code{what=="baseline"} or \code{what=="log_baseline"} to specify a threshold for the detection of peaks: if \code{res$baseline.mean-z.threshold*sqrt(baseline.var)>log(p.threshold)} then a peak is called; default is 1e-09.
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position; defaults to NULL; if provided, the function will output the interval in genomic coordinates.
#' @param what: a string, it can be either "baseline" or "log_baseline" or "effect"; default is "effect"
#' @export
#' @examples
#'\dontrun{
#' data(dat, package="multiseq")
#' res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
#' get.intervals(res, what="effect", region=dat$region)
#' }
get.intervals <- function(res, z.threshold=2, p.threshold=1e-09, region=NULL, what="effect"){
    if (is.null(region))
        if (!(is.null(res$region)))
            region=res$region
    
    if (what=="baseline" | what=="log_baseline"){
        if (is.null(res$baseline.mean) | is.null(res$baseline.var)){
            stop("no baseline in multiseq output res") 
        }else{
            get.intervals.utils(res$baseline.mean, res$baseline.var, what, z.threshold, log(p.threshold), region)
        }
    }else if (what=="effect"){
        if (is.null(res$effect.mean) | is.null(res$effect.var)){
            stop("no effect in multiseq output res")
        }else{
            get.intervals.utils(res$effect.mean, res$effect.var, what, z.threshold, p.threshold=0, region)
        }
    }else{
        stop("ERROR: invalid parameter 'what'")
    }
}


#' @title Write effect intervals to a \code{bed} file.
#' @param intervals: output of \code{\link{get.intervals}}
#' @param bedfile: output \code{bed} file
#' @examples
#' \dontrun{
#' data("dat",package="multiseq")
#' res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
#' intervals <- get.intervals(res, region=dat$region, what="effect")
#' write.bed(intervals, "out.bed")
#' }
#' @export
write.bed <- function(intervals, bedfile){
    if (is.null(intervals)){
        stop("Invalid input intervals")
    }else if (intervals$type!="sequence"){
        warning("No input intervals")
        return()
    }else if (is.null(intervals$start) | is.null(intervals$end) | is.null(intervals$chr)){
        warning("No input intervals")
        return()
    }
    dir.create(dirname(bedfile), showWarnings=FALSE, recursive=TRUE)
    cat(paste0(intervals$chr,
               "\t",
               intervals$start,
               "\t",
               intervals$end,
               "\t",
               ".",
               "\t",
               "1000",
               "\t",
               intervals$sign,
               "\n",
               collapse=""),
        sep="",
        file=bedfile)
}
       
    
#' @title Write \code{\link{multiseq}} results to a compressed file.
#' @description The output file has two columns: first column the effect (or the baseline) mean and second column is the effect (or the baseline) posterior standard deviation ^2.
#' @param res: \code{\link{multiseq}} output.
#' @param file: path to the output file.
#' @param what: if \code{what} is "log_baseline" then write \code{\link{multiseq}} \code{res$baseline.mean} and \code{res$baseline.var}; if \code{what} is "effect" then write \code{\link{multiseq}} \code{res$effect.mean} and \code{res$effect.var} output to a compressed file; defaults to "effect".
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' data(dat, package="multiseq")
#' res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
#' write.gz(res, "./results.mean.sd2.gz")
#' }
write.gz <- function(res, file="results.mean.sd2.gz", what="effect"){
    dir.create(dirname(file), showWarnings = FALSE, recursive=TRUE)
    gz1      <- gzfile(file, "w")
    if (what=="log_baseline"){
        if (is.null(res$baseline.mean) | is.null(res$baseline.var) ){
            close(gz1)
            stop("no baseline in multiseq output res")
        }else{
            write.table(cbind(res$baseline.mean, res$baseline.var), col.names=FALSE, row.names=FALSE, file=gz1)
        }
    }else if (what=="effect"){
        if (is.null(res$effect.mean) | is.null(res$effect.var)){
           close(gz1)
           stop("no effect in multiseq output res")
       }else{
           write.table(cbind(res$effect.mean, res$effect.var), col.names=FALSE, row.names=FALSE, file=gz1)
       }
    }else{
        stop("Error: wrong parameter 'what'")
    }
    close(gz1)
}
