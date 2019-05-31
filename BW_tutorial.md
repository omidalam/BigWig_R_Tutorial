Working with Bigwig files in R
================
Omid Gholamalamdari

# Importing BigWig files to R

I use rtacklayer package which is part of Bioconductor to import BigWigs
to R.

``` r
# How to install libs
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("rtracklayer", version = "3.8")
# BiocManager::install("GenomicRanges", version = "3.8")

library(rtracklayer)
library(GenomicRanges)
```

## Import function

Let’s downoad a K562 H3K4me3 BW from Encode and import it to R. import()
function of the rtracklayer imports a genomic file(i.e BW, bed,
bedgraph, …) as a GenomicRanges object into R. You can read more about
GenomicRanges
[here](https://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).
Examining the file shows a range and score. Score is the height in BW
plot over the specified range.

``` r
k562.h3k4me3<-import("./sample_file/k562_h3k4me3_20kb.bw")
head(k562.h3k4me3)
```

    ## GRanges object with 6 ranges and 1 metadata column:
    ##       seqnames        ranges strand |              score
    ##          <Rle>     <IRanges>  <Rle> |          <numeric>
    ##   [1]     chr1       1-20000      * |                  0
    ##   [2]     chr1   20001-40000      * |                  0
    ##   [3]     chr1   40001-60000      * |                  0
    ##   [4]     chr1   60001-80000      * | 0.0101521201431751
    ##   [5]     chr1  80001-100000      * | 0.0253790691494942
    ##   [6]     chr1 100001-120000      * | 0.0456807799637318
    ##   -------
    ##   seqinfo: 128 sequences from an unspecified genome

\#\#Binning function This is a function to bin genomic data, it gets a
BW (as GenomicRanges object) and bins it. The default is 20kb bins.

``` r
BW_bnr<-function(BW,bin=2e4){
  #gets a BE i GenomicRanges and bin it over fixed distances
  require(rtracklayer)
  require(GenomicRanges)
  bins<-tileGenome(seqinfo(BW),tilewidth = bin,cut.last.tile.in.chrom = TRUE)
  BW_Rle<-coverage(BW,weight="score")
  BW_bin<-binnedAverage(bins,BW_Rle,"score")
  return(BW_bin)
}
```

``` r
k562.h3k4me3.2e4<-BW_bnr(k562.h3k4me3,bin=2e4) 
# You can also export it and see how it looks.
export(k562.h3k4me3.2e4,"sample_file/k562_h3k4me3_20kb.bw")
```

Note that the range sizes are in 20kb increments

``` r
head(k562.h3k4me3)
```

    ## GRanges object with 6 ranges and 1 metadata column:
    ##       seqnames        ranges strand |              score
    ##          <Rle>     <IRanges>  <Rle> |          <numeric>
    ##   [1]     chr1       1-20000      * |                  0
    ##   [2]     chr1   20001-40000      * |                  0
    ##   [3]     chr1   40001-60000      * |                  0
    ##   [4]     chr1   60001-80000      * | 0.0101521201431751
    ##   [5]     chr1  80001-100000      * | 0.0253790691494942
    ##   [6]     chr1 100001-120000      * | 0.0456807799637318
    ##   -------
    ##   seqinfo: 128 sequences from an unspecified genome

Let’s import another bigwig file, a son\_k562\_tsa

``` r
k562.tsa<-import("sample_file/k562_c1r1_20k_mw20k_hg38.bw")
head(k562.tsa)
```

    ## GRanges object with 6 ranges and 1 metadata column:
    ##       seqnames        ranges strand |            score
    ##          <Rle>     <IRanges>  <Rle> |        <numeric>
    ##   [1]     chr1   10000-29999      * | 1.41149997711182
    ##   [2]     chr1   30000-49999      * | 1.43850004673004
    ##   [3]     chr1   50000-69999      * | 1.10230004787445
    ##   [4]     chr1   70000-89999      * | 1.06350004673004
    ##   [5]     chr1  90000-109999      * | 1.09889996051788
    ##   [6]     chr1 110000-129999      * | 1.12170004844666
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome

Note that the ranges are 20kb as before but the starting point is
diiferent. You need to make them uniform. You also need to make sure the
total number of bins are same.

``` r
k562.tsa
```

    ## GRanges object with 151570 ranges and 1 metadata column:
    ##            seqnames              ranges strand |            score
    ##               <Rle>           <IRanges>  <Rle> |        <numeric>
    ##        [1]     chr1         10000-29999      * | 1.41149997711182
    ##        [2]     chr1         30000-49999      * | 1.43850004673004
    ##        [3]     chr1         50000-69999      * | 1.10230004787445
    ##        [4]     chr1         70000-89999      * | 1.06350004673004
    ##        [5]     chr1        90000-109999      * | 1.09889996051788
    ##        ...      ...                 ...    ... .              ...
    ##   [151566]     chrX 155970000-155989999      * | 1.03840005397797
    ##   [151567]     chrX 155990000-156009999      * | 1.22580003738403
    ##   [151568]     chrX 156010000-156029999      * | 1.55460000038147
    ##   [151569]     chrX 156030000-156049999      * |  1.7257000207901
    ##   [151570]     chrX 156050000-156069999      * |                0
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome

``` r
k562.h3k4me3.2e4
```

    ## GRanges object with 155031 ranges and 1 metadata column:
    ##            seqnames            ranges strand |              score
    ##               <Rle>         <IRanges>  <Rle> |          <numeric>
    ##        [1]     chr1           1-20000      * |                  0
    ##        [2]     chr1       20001-40000      * |                  0
    ##        [3]     chr1       40001-60000      * |                  0
    ##        [4]     chr1       60001-80000      * | 0.0101521201431751
    ##        [5]     chr1      80001-100000      * | 0.0253790691494942
    ##        ...      ...               ...    ... .                ...
    ##   [155027]     chrY 57140001-57160000      * |                  0
    ##   [155028]     chrY 57160001-57180000      * |                  0
    ##   [155029]     chrY 57180001-57200000      * |                  0
    ##   [155030]     chrY 57200001-57220000      * |                  0
    ##   [155031]     chrY 57220001-57227415      * |                  0
    ##   -------
    ##   seqinfo: 128 sequences from an unspecified genome

As you can see above the number of ranges in each object is different
(i.e one has a longer genome, weird right?), the reason for such
descripancy is because bw files don’t store chromosome information (like
length) and binning things like I did but there are ways to get around
it. After you solved this problem you can merge the score field of these
Genomic ranges and create your matrix. Here is one approach to make the
bw files similar in length. But you need to tweak it a bit to work for
you.

``` r
higlass_file_prep<-function(bw_path,genome){
  library(rtracklayer)
  library(GenomeInfoDb)
  
  # Importing the proper negspy genome file
  gnm_file<-tempfile(pattern = "",fileext = ".txt")
  if (genome=="hg38"){
    chrom.sizes.path<-"https://raw.githubusercontent.com/omidalam/negspy/master/negspy/data/hg38/chromInfo.txt"
    download.file(url = chrom.sizes.path,destfile = gnm_file)
  }else if(genome=="hg19"){
    chrom.sizes.path<-"https://raw.githubusercontent.com/omidalam/negspy/master/negspy/data/hg19/chromInfo.txt"
    download.file(url = chrom.sizes.path, destfile = gnm_file)
  }
  gnm<-read.table(gnm_file)
  
  # Import the original bw
  bw<-import(bw_path)
  
  # Select the valid chromosomes in the original BW
  valids<-(bw[seqnames(bw)%in%as.character(gnm$V1)])
  valids<-valids[!is.na(valids$score)]
  # Clip the bedgraph based on the negspy genopme
  bedgraph<-tempfile(pattern = 'bg',fileext = ".bedGraph")
  bedgraph_clipped<-tempfile(pattern = 'bg',fileext = ".bedGraph")

  export(valids,con = bedgraph)
  clip_command<-paste("bedClip",bedgraph,gnm_file,bedgraph_clipped)
  system(clip_command)
  
  # Convert the bedGraph to bigwig
  bw_exp_path<-paste(tools::file_path_sans_ext(bw_path),"_negpy_",genome,".bw",sep='')
  bw_command<-paste("bedGraphToBigWig",bedgraph_clipped,gnm_file,bw_exp_path)
  system(bw_command)
}
```

Once you have all these information in proper format you can combine
them like this. But they should have correct number of ranges for each
chromosome. In the below example that I’m adding another column to the
original k562\_tsa data. As you can see it has two columns now.

``` r
k562.tsa$enrichment<-k562.tsa$score^2
head(k562.tsa)
```

    ## GRanges object with 6 ranges and 2 metadata columns:
    ##       seqnames        ranges strand |            score       enrichment
    ##          <Rle>     <IRanges>  <Rle> |        <numeric>        <numeric>
    ##   [1]     chr1   10000-29999      * | 1.41149997711182 1.99233218538666
    ##   [2]     chr1   30000-49999      * | 1.43850004673004 2.06928238444233
    ##   [3]     chr1   50000-69999      * | 1.10230004787445 1.21506539554402
    ##   [4]     chr1   70000-89999      * | 1.06350004673004  1.1310323493948
    ##   [5]     chr1  90000-109999      * | 1.09889996051788 1.20758112322621
    ##   [6]     chr1 110000-129999      * | 1.12170004844666 1.25821099868523
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome

After all these you can extract the all the data like this.

``` r
mcols(k562.tsa)
```

    ## DataFrame with 151570 rows and 2 columns
    ##                   score       enrichment
    ##               <numeric>        <numeric>
    ## 1      1.41149997711182 1.99233218538666
    ## 2      1.43850004673004 2.06928238444233
    ## 3      1.10230004787445 1.21506539554402
    ## 4      1.06350004673004  1.1310323493948
    ## 5      1.09889996051788 1.20758112322621
    ## ...                 ...              ...
    ## 151566 1.03840005397797 1.07827467210144
    ## 151567 1.22580003738403  1.5025857316507
    ## 151568 1.55460000038147 2.41678116118607
    ## 151569  1.7257000207901 2.97804056175495
    ## 151570                0                0
