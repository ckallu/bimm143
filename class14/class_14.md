Class\_14
================

Asthma SNPs Analysis in the MXL (Mexican Ancestry in Los ANgeles, CA) 1000 Patients Sequencing dad
--------------------------------------------------------------------------------------------------

``` r
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
mxl
```

    ##    Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s.
    ## 1                   NA19648 (F)                       A|A ALL, AMR, MXL
    ## 2                   NA19649 (M)                       G|G ALL, AMR, MXL
    ## 3                   NA19651 (F)                       A|A ALL, AMR, MXL
    ## 4                   NA19652 (M)                       G|G ALL, AMR, MXL
    ## 5                   NA19654 (F)                       G|G ALL, AMR, MXL
    ## 6                   NA19655 (M)                       A|G ALL, AMR, MXL
    ## 7                   NA19657 (F)                       A|G ALL, AMR, MXL
    ## 8                   NA19658 (M)                       A|A ALL, AMR, MXL
    ## 9                   NA19661 (M)                       A|G ALL, AMR, MXL
    ## 10                  NA19663 (F)                       A|A ALL, AMR, MXL
    ## 11                  NA19664 (M)                       G|A ALL, AMR, MXL
    ## 12                  NA19669 (F)                       A|A ALL, AMR, MXL
    ## 13                  NA19670 (M)                       A|A ALL, AMR, MXL
    ## 14                  NA19676 (M)                       G|G ALL, AMR, MXL
    ## 15                  NA19678 (F)                       A|A ALL, AMR, MXL
    ## 16                  NA19679 (M)                       A|G ALL, AMR, MXL
    ## 17                  NA19681 (F)                       A|G ALL, AMR, MXL
    ## 18                  NA19682 (M)                       A|G ALL, AMR, MXL
    ## 19                  NA19684 (F)                       A|G ALL, AMR, MXL
    ## 20                  NA19716 (F)                       G|A ALL, AMR, MXL
    ## 21                  NA19717 (M)                       A|G ALL, AMR, MXL
    ## 22                  NA19719 (F)                       G|G ALL, AMR, MXL
    ## 23                  NA19720 (M)                       G|G ALL, AMR, MXL
    ## 24                  NA19722 (F)                       G|A ALL, AMR, MXL
    ## 25                  NA19723 (M)                       G|G ALL, AMR, MXL
    ## 26                  NA19725 (F)                       A|G ALL, AMR, MXL
    ## 27                  NA19726 (M)                       A|A ALL, AMR, MXL
    ## 28                  NA19728 (F)                       A|A ALL, AMR, MXL
    ## 29                  NA19729 (M)                       A|G ALL, AMR, MXL
    ## 30                  NA19731 (F)                       A|A ALL, AMR, MXL
    ## 31                  NA19732 (M)                       A|G ALL, AMR, MXL
    ## 32                  NA19734 (F)                       G|A ALL, AMR, MXL
    ## 33                  NA19735 (M)                       G|G ALL, AMR, MXL
    ## 34                  NA19740 (F)                       A|A ALL, AMR, MXL
    ## 35                  NA19741 (M)                       A|A ALL, AMR, MXL
    ## 36                  NA19746 (F)                       A|A ALL, AMR, MXL
    ## 37                  NA19747 (M)                       G|A ALL, AMR, MXL
    ## 38                  NA19749 (F)                       A|G ALL, AMR, MXL
    ## 39                  NA19750 (M)                       A|G ALL, AMR, MXL
    ## 40                  NA19752 (F)                       A|G ALL, AMR, MXL
    ## 41                  NA19755 (F)                       A|A ALL, AMR, MXL
    ## 42                  NA19756 (M)                       G|A ALL, AMR, MXL
    ## 43                  NA19758 (F)                       A|G ALL, AMR, MXL
    ## 44                  NA19759 (M)                       G|A ALL, AMR, MXL
    ## 45                  NA19761 (F)                       G|A ALL, AMR, MXL
    ## 46                  NA19762 (M)                       A|A ALL, AMR, MXL
    ## 47                  NA19764 (F)                       A|A ALL, AMR, MXL
    ## 48                  NA19770 (F)                       A|G ALL, AMR, MXL
    ## 49                  NA19771 (M)                       A|A ALL, AMR, MXL
    ## 50                  NA19773 (F)                       A|A ALL, AMR, MXL
    ## 51                  NA19774 (M)                       A|G ALL, AMR, MXL
    ## 52                  NA19776 (F)                       A|G ALL, AMR, MXL
    ## 53                  NA19777 (M)                       A|A ALL, AMR, MXL
    ## 54                  NA19779 (F)                       G|A ALL, AMR, MXL
    ## 55                  NA19780 (M)                       A|A ALL, AMR, MXL
    ## 56                  NA19782 (F)                       G|A ALL, AMR, MXL
    ## 57                  NA19783 (M)                       A|G ALL, AMR, MXL
    ## 58                  NA19785 (F)                       A|A ALL, AMR, MXL
    ## 59                  NA19786 (M)                       G|A ALL, AMR, MXL
    ## 60                  NA19788 (F)                       A|G ALL, AMR, MXL
    ## 61                  NA19789 (M)                       G|G ALL, AMR, MXL
    ## 62                  NA19792 (M)                       A|A ALL, AMR, MXL
    ## 63                  NA19794 (F)                       G|A ALL, AMR, MXL
    ## 64                  NA19795 (M)                       A|G ALL, AMR, MXL
    ##    Father Mother
    ## 1       -      -
    ## 2       -      -
    ## 3       -      -
    ## 4       -      -
    ## 5       -      -
    ## 6       -      -
    ## 7       -      -
    ## 8       -      -
    ## 9       -      -
    ## 10      -      -
    ## 11      -      -
    ## 12      -      -
    ## 13      -      -
    ## 14      -      -
    ## 15      -      -
    ## 16      -      -
    ## 17      -      -
    ## 18      -      -
    ## 19      -      -
    ## 20      -      -
    ## 21      -      -
    ## 22      -      -
    ## 23      -      -
    ## 24      -      -
    ## 25      -      -
    ## 26      -      -
    ## 27      -      -
    ## 28      -      -
    ## 29      -      -
    ## 30      -      -
    ## 31      -      -
    ## 32      -      -
    ## 33      -      -
    ## 34      -      -
    ## 35      -      -
    ## 36      -      -
    ## 37      -      -
    ## 38      -      -
    ## 39      -      -
    ## 40      -      -
    ## 41      -      -
    ## 42      -      -
    ## 43      -      -
    ## 44      -      -
    ## 45      -      -
    ## 46      -      -
    ## 47      -      -
    ## 48      -      -
    ## 49      -      -
    ## 50      -      -
    ## 51      -      -
    ## 52      -      -
    ## 53      -      -
    ## 54      -      -
    ## 55      -      -
    ## 56      -      -
    ## 57      -      -
    ## 58      -      -
    ## 59      -      -
    ## 60      -      -
    ## 61      -      -
    ## 62      -      -
    ## 63      -      -
    ## 64      -      -

Find percent of individuals with A|A SNP variant
------------------------------------------------

``` r
genotypes <- round (table (mxl[,2])/nrow(mxl)* 100, 2)
genotypes
```

    ## 
    ##   A|A   A|G   G|A   G|G 
    ## 34.38 32.81 18.75 14.06

There are 34.38% of AA genotype individuals

Interpreting Base Qualities in R
--------------------------------

``` r
#install.packages("seqinr")
#install.packages("gtools")
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Population Gene Expression
--------------------------

``` r
gene_exp <- read.table("https://bioboot.github.io/bimm143_F18/class-material/rs8067378_ENSG00000172057.6.txt")
#gene_exp
summary(gene_exp)
```

    ##      sample     geno          exp        
    ##  HG00096:  1   A/A:108   Min.   : 6.675  
    ##  HG00097:  1   A/G:233   1st Qu.:20.004  
    ##  HG00099:  1   G/G:121   Median :25.116  
    ##  HG00100:  1             Mean   :25.640  
    ##  HG00101:  1             3rd Qu.:30.779  
    ##  HG00102:  1             Max.   :51.518  
    ##  (Other):456

``` r
inds.aa <- gene_exp$geno == "A/A"
inds.ag <- gene_exp$geno == "A/G"
inds.gg <- gene_exp$geno == "G/G"
tab_gene_exp <- boxplot(exp ~ gene_exp$geno, data = gene_exp)
```

![](class_14_files/figure-markdown_github/unnamed-chunk-4-1.png)
