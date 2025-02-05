# GraBLD

Gradient boosted and LD adjusted (GraBLD) is an R based software package that performs polygenic traits prediction using gene scores. The package is best suited for moderate to large sample sizes and takes advantage of the consortium summary statistics to boost the prediction performance. 

## 1.Quick Start  

The current release is: version 1.2. See "release" tab. 

To install our R package “GraBLD”, you can either run in R directly:

    install.packages("devtools") # if you have not installed already
    devtools::install_github("WeiAkaneDeng/GraBLD")

To load the library in R, simply run in R:

    library(“GraBLD”)
    
For complex traits that are highly polygenic, the best way to proceed for a genome-wide PRS is to use the command line options as demonstrated in the tips below.

## 2. Tips for using GraBLD on large datasets

The package contains two main functions, one is to calculate the LD adjustments and the other to calculate the gradient boosted polygenic score weights. The two are combined in the end to give the actual weights of the polygenic score. Though the method works for virtually any size of data (large or small sample size, and large or small number of SNPs), in practice, the performance will depend on the underlying effect size distribution of the phenotype of interest.  

A very important point here is to calculate LD adjustments using only the SNPs that will be entering the gene scores (that will be used for prediction), otherwise you risk over penalizing each SNP. Naturally, we do not recommend using existing databases, such as that for LDscore for this purpose. For the purpose of predictive PRS, it is best to calculate the LD adjustment in-sample rather than out-of-sample. 

Logistically, the first step is to construct the annotation files that will be used to update the gradient boosted weights. Then, compare this file with the genotype data in-sample, we can extract a list of SNPs for which LD adjustments and gradient boosted weights will be calculated on. The final step involves combining these two sets of weights for polygenic risk score construction using *GraBLD.score()*.


### 2.1 A note on annotation file

The main driver of the gain performance is the summary statistics of a matching trait from large consortium, which is the second column of the annotation file. But the GraB algorithm does not need to be restricted to just one set of summary statistics, those from similar classes can be included. For example, for prediction of BMI, the second column can be the summary statistics for BMI from [GIANT](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files), and the third column can be the summary statistics for T2D from [DIAGRAM](https://www.diagram-consortium.org/), and fourth can be HDL from [GLGC](http://lipidgenetics.org). Additional functional information can also be included to updates the weights of each SNP included in the gene score.

An important issue in using the external summary statistics with those from target data is, the phenotype on which the summary stats are calculated should be as similar as possible, both in terms of range and transformation applied. For example, if the BMI GWAS summary data from GIANT performed log transformation and was adjusted for sex, age, and waist, then the BMI in target data should be processed similarly to obtain target specific summary statistics. 

Note that the computation of the target specific summary statistics must be done after sample splitting (if no additional data available). For prediction within sample, I recommend 50/50 split of the target data where each 50% portion is used to generate a set of GWAS summary. Then we simultanesouly (but separately) update the two sets of GWAS summary from target by the annotation (external GWAS). The GraBLD weights are then available to each 50/50 without sample overlap and can be used to perform prediction in the entire target sample without obvious overfitting issues. However, if the goal is to validate the PRS, then any random sample splitting is possible, provided that the samples used to generate target specific summary statistics are not part of the prediction/testing sample.



### 2.2 LD adjustments

For large datasets in PLINK binary format (data_name.bim, data_name.fam, data_name.bed), it is recommended to run from the command line with each chromosome separately: 

    for((i = 1; i <= chr; i++))
    do
    Rscript PerformLDadj.R size data_name ${i} &
    done

where the R script "PerformLDadj.R" might look something like this, while additional options can be added to the argument list: 
    
    #!/bin/sh
    rm(list = ls())
    library('GraBLD')
    args = (commandArgs(TRUE))
    size = eval(parse(text=args[1]))
    source_data = args[2]
    chr = eval(parse(text=args[3]))
    library("BGData")
    library("BEDMatrix")
    LD_OUT <- LDadj(source_data = source_data, chr = chr, size = size, write = TRUE)
    print("Analysis completed for chromosome", chr, ".")


### 2.3 Gradient boosted weights

For large datasets, it is recommended to run from the command line with 
    
    source_data="chr1" with chr1.bim, chr1.bed, chr1.fam
    trait_name="BMI" # name of the trait
    annotations_file="annotation.txt" # or the actual file name of the annotation data
    validation=5 # or value of your choosing
    interaction_depth=5 # or value of your choosing
    shrinkage_parameter=0.01 # or value of your choosing
    bag_fraction=0.5 # or value of your choosing
    maximum_tree=2000 # or value of your choosing
    
    for (( i = 1; i <= $validation; i++))
    do
    Rscript calculate_gbm.R $source_data
       $trait_name $annotations_file ${i}
       $validation $interaction_depth
       $shrinkage_parameter $bag_fraction
       $maximum_tree &
    done
    
where the R script calculate_gbm.R might look something like this, while additional options can be added to the argument list: 

    #!/bin/sh
    rm(list = ls())
    library('GraBLD')
    args = (commandArgs(TRUE))
    source_data = args[1]
    trait_name = args[2]
    annotations_file = args[3]
    steps = eval(parse(text=args[4]))
    validation = eval(parse(text=args[5]))
    p1 = eval(parse(text=args[6]))
    p2 = eval(parse(text=args[7]))
    p3 = eval(parse(text=args[8]))
    p4 = eval(parse(text=args[9]))
    betas = load_beta(trait_name)
    annotation = load_database(annotations_file, pos = 2:3) # taking the 2nd and 3rd columns of "annotations_file"
    geno <- load_geno(source_data)
    GraB(betas = betas, annotations = annotation,
      trait_name = trait_name, steps = steps, validation = validation,
      interval = 200, sig = 1e-05, interact_depth = p1, shrink = p2,
      bag_frac = p3, max_tree = p4, WRITE = TRUE)



## 3.Combining the two 

The two steps above produce separate files, one for the LD adjustments and the other the boosted weights, and they can be combined by the function "GraBLD.score" by the SNP column.

    LD_val <- read.table("OUTPUTS_FROM_STEP1.txt") ## if each chromosome is computed separately, you will need to append them into a single file 
    gbm_val <- read.table("OUTPUTS_FROM_STEP2.txt")
    gs <- GraBLD.score(geno_raw = YOUR_GENO_DATA, LDadjVal = LD_val, gbmVal = gbm_val, Pheno = NULL)

Notice that if "Pheno" is supplied, both the polygenic gene score as well as the prediction R-squared (adjusted) are returned (for quantitative trait), otherwise only the polygenic gene score is returned.
