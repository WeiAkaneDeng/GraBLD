#'
#' Remove missing values
#'
#' The function records and removes the missing values from the input data matrix.
#'
#' @param data a numeric data matrix
#' @param NA_value the missing value code in the data, the default is \code{NaN}.
#'    Could also be \code{NA} or \code{-9} commonly used in PLINK data.
#'
#' @return a data matrix of genotypes.
#'
#'

get_NAs = function(data, NA_value = "NaN") {

    printNAs = as.numeric()

    for (i in 1:dim(data)[2]) {
        printNAs = c(printNAs, which(is.na(data[, i])), which(data[,
            i] == NA_value))
    }
    return(printNAs)
}


#'
#' Load univariate regression beta coefficients.
#'
#' The function loads the univariate regression beta coefficients into the workspace.
#'
#' @param trait_name a character for the name of the quantitative trait, assuming the file is named as \code{trait_name_univariate_beta.txt}.
#'
#' @param file_name a character to specify the full directory of the univariate regression coefficient file instead of \code{trait_name_univariate_beta.txt}.
#'
#' @details The univariate regression coefficients are assumed to be computed univariately using: \cr
#'    \code{coef(summary(lm(pheno_data ~ geno[,j])))[2,1]}.
#'    Both genotype data and phenotype data over individuals need to be standardized with
#'    mean = 0 and variance = 1. The data file should have two columns indicating chromosome number
#'    and the regression coefficients without headers.
#'
#' @return a data matrix of univariate regression coefficients.
#'
#' @importFrom utils read.table
#'

load_beta <- function(trait_name, file_name = NULL) {

    if (is.null(trait_name)) {

        beta_value = as.matrix(utils::read.table(paste(file_name),
            header = F))

    } else {

        beta_value = as.matrix(utils::read.table(paste(trait_name,
            "_univariate_beta.txt", sep = ""), header = F))

    }
    return(beta_value)
}


#'
#' Load annotations.
#'
#' The function loads the annotation predictor variables into the workspace.
#'
#' @param annotation_file the directory to a \code{data.frame} of annotation variables used to update the \code{beta} values
#'   through gradient boosted regression tree models. The first column of the matrix must be the SNP IDs and the remaining columns could be
#'   additional annotation information. The SNP IDs must be in the same order as those in \code{beta}.
#'
#' @param pos an integer indicating which columns of the data matrix \code{annotations} is the corresponding consortium value and
#'    additionally which columns should also be included.
#'
#' @details The annotation matrix provides the necessary predictor variables used to update the weights of
#'    polygenic gene score via gradient boosted regression tree. The data.frame should have at least two columns,
#'    the first column is SNP_ID; the rest are the adjusted consortia regression coefficient or summary statistics.
#'    It is recommended to adjust the consortia regression coefficient by the minor allele frequency of the SNP: \cr
#'    \preformatted{
#'    SNP_SD = sqrt(2 * as.numeric(MAF[,5]) * (1 - as.numeric(MAF[,5])))
#'    beta_adj = as.numeric(beta) * SNP_SD
#'    }
#'    For any one trait, at least one column of corresponding adjusted beta from the consortium is required.
#'    For instance, if we work on BMI, at least the adjusted regression coefficient for association with BMI in
#'    a consortium study should be provided. Additional annotations such as related regression coefficients of other
#'    traits, or SNP functional annotations can also be included.
#'
#' @return a data frame of predictor variables that can be used to update SNPs weights.
#'
#' @importFrom data.table fread
#'


# load prediction variables
load_database = function(annotation_file, pos = 2) {
    annotations = as.matrix(data.table::fread(annotation_file,
        sep = "\t", header = T))
    return(annotations[, c(1, pos)])
}


#-------------------------------------------------------------

#'
#' Regression Tree Formula
#'
#' The function generates the formula to be used in regression trees model.
#'
#' @param name a character to indicate the regression coefficients from which complex trait.
#'    It should match the trait name used in the other analyses.
#'
#' @param var_names the names of the variables to be included.
#'
#' @return a formula specifying the response variable (univariate regression coefficients from target population)
#'    and the predictors are used to tune the weights in boosted regression trees model.
#'
#' @importFrom stats as.formula
#'
#'
#' @examples
#' data(annotation_data)
#' var_names <- colnames(annotation_data)
#' formulas = get_formulas(name = 'BMI', normalized = FALSE, var_names = NULL)
#' print(formulas)
#'

get_formulas = function(name, var_names) {

    formulas = paste("target_beta_use ~ ", var_names[1], sep = "")

        if (length(var_names) >= 2) {
        for (i in 2:length(var_names)) {
            formulas = paste(formulas, " + ", var_names[i],
                sep = "")
        }
    }
    return(formulas)
}


#-------------------------------------------------------------

#'
#' Normalize the regression coefficients and annotation data
#'
#' The function process the unviariate beta regression as well as the annotation data matrix and combine
#'    the normalized data used for estimating the optimal boosted regression trees model.
#'
#' @param betas a data matrix of regression coefficients from association analysis in the target population.
#'   The data must contain a minimum two columns, the SNP IDs (denoted by "SNP") and the column with the
#'   regression coefficient ("target_beta"). The chromosome or other information could be present as
#'   additional columns, but will not be used.
#'
#'   The betas were generated from the model: \cr
#'     \code{coef(summary(lm(pheno_data ~ geno[,j])))[2,1]}
#'     and appropriately adjusting for necessary covariates, such as age, sex, and principal component covariates.
#'
#'   Both genotype data and phenotype data over individuals need to be standardized
#'   to have \code{mean} = 0 and \code{variance} = 1.
#'
#' @param annotations a data matrix of annotation variables used to update the \code{beta} values
#'   through gradient boosted regression tree models. Usually, this can be taken from the
#'   summary-level test statistics of matching traits from genome-wide consortia available
#'   online. The data must contain a minimum two columns, the SNP IDs (denoted by "SNP") and a column with the
#'   summary statistics ("annot_primary").
#'   Additional columns for summary statistics or functional information can be included. Any column with names
#'   containing "annot" will be included in the gradient boosted regression tree model, but will not be treated
#'   as the primary annotation, which is used to update the sign of the univariate regression coefficient in the
#'   target population.
#'
#'   For example, the columns can be "SNP", "annot_primary", "annot_trait1", "annot_trait2",
#'   "annot_functional", etc. The "annot_primiary" is usually the consortium univariate regression
#'   coefficient of the same trait.
#'
#' @param abs_effect a vector of integers indicating which columns of the data matrix \code{annotations}
#'   should be used as absolute effect by taking the absolute sign. For example, when only the strength of
#'   the effect rather than the direction of the effect is informative for improving the polygenic score weights.
#'   The default is "NULL", indicating no absolute effect is invoked.
#'
#' @param normalize a logic indicating whether the univariate beta regression coefficients in \code{beta} should be normalized with respect to the consortium values in \code{annotations}.
#'
#' @return a data matrix that can be directly used to estimate the optimal boosted regression trees model.
#'
#' @examples
#' data(annotation_data)
#' get_data_num(betas = univariate_beta, annotations = annotation_data)
#'


get_data_num <- function(betas, annotations, abs_effect = NULL, normalize = FALSE) {

    if (is.null(betas) | is.null(annotations)) {
        stop("Please provide the univariate regression matrix as well as the annotation data matrix.")
    }

    if (sum(grepl("SNP|target_beta", names(betas))) < 2){
        stop("Please make sure the betas has column name SNP and target_beta.")
    }
    if (sum(grepl("SNP|annot_primary", names(annotations))) < 2){
        stop("Please make sure the betas has column name SNP and annot_primary.")
    }

    print(paste("I find", sum(grepl("annot", names(annotations))), "annotation columns and",
                sum(grepl("annot_primary", names(annotations))), "of thme is the primary."))

    if (!is.null(abs_effect)){
        if (sum(grepl("annot_primary", names(annotations)[abs_effect]))==1){
            print("The sign of the primary cannot be ignored, I will proceed with the remaining columns")
        }
        for (i in which(names(annotation) %in% names(annotations)[abs_effect])) {
            norm_data_num[, i] <- abs(data_num[, i])
        }
        }

    suppressMessages(datas <- plyr::join(betas, annotation))

    data_num = datas[, grepl("SNP|annot|beta", names(datas))]

    for (i in which(grepl("annot|beta", names(datas)))) {
        data_num[, i] <- as.numeric(as.character(datas[, i]))
    }

    if (normalize == TRUE) {
        norm_data_num <- data_num
        sign_index <- sign(data_num[, names(data_num) %in% "annot_primary"])
        norm_data_num$target_beta_update <- sign_index*(norm_data_num$target_beta-data_num[, names(data_num) %in% "annot_primary"])

        return(norm_data_num)

    } else {

        return(data_num)
    }
}



#---------------------------------------------------------------------------------

#'
#' Prediction using the Optimal Gradient Boosted Regression Trees using \code{\link[gbm]{gbm}}.
#'
#'
#' The function returns the prediction of the polygenic gene score weights based on the optimal gradient boosted
#'    regression trees model.
#'
#' @param trait_name a character for the name of the quantitative trait, assuming the file is named as \code{trait_name_univariate_beta.txt}.
#'
#' @param betas a data matrix of regression coefficients from association analysis in the target population.
#'   The data must contain a minimum two columns, the SNP IDs (denoted by "SNP") and the column with the
#'   regression coefficient ("target_beta"). The chromosome or other information could be present as
#'   additional columns, but will not be used.
#'
#'   The betas were generated from the model: \cr
#'     \code{coef(summary(lm(pheno_data ~ geno[,j])))[2,1]}
#'     and appropriately adjusting for necessary covariates, such as age, sex, and principal component covariates.
#'
#'   Both genotype data and phenotype data over individuals need to be standardized
#'   to have \code{mean} = 0 and \code{variance} = 1.
#'
#' @param annotations a data matrix of annotation variables used to update the \code{beta} values
#'   through gradient boosted regression tree models. Usually, this can be taken from the
#'   summary-level test statistics of matching traits from genome-wide consortia available
#'   online. The data must contain a minimum two columns, the SNP IDs (denoted by "SNP") and a column with the
#'   summary statistics ("annot_primary").
#'   Additional columns for summary statistics or functional information can be included. Any column with names
#'   containing "annot" will be included in the gradient boosted regression tree model, but will not be treated
#'   as the primary annotation, which is used to update the sign of the univariate regression coefficient in the
#'   target population.
#'
#'   For example, the columns can be "SNP", "annot_primary", "annot_trait1", "annot_trait2",
#'   "annot_functional", etc. The "annot_primiary" is usually the consortium univariate regression
#'   coefficient of the same trait.
#'
#' @param abs_effect a vector of integers indicating which columns of the data matrix \code{annotations}
#'   should be used as absolute effect by taking the absolute sign. For example, when only the strength of
#'   the effect rather than the direction of the effect is informative for improving the polygenic score weights.
#'   The default is "NULL", indicating no absolute effect is invoked.
#'
#' @param normalize a logic indicating whether the univariate beta regression coefficients in \code{beta} should be updated with respect to the consortium values in \code{annotations} in terms of directions of effect.
#'
#' @param steps an integer indicating the current cross-fold
#'
#' @param validation an integer indicating the total number of cross folds. The default and recommended
#'    number of cross-fold is 5.
#'
#' @param interval an integer indicating the number of iterated cycles to calculate the best trees using \code{\link[gbm]{gbm}}.
#'    The default value is 200.
#'
#' @param sig the significance level for including a predictor to build the regression trees model.
#'
#' @param interact_depth an integer for the maximum depth of variable interactions used in \code{\link[gbm]{gbm}}.
#'    See ?\code{\link[gbm]{gbm}}. The default here is 5.
#'
#' @param shrink a shrinkage parameter or the learning rate of the tree models in \code{\link[gbm]{gbm}}. See ?\code{\link[gbm]{gbm}}.
#'    The default is 0.001.
#'
#' @param bag_frac a numeric between 0 and 1, controls the fraction of the training set
#'    observations randomly selected to propose the next tree. See ?\code{\link[gbm]{gbm}}. The default value is 0.5.
#'
#' @param max_tree an integer indicating the total number of trees to fit in \code{\link[gbm]{gbm}}. See ?\code{\link[gbm]{gbm}}.
#'    The default used here is 2000.
#'
#' @param verbose a logic indicating whether the adjusted prediction r-squared for each tested model
#'   with different number of trees should be returned.
#'
#' @param WRITE a logic indicating whether the results of the GraBLD weights should be written to a file with
#'    file name \code{trait_name_gbm_beta.txt}.
#'
#' @return a numeric vector of updated weights with length matching the number of SNPs.
#'
#' @importFrom gbm gbm
#' @importFrom stats as.formula, lm
#' @importFrom utils write.table
#'
#' @examples
#' data(univariate_beta)
#' data(annotation_data)
#' GraB(betas = univariate_beta, annotations = annotation_data,
#'    trait_name = 'BMI', steps = 2, validation = 5)
#'
#' @details
#' For large datasets, it is recommended to run from the command line with \cr
#'
#' \preformatted{
#' validation=5
#' for (( i = 1; i <= $validation; i++))
#' do
#' Rscript calculate_gbm.R geno_data
#'    trait_name annotations_file pos ${i}
#'    $validation interaction_depth
#'    shrinkage_parameter bag_fraction
#'    maximum_tree &
#' done
#' }
#' where the R script \code{calculate_gbm.R} might look something like this, while additional options can be
#' added to the argument list: \cr
#'
#' \preformatted{#!/bin/sh
#' rm(list = ls())
#' library('GraBLD')
#' args = (commandArgs(TRUE))
#' geno_data = args[1]
#'  trait_name = args[2]
#'  annotations_file = args[3]
#'  pos = eval(parse(text=args[4]))
#'  steps = eval(parse(text=args[5]))
#'  validation = eval(parse(text=args[6]))
#'  p1 = eval(parse(text=args[7]))
#'  p2 = eval(parse(text=args[8]))
#'  p3 = eval(parse(text=args[9]))
#'  p4 = eval(parse(text=args[10]))
#'  betas = load_beta(trait_name)
#'  annotation = load_database(annotations_file, pos = 2:3)
#'  geno <- load_geno(geno_data)
#'  GraB(betas = betas, annotations = annotation,
#'    trait_name = trait_name, steps = steps, validation = validation,
#'    interval = 200, sig = 1e-05, interact_depth = p1, shrink = p2,
#'    bag_frac = p3, max_tree = p4, WRITE = TRUE)}
#'
#'
#' @references Greg Ridgeway with contributions from others (2015). gbm: Generalized Boosted Regression Models. R package version 2.1.1. \url{https://CRAN.R-project.org/package=gbm}
#'
#'Guillaume Pare, Shihong Mao, Wei Q Deng (2017)
#' A machine-learning heuristic to improve gene score prediction of
#' polygenic traits Short title: Machine-learning boosted gene scores,
#' \emph{bioRxiv 107409}; doi: \url{https://doi.org/10.1101/107409};
#' \url{http://www.biorxiv.org/content/early/2017/02/09/107409}
#'


GraB <- function(betas, annotations, abs_effect = NULL, trait_name = NULL,
    steps = 1, validation = 5, verbose = FALSE, interval = 200,
    sig = 1e-05, interact_depth = 5, shrink = 0.001, bag_frac = 0.5,
    max_tree = 2000, WRITE = FALSE) {


    if (is.null(betas) | is.null(annotations)) {
        stop("Please provide the univariate regression matrix as well as the annotation data matrix.")
    }

    data_num_pr = get_data_num(betas = betas, annotations = annotations, abs_effect = abs_effect, normalize = normalize)

    sign_index <- sign(data_num_pr$target_beta)

    if ("target_beta_update" %in% names(data_num_pr)){
        data_num_pr$target_beta_use <-  data_num_pr$target_beta_update
    } else {
    data_num_pr$target_beta_use <- data_num_pr$target_beta
    }

    nb_SNPs = floor(dim(data_num_pr)[1]/validation)
    results = as.numeric()

    starting = (steps - 1) * nb_SNPs + 1
    ending = steps * nb_SNPs
    if (steps == validation) {
        ending = dim(data_num_pr)[1]
    }

    var_names <- colnames(data_num_pr)[grepl("annot", colnames(data_num_pr))]
    formulas = get_formulas(name = trait_name, var_names = var_names)

    data_predict <- as.data.frame(data_num_pr[c(starting:ending),
        ])
    data_train = as.data.frame(data_num_pr[-c(starting:ending),
        ])


    lm1_variables <- suppressWarnings(stats::lm(stats::as.formula(formulas),
        data = data_train))
    var_index <- which(summary(lm1_variables)[["coefficients"]][,
        4] < sig)
    if (var_index[1] == 0) {
        print("The significance threshold selected does not support additional variables to be included, only the primary will be included in the final model.")
        formula_final = paste("target_beta_use ~ annot_primary")
            } else {
        formula_final = paste("target_beta_use ~ ", paste(var_names[var_index[-1]-1], collapse = "+"))
            }

    gbm1_diff2c <- gbm::gbm(stats::as.formula(formula_final),
        data = data_train, distribution = "gaussian", interaction.depth = interact_depth,
        shrinkage = shrink, bag.fraction = bag_frac, n.trees = max_tree,
        n.cores = 1)

    result = as.numeric()

    cycles = max_tree/interval

    for (i in 1:cycles) {
        tree = i * interval
        gbm1_diff2c_predict <- gbm::predict.gbm(gbm1_diff2c,
            data_predict, n.trees = tree)
        result = c(result, summary(stats::lm(data_predict$target_beta_use ~ gbm1_diff2c_predict))$adj.r.squared)
    }

    results = cbind(results, result)

    GRS_adj_GIANT_all = data_predict$target_beta[starting:ending]
    sign_index_predict = sign_index[starting:ending]
    SNP_name <- as.character(data_num_pr$SNP[starting:ending])

    best_nb_tree = which(max(result) == result) * interval

    gbm1 <- gbm::gbm(stats::as.formula(formulas), data = data_train,
        distribution = "gaussian", interaction.depth = interact_depth,
        shrinkage = shrink, bag.fraction = bag_frac, n.trees = best_nb_tree,
        n.cores = 1)
    gbm1_predict_all <- gbm::predict.gbm(gbm1, data_predict,
        n.trees = best_nb_tree) * sign_index_predict + GRS_adj_GIANT_all


    if (verbose == TRUE) {

        if (WRITE == TRUE) {

            utils::write.table(SNP_name, paste(trait_name,
                "_SNPs.txt", sep = ""), col.names = F, row.names = F,
                quote = F, sep = "\t")
            utils::write.table(gbm1_predict_all, paste(trait_name,
                "_gbm_", steps, ".txt", sep = ""), col.names = F,
                row.names = F, quote = F, sep = "\t")
            utils::write.table(result, paste(trait_name, "_selectTrees.txt",
                sep = ""), col.names = F, row.names = F, quote = F,
                sep = "\t")

        } else {

            return(list(weights = gbm1_predict_all, selectTrees = results,
                SNP_name = SNP_name))

        }

    } else {

        if (WRITE == TRUE) {

            utils::write.table(SNP_name, paste(trait_name,
                "_SNPs.txt", sep = ""), col.names = F, row.names = F,
                quote = F, sep = "\t")
            utils::write.table(gbm1_predict_all, paste(trait_name,
                "_gbm_", steps, ".txt", sep = ""), col.names = F,
                row.names = F, quote = F, sep = "\t")

        } else {

            return(weights = gbm1_predict_all)

        }
    }
}

#-------------------------------------------------------------


