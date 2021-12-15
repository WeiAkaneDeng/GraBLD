#'
#' Gradient Boosted and LD adjusted Prediction.
#'
#' The function returns the prediction r-squared on the target population using polygenic gene score based on the GraBLD heuristic.
#'
#' @param source_data the name of the file which the genotype data are to be read from.
#'    Each row of the matrix appears as one line of the file. Could be an absolute path to the
#'    file or the name of the file assuming in the present directory \code{\link{getwd}()}.
#'
#' @param PLINK a logic indicating whether the supplied file is of the .raw format from PLINK, if not, the first six
#'    columns will not be removed and all columns will be read in.
#'
#' @param chr an integer indicating the maximum number of chromosomes to be read in, this option is used in combination with
#'    \code{source_data}, to perform analysis by each chromosome. In this case, the file name should follow:
#'    ``source_data_i.raw'' for all \eqn{i \le chr}. For example, \code{Geno_Matrix_23.raw}.
#'
#' @param SNP_names a vector of characters for the names of SNPs matching the corresponding genotype data, these are used in the output of GraBLD weights.
#'
#' @param geno_raw the genotype matrix assuming each row is the individual and each column is the SNP. Can be skipped if \code{source_data} was provided.
#'
#' @param max_size an integer for the maximum size of SNPs to be standardized at a time, the default is 100000.
#'    This option can be ignored for data with less than 1 million SNPs.
#'
#' @param NAval the missing value used in the output data matrix. The default is NA, for use with PLINK, set \code{NAval} = -9.
#'
#' @param Pheno a numeric vector of quantitative traits with the same length as the number of rows in the genotype matrix.
#'
#' @param trait_name a character indicating the name of the quantitative trait.
#'
#' @param WRITE a logic indicating whether the results of the GraBLD weights should be written to a file with file name \code{trait_name_gbm_beta.txt}.
#'
#' @param LDadjVal a data.frame of LD adjusted weights. This can be taken from the output of \code{\link{LDadj}()}. The column names must be "SNP" and "LDadj".
#'
#' @param gbmVal a data.frame of gradient boosted weights. This can be taken from the output of \code{\link{GraB}()}. The column names must be "SNP" and "GraB".
#'
#' @return if \code{Pheno} is supplied, both the polygenic gene score as well as the prediction R-squared (adjusted) are returned, otherwise only the polygenic gene score is returned.
#'
#' @importFrom data.table fread
#' @importFrom gbm gbm
#' @importFrom plyr join
#' @importFrom stats lm sd as.formula
#' @importFrom utils write.table
#' @export GraBLD.score
#'
#' @examples
#' data(geno)
#' data(univariate_beta)
#' data(annotation_data)
#' LD_val <- LDadj(geno_raw = geno[grepl("rs", names(geno))],
#' SNP_names = names(geno)[grepl("rs", names(geno))],
#' chr = 1, size = 200)
#' gbm_val <- list()
#' for (j in 1:5){
#' gbm_val[[j]] <- GraB(betas = univariate_beta,
#' annotations = annotation_data, trait_name = 'BMI',
#' steps = j, validation = 5)
#' }
#' data(BMI)
#' gs <- GraBLD.score(geno_raw = geno[grepl("rs", names(geno))],
#'  SNP_names = names(geno)[grepl("rs", names(geno))],
#'  LDadjVal=LD_val, gbmVal = do.call(rbind, gbm_val),
#'  trait_name='BMI', Pheno = BMI[,1])
#'
#' @references Paré, G., Mao, S. & Deng, W.Q. A machine-learning heuristic to improve gene score prediction of polygenic traits. \emph{Scientific reports} \strong{7}, 12665 (2017). \doi{10.1038/s41598-017-13056-1}.
#'

GraBLD.score <- function(source_data = NULL, SNP_names = NULL,
                         chr = NULL, geno_raw, PLINK = TRUE,
                         max_size = 1e+05, NAval = NA,
                         Pheno = NULL, LDadjVal, gbmVal,
                         trait_name = NULL, WRITE = FALSE) {

    if (is.null(LDadjVal) | is.null(gbmVal)) {
        stop("Please provide the LD adjustments and the gbm tuned weights.")
    }

    if (dim(gbmVal)[1]!=dim(LDadjVal)[1]) {
        stop("Please ensure the length of the LD adjustments matches that of the gbm tuned weights.")
    } else {
        suppressMessages(datas <- plyr::join(gbmVal, LDadjVal))
        gbm_adj = datas$GraB/datas$LDadj
    }

    geno_norm = as.numeric()

    if (is.null(chr) & is.null(source_data)) {

        geno_norm = full_normal_geno(geno_raw, NAval = NAval,
            max_size = max_size)

    } else {

        for (i in 1:chr) {
            geno_pre_chr = load_geno(source_data, PLINKbin = PLINK)
            geno_norm_chr = full_normal_geno(geno_pre_chr,
                NAval = NAval, max_size = max_size)
            geno_norm = cbind(geno_norm, geno_norm_chr)
        }
    }

    if (dim(geno_norm)[2]!= dim(gbmVal)[1]) {
        stop("Please ensure the number of SNPs matches the length of the GraBLD weights.")
    }

    suppressMessages(datas <- plyr::join(gbmVal, LDadjVal))

    if (is.null(SNP_names)){
     gbm_adj = (datas$GraB/datas$LDadj)
      warning("SNP names not supplied, assuming they have been matched to the same order as the GraBLD weights.")
    } else {
    gbm_adj = (datas$GraB/datas$LDadj)[match(SNP_names, datas$SNP)]
    }


    gbm_adj[is.na(gbm_adj)] <- 0
    gene_score = as.matrix(geno_norm) %*% as.vector(gbm_adj)

    if (WRITE == TRUE) {

        gbm_adj_title = c("gbm_beta", gbm_adj)

        if (is.null(SNP_names)) {
            if (length(SNP_names) != length(gbm_adj_title)) {
                warnings("the number of SNPs supplied does not match the length of the GraBLD weights,
               only the GraBLD weights are outputted.")
                utils::write.table(gbm_adj_title, paste(trait_name,
                  "_gbm_beta.txt", sep = ""), col.names = F,
                  row.names = F, quote = F, sep = "\t")
            } else {
                utils::write.table(cbind(SNP_names, gbm_adj_title),
                  paste(trait_name, "_gbm_beta.txt", sep = ""),
                  col.names = F, row.names = F, quote = F,
                  sep = "\t")
            }
        } else {
            utils::write.table(gbm_adj_title, paste(trait_name,
                "_gbm_beta.txt", sep = ""), col.names = F,
                row.names = F, quote = F, sep = "\t")
        }
    }

    if (is.null(Pheno)) {

        return(GraBLD.gs = gene_score)

    } else {
        Pheno_norm = (Pheno - mean(Pheno, na.rm = T))/stats::sd(Pheno,
            na.rm = T)
        adjR2 = summary(stats::lm(Pheno_norm ~ gene_score))$adj.r.squared
        return(list(GraBLD.gs = gene_score, adj.R2 = adjR2))
    }
}

