##'
#' Loading the genotype data.
#'
#' The function automatically loads the genotype data matrix assuming
#' each row is an individual and each column is the genotype of a
#' biallelic SNP. If standard PLINK .raw file was provided,
#' the function returns only the genotype matrix by removing
#' the first 6 columns.
#'
#' @param source_data the name of the file which the genotype data are to be read from.
#'  Could be an absolute path to the file or the name of the file assuming in the
#'  present directory \code{\link{getwd}()}. In this case, the file name should be :
#'    \code{source_data.raw}.
#'    or
#'    \code{source_data.bed}, \code{source_data.bim}, and \code{source_data.fam}.
#'
#' @param PLINKbin a logic indicating whether the supplied file is of the .raw format from PLINK, if not, the first six columns will not be removed and all columns will be read in.
#'
#'
#' @return a data matrix of genotypes.
#'
#' @importFrom data.table fread
#' @import BGData
#' @import BEDMatrix
#' @importFrom utils installed.packages install.packages
#' @export load_geno


load_geno <- function(source_data, PLINKbin = TRUE) {

    if (PLINKbin) {

        if("BEDMatrix" %in% rownames(utils::installed.packages()) == FALSE) {
            print("BEDMatrix not installed, trying to intall now ...")
            utils::install.packages("BEDMatrix", repos='http://cran.us.r-project.org')
        }

        if("BGData" %in% rownames(utils::installed.packages()) == FALSE) {
            print("BGData not installed, trying to intall now ...")
            utils::install.packages("BGData", repos='http://cran.us.r-project.org', dependencies=T)
        }

    bedFiles <- BEDMatrix::BEDMatrix(source_data)
    bg <- BGData::as.BGData(bedFiles, alternatePhenotypeFile = paste(paste0(source_data), ".fam",sep=""))
    Data <- BGData::geno(bg)

    } else {
        raw_Data = as.matrix(data.table::fread(paste(source_data,".raw", sep = ""), header = T))
        Data <- matrix(as.numeric(raw_Data[, 7:dim(raw_Data)[2]]),
                       nrow = dim(raw_Data)[1], byrow = FALSE)
        colnames(Data) <- names(raw_Data[, 7:dim(raw_Data)[2]])
    }

    return(Data)
}


#-------------------------------------------------------------------------------------

#'
#' Standardizing the genotype data.
#'
#' The function standardizes a genotype data such that for each SNP the mean is 0 and standard
#'    deviation is 1.
#'
#' @param geno_raw the genotype matrix assuming each row is the individual and each column is the SNP.
#'
#' @param max_size an integer for the maximum size of SNPs to be standardized at a time, the default is 100000.
#'    This option can be ignored for data with less than 1 million SNPs.
#'
#' @param NAval the missing value used in the output data matrix. The default is NA, for use with
#'    PLINK, set \code{NAval} = -9.
#'
#' @return a data matrix of standardized genotypes.
#'
#' @importFrom stats sd complete.cases
#'
#' @examples
#' data(geno)
#' norm_geno <- full_normal_geno(geno_raw = geno)
#' @export full_normal_geno


full_normal_geno <- function(geno_raw, max_size = 1e+05, NAval = NA) {

    Geno_norm <- as.numeric()
    cycles <- floor(dim(geno_raw)[2]/max_size) + 1

    for (i in 1:cycles) {
        starting = (i - 1) * max_size + 1
        ending = i * max_size
        if (i == cycles) {
            ending = dim(geno_raw)[2]
        }
        raw_Data = geno_raw[, starting:ending]

        if (sum(stats::complete.cases(raw_Data))<dim(raw_Data)[1]){
        ## imputing missing value by mean
        for (g in which(!stats::complete.cases(raw_Data))) {
            raw_Data[g, which(is.na(raw_Data[g,]))] = mean(raw_Data[g,
                ], na.rm = TRUE)
        }
        }

        Geno = raw_Data

        n = dim(Geno)[1]
        m = dim(Geno)[2]

        for (i in 1:m) {
            Geno[, i] = (Geno[, i] - mean(Geno[, i], na.rm=T))/min(stats::sd(Geno[,
                i], na.rm=T), 10000)
            ## just in case bad values,
            ## the variance will make this observation negligible
        }

        Geno_norm = cbind(Geno_norm, Geno)
        #Geno_norm = Geno
    }
    Geno_norm[is.na(Geno_norm)] = NAval
    return(Geno_norm)
}


#-------------------------------------------------------------------------------------

#'
#' Regional LD adjustment
#'
#' The function calculates the LD adjustment of each SNP for
#' one given size of genotype data block.
#'
#' @param geno_data the standardized genotype matrix assuming each row is the individual and each column is the SNP.
#'
#' @param size an integer for the number of SNPs that should
#' be included in a block. Usually, for genome-wide datasets,
#' there are about 2 million SNPs and +/- 300 SNPs is roughly
#' equivalent to 1Mb physical distance, thus 300 is set as
#' the default size.
#'
#' @param position a character indicating the position of the SNP in a block of SNPs, one of\cr
#'    \code{start} (the start of the block), \code{mid} (the middle of the block),
#'    \code{end} (the end of the block), or \code{together} (combining SNPs in a few blocks),
#'    \code{short} (in a short block), or \code{all} (using all regions).
#'
#' @return a numeric vector of LD adjustments with the same length as the number of SNPs in the genotype matrix.
#'
#' @references Pare, Guillaume, Shihong Mao, and Wei Q. Deng.
#' A method to estimate the contribution of regional genetic
#' associations to complex traits from summary association statistics.
#' \emph{Scientific reports} \strong{6} (2016): 27644.
#'
#' @examples
#' data(geno)
#' norm_geno <- full_normal_geno(geno_raw = geno[,-c(1:6)])
#' regLD <- regionalLD(geno_data = norm_geno, position = 'mid', size = 300)
#' print(regLD)
#'
#' @export regionalLD
#'


regionalLD <- function(geno_data, position = "start", size = 300) {

    if (size < 1) {
        stop("The size of the region must be greater than 1.")
    }
    geno_data <- as.matrix(geno_data)
    size = as.integer(size)

    if (position == "start") {

        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        for (i in 1:(size + 1)) {
            value = sum((R2_matrix[i, 1:(i + size)])^2)
            LDadj = c(LDadj, value)
        }
    }

    if (position == "mid") {

        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        if ((snps - size) < (size + 2)) {
            stop("If the number of SNPs is small, please adjust the size accordingly.")
        }

        for (i in (size + 2):(snps - size)) {
            value = sum((R2_matrix[i, (i - size):(i + size)])^2)
            LDadj = c(LDadj, value)
        }
    }

    if (position == "end") {

        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        for (i in (snps - size + 1):snps) {
            value = sum((R2_matrix[i, (i - size):snps])^2)
            LDadj = c(LDadj, value)
        }
    }

    if (position == "together") {
        size = as.integer(size)
        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        for (i in 1:(size + 1)) {
            value = sum((R2_matrix[i, 1:(1 + size)])^2)
            LDadj = c(LDadj, value)
        }

        if ((snps - size) > (size + 2)) {
            for (i in (size + 2):(snps - size)) {
                value = sum((R2_matrix[i, (i - size):(i + size)])^2)
                LDadj = c(LDadj, value)
            }

            starting = snps - size + 1
            if (snps <= 2 * size) {
                starting = size + 1
            }
        } else {
            starting = size + 2
        }
        for (i in starting:snps) {
            value = sum((R2_matrix[i, (i - size):snps])^2)
            LDadj = c(LDadj, value)
        }

    }

    if (position == "short") {

        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        for (i in 1:(size + 1)) {
            ending = i + size
            if (ending > snps) {
                ending = snps
            }
            value = sum((R2_matrix[i, 1:ending])^2)
            LDadj = c(LDadj, value)
        }

        starting = size + 2
        for (i in starting:snps) {
            value = sum((R2_matrix[i, (i - size):snps])^2)
            LDadj = c(LDadj, value)
        }
    }


    if (position == "all") {

        individual = dim(geno_data)[1]
        R2_matrix = t(geno_data) %*% geno_data/individual
        snps = dim(R2_matrix)[2]
        LDadj = as.numeric()

        for (i in 1:snps) {
            value = sum((R2_matrix[i, 1:snps])^2)
            LDadj = c(LDadj, value)
        }


    }

    return(LDadj)
}



#'
#' Genome-wide LD adjustment
#'
#' The function compute the LD adjustment on the entire
#' genotype data matrix by first dividing them into SNP blocks.
#'
#' @param geno_data the standardized genotype matrix assuming each row is the individual and each column is the SNP.
#'
#' @param size an integer for the number of SNPs that should be
#' included in a block. Usually, for genome-wide datasets,
#' there are about 2 million SNPs and +/- 300 SNPs
#' is roughly equivalent to 1Mb physical distance,
#' thus 300 is set as the default size.
#'
#' @return a numeric vector of LD adjustments for all SNPs in the genotype matrix.
#'
#' @references Pare, Guillaume, Shihong Mao, and Wei Q. Deng.
#' A method to estimate the contribution of regional genetic
#' associations to complex traits from summary association statistics.
#' \emph{Scientific reports} \strong{6} (2016): 27644.
#'
#' @export genomewideLD
#'

genomewideLD <- function(geno_data, size = 300) {

    if (size < 1) {
        stop("The size of the region must be greater than 1.")
    }
    geno_data <- as.matrix(geno_data)
    size = as.integer(size)

    bin = 6 * as.integer(size)
    starting = 1
    ending = starting + bin
    geno_size = dim(geno_data)[2]
    LDadj = as.numeric()
    if (size >= geno_size) {
        LDadj = regionalLD(geno_data, size = size, position = "all")
    } else if (2 * size >= geno_size) {
        LDadj = regionalLD(geno_data, size = size, position = "short")
    } else if (ending >= geno_size & 2 * size <= geno_size) {
        LDadj = regionalLD(geno_data, size = size, position = "together")
    } else {
        while (ending < geno_size) {
            geno_data_inter = geno_data[, starting:ending]
            if (starting == 1) {
                LDadj = c(LDadj, regionalLD(geno_data_inter,
                  size = size, position = "start"))
            }
            LDadj = c(LDadj, regionalLD(geno_data_inter, size = size,
                position = "mid"))

            starting = ending - size * 2
            ending = starting + bin
            if (ending > geno_size) {
                ending = geno_size
                geno_data_inter = geno_data[, starting:ending]
                LDadj = c(LDadj, regionalLD(geno_data_inter,
                  size = size, position = "mid"))
                LDadj = c(LDadj, regionalLD(geno_data_inter,
                  size = size, position = "end"))
            }
        }
    }
	LDadj[which(LDadj < 1)] = 1 #Monomorphic SNPs
    return(LDadj)
}
