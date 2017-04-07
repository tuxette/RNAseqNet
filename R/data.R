
##########################
## Lung data (X)
#############################
#'
#' RNA-seq expression from lung tissue
#' 
#' @description
#' This data set is a small subset of the full data set from GTEx.  It contains 
#' RNA-seq expression measured from lung tissue. The RNA-seq expression have been normalized with the TMM method.
#' 
#' @name lung
#' @docType data
#' @format a data frame with 221 rows and 100 variables (genes)
#' 
#' @example 
#' data(lung)
#' 
#' @source The raw data were download from \url{https://gtexportal.org/}. 
#' The TMM normalisation of RNA-seq expression was realized with the R package edgeR.
#' 
NULL


##########################
## Thyroid data (Y)
#############################
#'
#' RNA-seq expression from thyroid tissue
#' 
#' @description
#' This data set is a small subset of the full data set from GTEx.  It contains 
#' RNA-seq expression measured from thyroid tissue. The RNA-seq expression have been normalized with the TMM method.
#' 
#' @docType data
#' @name thyroid
#' @format a data frame with 221 rows and 50 variables (genes)
#' 
#' @example 
#' data(thyroid)
#' 
#' @source The raw data were download from \url{https://gtexportal.org/}.
#'  The TMM normalisation of RNA-seq expression was realized with the R package edgeR.
#'  
NULL