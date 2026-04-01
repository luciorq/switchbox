#' Training Expression Data
#'
#' A numeric matrix of gene expression data for 78 breast cancer patients and
#' 70 genes from the MammaPrint assay. Samples are in columns and genes in
#' rows.
#'
#' @format A numeric matrix with 70 rows (genes) and 78 columns (samples).
#'
#' @source
#' Glas et al., BMC Genomics, 2006. Originally obtained from the
#' MammaPrintData Bioconductor package as described in Marchionni et al.,
#' BMC Genomics, 2013.
#'
#' @examples
#' data(trainingData)
#' dim(matTraining)
#' table(trainingGroup)
"matTraining"


#' Training Class Labels
#'
#' A factor with two levels (`"bad"` and `"good"`) indicating the prognosis
#' group for each of the 78 training samples.
#'
#' @format A factor of length 78 with levels `"bad"` and `"good"`.
#'
#' @source
#' Glas et al., BMC Genomics, 2006.
#'
#' @examples
#' data(trainingData)
#' table(trainingGroup)
"trainingGroup"


#' Test Expression Data
#'
#' A numeric matrix of gene expression data for 307 breast cancer patients
#' and 70 genes from the MammaPrint validation cohort. Samples are in columns
#' and genes in rows.
#'
#' @format A numeric matrix with 70 rows (genes) and 307 columns (samples).
#'
#' @source
#' Buyse et al., JNCI, 2006. Originally obtained from the MammaPrintData
#' Bioconductor package as described in Marchionni et al., BMC Genomics, 2013.
#'
#' @examples
#' data(testingData)
#' dim(matTesting)
#' table(testingGroup)
"matTesting"


#' Test Class Labels
#'
#' A factor with two levels (`"bad"` and `"good"`) indicating the prognosis
#' group for each of the 307 test samples.
#'
#' @format A factor of length 307 with levels `"bad"` and `"good"`.
#'
#' @source
#' Buyse et al., JNCI, 2006.
#'
#' @examples
#' data(testingData)
#' table(testingGroup)
"testingGroup"
