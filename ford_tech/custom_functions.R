#' Process BC or F2 mapping population to a mapmaker compatible format.
#' @param x A data frame of samples as rows and markers as columns.
#' @param type A string indicating the type of mapping population; currently BC
#'  or F2.
#' @param homo1 Coding used for homozygote 1.
#' @param homo2 Coding used for homozygote 2.
#' @param het Coding used for heterozygote 1.
#' @param geno_start An integer indicating the column index for the beginning
#' of the markers.
#' @param trait_id A character vector indicating the column names for traits
#' if any.
#' @param ext A string for file extension; either .raw or .txt.
proc_mapmaker <- function(x,
                          type = c("f2 backcross", "f2 intercross"),
                          homo1 = 1,
                          homo2 = 3,
                          het = 2,
                          geno_start = 1,
                          trait_id = NULL,
                          out_name = "mydata",
                          ext = c('.raw', '.txt')) {

  type <- match.arg(type)
  ext <- match.arg(ext)

  # Remove sample meta data from the imported genotype data
  if (geno_start > 1) {

    geno_mat <- x[, c(geno_start:ncol(x))]

  } else if (geno_start  == 1) {

    geno_mat <- x

  }

  # Process traits if present
  if (!is.null(trait_id)) {

    # Subset only marker data without traits
    geno_mat <- geno_mat[, !colnames(geno_mat) %in% trait_id]

    # Subset only traits data without marker data
    trait <- x[, colnames(x) %in% trait_id]

    ntrait <- length(trait_id) # Get number of traits

    trait_ns <- paste0('*', trait_id) # Add asterisk to trait names

    res <- c() # Empty vector to store trait values

    # Convert each trait data to MapMaker format
    for (i in seq_len(ntrait)) {

      t1 <- paste(trait[, i], collapse = " ")

      res[i] <- paste(trait_ns[i], t1)

    }

  } else {

    geno_mat <- geno_mat

    ntrait <- 0

    res <- NULL

  }

  # Data type statement
  data_type <- paste('data type', type)

  filecon <- file(paste0(out_name, ext))

  # Get the number of samples
  nnrow <- nrow(geno_mat)

  # Get the number of markers
  nmarker <- ncol(geno_mat)

  numbs <- paste(nnrow, nmarker, ntrait)

  # Recode genotype calls
  if (type == "f2 backcross") {

    geno_mat[geno_mat == homo1 | geno_mat == homo2] <- "A"
    geno_mat[geno_mat == het] <- "H"

  } else {

    geno_mat[geno_mat == homo1] <- "A"
    geno_mat[geno_mat == het] <- "H"
    geno_mat[geno_mat == homo2] <- "B"

  }

  # Get marker names and append asterisk to each name
  mks <- paste0('*', colnames(geno_mat))

  loc <- c() # Empty vector to store marker genotypes

  # Convert each marker data to MapMaker format
  for (n in seq_len(ncol(geno_mat))) {

    marker <- paste0(geno_mat[, n], collapse = "")

    loc[n] <- paste0(mks[n], " ", marker)

  }

  loc <- c(loc, res) # Combine marker and trait data
  writeLines(c(data_type, numbs, "", loc), filecon)
  close(filecon)

}


#' Re-order markers after linkage mapping to match marker order in linkage map
#' @param geno_dat A data frame of original genotype file used for linkage
#' mapping.
#' @param linakge_map A data frame of linkage map output from the `onemap`
#' package.
#' @param trait_id A character vector indicating the column names for traits
#' if any.
#' @param ... Other valid arguments that can be passed to the `proc_marker()`
#' function.
reorder_markers <- function (geno_dat,
                             linkage_map,
                             trait_id = NULL,
                             ...
) {

  ordered_markers <- linkage_map$V2

  # Subset only marker data without traits
  geno_mat <- geno_dat[, colnames(geno_dat) %in% ordered_markers]

  nmarkers <- ncol(geno_mat)
  n.ind <- nrow(geno_mat)

  res <- as.data.frame(matrix(NA, nrow = n.ind, ncol = nmarkers))
  colnames(res) <- ordered_markers
  #
  for (i in seq_len(nmarkers)) {

    res[, i] <- geno_mat[, ordered_markers[i]]

  }

  # Process traits if present
  if (!is.null(trait_id)) {

    # Subset only traits data without marker data
    trait <- geno_dat[, colnames(geno_dat) %in% trait_id]


  } else {

    trait <- NULL

  }

  # Column bind geno and trait data
  geno_mat <- cbind(res, trait)

  # Convert data to Mapmaker format and save in working directory
  do.call(proc_mapmaker, args = list(x = geno_mat,
                                     trait_id =  trait_id,
                                     ...))
  return(geno_mat)

}
