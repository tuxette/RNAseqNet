## User guide launcher
#' @title View RNAseqNet User's Guide
#' @name RNAseqNetUsersGuide
#' @export
#'
#' @importFrom utils browseURL
#'
#' @description Find the location of the RNAseqNet User's Guide and optionnaly
#' opens it
#' @param html logical. Should the document returned by the function be the
#' compiled PDF or the Rmd source. Default to \code{TRUE}
#' @param view logical. Should the document be opened using the default HTML
#' viewer? Default to \code{html}. It has no effect if \code{html = FALSE}
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @details The function \code{vignette("RNAseqNet")} will find the short
#' RNAseqNet vignette that describes how to obtain the RNAseqNet User's Guide.
#' The User's Guide is not itself a true vignette because it is not
#' automatically generated during the package build process. However, the
#' location of the Rmarkdown source is returned by the function if
#' \code{html = FALSE}.
#' If the operating system is not Windows, then the HTML viewer used is that
#' given by \code{Sys.getenv("R_BROWSER")}. The HTML viewer can be changed using
#'  \code{Sys.setenv(R_BROWSER = )}.
#'
#' @return Character string giving the file location. If \code{html = TRUE} and
#' \code{view = TRUE}, the HTML document reader is started and the User's Guide
#' is opened in it.
#'
#' @examples
#' RNAseqNetUsersGuide(view = FALSE)
#' RNAseqNetUsersGuide(html = FALSE)
#' \dontrun{RNAseqNetUsersGuide()}

RNAseqNetUsersGuide <- function(html = TRUE, view = html) {
  if (html) {
    f <- system.file("doc", "RNAseqNetUsersGuide.html", package = "RNAseqNet")
    if (view) {
      if (.Platform$OS.type == "windows")
        shell.exec(f)
      else browseURL(paste0("file://", f))
    }
  } else {
    f <- system.file("doc", "RNAseqNetUsersGuide.Rmd", package = "RNAseqNet")
    if (view) {
      warning("'RNAseqNetUserGuide.Rmd' can not be viewed.
              However, the location of the file is returned by the function.",
              call. = FALSE)
    }
  }
  return(f)
}
