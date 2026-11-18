#' Launch the TaxonSampling Shiny application
#' @param ... Additional arguments passed to the application.
#' @export
launch_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Please install the 'shiny' package.")
  }

  app_dir <- system.file("shiny", package = "TaxonSampling")

  if (!nzchar(app_dir)) {
    stop("TaxonSampling Shiny application not found.")
  }

  shiny::runApp(app_dir, ...)
}
