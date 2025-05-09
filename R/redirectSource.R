#' redirectSource
#'
#' redirectSource will call a source specific redirect function if it exists
#' (called e.g. redirectTau), in which case the arguments are passed on to that
#' function. If such a function is not available \code{\link{redirect}} is called.
#' @param ... Additional arguments, passed on to the specific redirect function.
#' @inheritParams redirect
#' @return The result of the specific redirect function or \code{\link{redirect}}.
#' @author Pascal Sauer
#' @examples \dontrun{
#' f <- function() {
#'   redirectSource("Tau", target = "~/TauExperiment")
#'   # the following call will change directory
#'   # into ~/TauExperiment instead of <getConfig("sourcefolder")>/Tau
#'   readSource("Tau")
#' }
#' f()
#' # Tau is only redirected in the local environment of f,
#' # so it will use the usual source folder here
#' readSource("Tau")
#' }
#' @export
redirectSource <- function(type, target, ..., linkOthers = TRUE, local = TRUE) {
  if (isTRUE(local)) {
    local <- parent.frame()
  }

  specificRedirect <- get0(paste0("redirect", type), mode = "function")
  if (is.null(specificRedirect)) {
    if (...length() > 0) {
      warning("redirectSource calls madrat::redirect, so additional arguments are ignored.")
    }
    return(redirect(type = type, target = target, linkOthers = linkOthers, local = local))
  } else {
    return(specificRedirect(target = target, ..., linkOthers = linkOthers, local = local))
  }
}
