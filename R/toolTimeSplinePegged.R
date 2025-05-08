#' Smooth a magclass time series with a spline
#'
#' toolTimeSplinePegged fits a smoothing spline to each series in a magclass object.
#' Optionally, you can fix (“peg”) specified years exactly at their original values.
#'
#' @param x   A magclass object.
#' @param dof Degrees-of-freedom per 100 years (higher → more degrees of freedom, less smoothing; default 5).
#' @param pegged_years Integer vector of years (e.g. `c(2020, 2050, 2100)`) to hold fixed; NULL for none.
#' @param anchor_factor Numeric multiplier for anchor weights (default 10); larger values more strongly enforce pegging.
#'
#' @return A magclass object of the same shape, with each time series spline-smoothed.
#' @author Kristine Karstens, Felicitas Beier, Michael Crawford
#' @importFrom stats smooth.spline predict
#' @export

toolTimeSplinePegged <- function(x,
                                 dof = NULL,
                                 pegged_years = NULL, # NULL is no anchoring (old behavior)
                                 anchor_factor = 10) {
  ## 1) Input checks
  if (!is.magpie(x)) {
    stop("Input must be a magclass (MAgPIE) object!")
  }

  ## 2) Time axis & df calculation
  years <- getYears(x, as.integer = TRUE)
  nyr <- length(years)
  if (nyr < 2) {
    warning("Less than two time steps: nothing to smooth.")
    return(x)
  }
  timespan <- years[nyr] - years[1]

  if (is.null(dof)) {
    dof_par <- 5
  } else if (!is.numeric(dof) || dof < 1) {
    warning("Invalid dof; resetting to 5.")
    dof_par <- 5
  } else {
    dof_par <- dof
  }
  if (dof_par > 30) {
    warning("High dof vs. timespan may reduce smoothing effect.")
  }
  df_val <- timespan * dof_par / 100

  ## 3) Build weight vector
  if (is.null(pegged_years)) {
    # old behavior: no anchors
    wts <- rep(1, nyr)
    pegged_years_all <- NULL
  } else {
    # parse user‐supplied anchors (allow "yYYYY" or numeric)
    yrs_user <- as.integer(sub("^y", "", as.character(pegged_years), ignore.case = TRUE))
    if (!requireNamespace("magpiesets", quietly = TRUE)) {
      warning("magpiesets not available; historical anchoring skipped.")
      hist_yrs <- integer(0)
    } else {
      past_str <- magpiesets::findset("past") # e.g. "y1965",…
      hist_yrs <- as.integer(sub("^y", "", past_str, ignore.case = TRUE))
    }
    # combine user anchors + historical period + endpoints
    pegged_years_all <- unique(c(years[1], years[nyr], hist_yrs, yrs_user))
    # keep only years present in data
    pegged_years_all <- intersect(pegged_years_all, years)
    if (!all(yrs_user %in% years)) {
      stop("One or more user-supplied anchors not in data years.")
    }

    wts <- rep(1, nyr)
    wts[years %in% pegged_years_all] <- nyr * anchor_factor
  }

  ## 4) Preserve non-negativity flag
  negative_flag <- any(as.array(x) < 0)

  ## 5) Per-series spline (uses fit$y so no predict() call)
  tmpspline <- function(ts, df) {
    fit <- stats::smooth.spline(
      x            = years,
      y            = ts,
      w            = wts,
      df           = df,
      control.spar = list(high = 2)
    )
    fit$y
  }

  ## 6) Apply over time-series (dim 2 inner) – same as original
  arr_in <- as.array(x)
  arr_out <- apply(arr_in, c(1, 3), tmpspline, df = df_val)

  ## 7) Reconstruct magpie object
  dimnames(arr_out)[[1]] <- as.character(years)
  names(dimnames(arr_out))[1] <- getSets(x, fulldim = FALSE)[2]
  out <- as.magpie(arr_out, spatial = 2, temporal = 1)

  ## 8) Enforce non-negativity
  if (!negative_flag) out[out < 0] <- 0

  ## 9) Comment and return
  anchor_text <- if (is.null(pegged_years_all)) {
    "none"
  } else {
    paste(pegged_years_all, collapse = ",")
  }
  comment <- paste0(
    getComment(x),
    "; toolTimeSpline smoothed (anchors: ", anchor_text,
    "; df=", round(df_val, 2), ") [", date(), "]"
  )
  getComment(out) <- comment

  return(out)
}
