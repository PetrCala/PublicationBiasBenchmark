#' @title MAIVE: Meta-Analysis Instrumental Variable Estimator
#'
#' @author Petr Cala \email{cala.p@@seznam.cz}
#'
#' @description
#' Implements the MAIVE method for publication bias correction using
#' instrumental variable estimation with variance instrumentation. MAIVE
#' addresses spurious precision in meta-analysis by instrumenting standard
#' errors with inverse sample sizes, providing consistent estimates even
#' when precision is manipulated through p-hacking.
#'
#' The method implements several estimators:
#' \itemize{
#'   \item PET (Precision-Effect Test): Linear precision-effect model
#'   \item PEESE (Precision-Effect Estimate with Standard Error): Quadratic model
#'   \item PET-PEESE: Conditional selection based on PET significance
#'   \item EK (Endogenous Kink): Flexible bias function with kink point
#'   \item WAIVE: Robust variant with outlier downweighting
#' }
#'
#' @param method_name Method identifier (automatically passed by framework)
#' @param data Data frame with yi (effect sizes), sei (standard errors),
#'   ni (sample sizes), and optionally study_id for clustering
#' @param settings List of method settings from method_settings.MAIVE()
#'
#' @return Single-row data frame with standardized output columns:
#'   \describe{
#'     \item{method}{Method identifier}
#'     \item{estimate}{Meta-analytic effect size estimate}
#'     \item{standard_error}{Standard error of estimate}
#'     \item{ci_lower, ci_upper}{95% confidence interval bounds (Anderson-Rubin if available)}
#'     \item{p_value}{Two-tailed p-value}
#'     \item{BF}{Bayes factor (NA for MAIVE)}
#'     \item{convergence}{Logical convergence indicator}
#'     \item{note}{Error messages if any}
#'     \item{first_stage_f}{First-stage F-statistic for instrument strength}
#'     \item{hausman_stat}{Hausman test statistic comparing IV vs OLS}
#'     \item{bias_p_value}{P-value for publication bias test}
#'     \item{used_ar_ci}{Whether Anderson-Rubin CI was used}
#'     \item{ar_ci_available}{Whether AR CI was computed successfully}
#'   }
#'
#' @details
#' MAIVE uses inverse sample sizes (1/N) as instruments for variances (SE^2)
#' in the first stage, then uses the instrumented variances in second-stage
#' PET/PEESE models. This approach provides consistent estimation when
#' precision is endogenous due to p-hacking or selective reporting.
#'
#' The Anderson-Rubin confidence interval is robust to weak instruments and
#' is automatically computed for unweighted IV estimators when feasible
#' (n < 5000). For weighted estimators or large samples, standard CIs are used.
#'
#' WAIVE extends MAIVE by downweighting: (1) negative residuals (spurious
#' precision) using exponential decay, and (2) extreme residuals (|z| > 2)
#' as potential outliers. This provides additional robustness against
#' publication bias and outliers.
#'
#' Available settings (see method_settings.MAIVE()):
#' \describe{
#'   \item{default}{PET-PEESE with IV, unweighted, levels first-stage}
#'   \item{PET}{PET with IV, unweighted}
#'   \item{PEESE}{PEESE with IV, unweighted}
#'   \item{EK}{Endogenous Kink model with IV}
#'   \item{weighted}{PET-PEESE with MAIVE-adjusted weighting}
#'   \item{WAIVE}{Robust WAIVE variant with outlier downweighting}
#'   \item{log_first_stage}{PET-PEESE with log-linear first stage}
#'   \item{no_IV}{Standard PET-PEESE without instrumentation (baseline)}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [run_method()], [method_settings()], [method_extra_columns()]
#'
#' @examples
#' \dontrun{
#' # Generate test data
#' data <- simulate_dgm("Stanley2017", condition_id = 1)
#'
#' # Apply default MAIVE (PET-PEESE with IV)
#' result <- run_method("MAIVE", data, "default")
#'
#' # Apply WAIVE variant
#' result_waive <- run_method("MAIVE", data, "WAIVE")
#'
#' # Apply weighted MAIVE
#' result_weighted <- run_method("MAIVE", data, "weighted")
#'
#' # View available configurations
#' method_settings("MAIVE")
#' }
#'
#' @export
method.MAIVE <- function(method_name, data, settings) {
  # ============================================================================
  # SECTION 1: INPUT EXTRACTION AND VALIDATION
  # ============================================================================

  # Extract required columns from PublicationBiasBenchmark format
  yi <- data$yi
  sei <- data$sei
  ni <- data$ni

  # Optional: Extract study IDs for clustering
  study_id <- data[["study_id"]]

  # Input validation

  # Check for required columns first (before checking content)
  if (is.null(yi) || length(yi) == 0) {
    stop("Effect sizes (yi) are required for MAIVE. ",
      "The data must include a 'yi' column with effect size estimates.",
      call. = FALSE
    )
  }

  if (is.null(sei) || length(sei) == 0) {
    stop("Standard errors (sei) are required for MAIVE. ",
      "The data must include a 'sei' column with standard errors.",
      call. = FALSE
    )
  }

  if (is.null(ni) || length(ni) == 0) {
    stop("Sample sizes (ni) are required for MAIVE. ",
      "The data must include a 'ni' column with positive sample sizes. ",
      "MAIVE uses inverse sample sizes (1/N) as instruments for variance.",
      call. = FALSE
    )
  }

  if (length(yi) < 3) {
    stop("MAIVE requires at least 3 effect size estimates for reliable estimation",
      call. = FALSE
    )
  }

  if (any(is.na(yi))) {
    stop("Effect sizes (yi) contain missing values", call. = FALSE)
  }

  if (any(is.na(sei))) {
    stop("Standard errors (sei) contain missing values", call. = FALSE)
  }

  if (any(sei <= 0)) {
    stop("All standard errors must be positive (sei > 0)", call. = FALSE)
  }

  if (any(is.na(ni))) {
    stop("Sample sizes (ni) contain missing values", call. = FALSE)
  }

  if (any(ni <= 0)) {
    stop("All sample sizes must be positive (ni > 0)", call. = FALSE)
  }


  # ============================================================================
  # SECTION 2: DATA PREPARATION (Format conversion)
  # ============================================================================

  # Convert from PublicationBiasBenchmark format to MAIVE format
  # PublicationBiasBenchmark: yi, sei, ni, study_id
  # MAIVE expects:            bs, sebs, Ns, study_id

  maive_data <- data.frame(
    bs = yi, # Effect sizes
    sebs = sei, # Standard errors
    Ns = ni # Sample sizes
  )

  # Add study_id if provided (needed for clustering)
  if (!is.null(study_id)) {
    maive_data$study_id <- study_id
  }


  # ============================================================================
  # SECTION 3: EXTRACT AND VALIDATE SETTINGS
  # ============================================================================

  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x

  # Extract settings with defaults
  use_waive <- settings$use_waive %||% FALSE
  maive_method <- settings$method %||% 3
  maive_weight <- settings$weight %||% 0
  maive_instrument <- settings$instrument %||% 1
  maive_studylevel <- settings$studylevel %||% 0
  maive_SE <- settings$SE %||% 0
  maive_AR <- settings$AR %||% 1
  maive_first_stage <- settings$first_stage %||% 0

  # Collect any automatic adjustments for the note field
  adjustment_notes <- character(0)

  # --- studylevel auto-selection / downgrades --------------------------------
  # If study_id is missing, clustering/fixed effects are impossible.
  if (is.null(study_id) && maive_studylevel > 0) {
    maive_studylevel <- 0
    adjustment_notes <- c(adjustment_notes, "studylevel auto-adjusted to 0 (no study_id)")
  }

  # If studylevel is left at 0 but we have repeated study_id values, clustering
  # is typically sensible (e.g., PRE DGMs with multiple estimates per study).
  if (!is.null(study_id) && maive_studylevel == 0) {
    n_clusters <- length(unique(study_id))
    if (n_clusters < length(study_id)) {
      maive_studylevel <- 2
      adjustment_notes <- c(adjustment_notes, "studylevel auto-adjusted to 2 (cluster) due to repeated study_id")
    }
  }

  # If clustering is requested but there is only one cluster, clubSandwich will
  # fail (and clustering is not meaningful anyway). Drop cluster component.
  if (!is.null(study_id) && maive_studylevel %/% 2 == 1) {
    n_clusters <- length(unique(study_id))
    if (n_clusters < 2) {
      dummy_component <- maive_studylevel %% 2
      maive_studylevel <- dummy_component
      adjustment_notes <- c(adjustment_notes, "cluster component dropped (only 1 cluster in study_id)")
    }
  }

  # --- instrumentation feasibility checks ------------------------------------
  # MAIVE instruments SE^2 with 1/N. If N does not vary, the instrument has no
  # variation and IV is not identified; older MAIVE versions can crash in this
  # situation (aliased slopes -> vcov indexing errors).
  if (maive_instrument == 1) {
    n_unique_n <- length(unique(ni))
    if (n_unique_n < 2) {
      maive_instrument <- 0
      maive_AR <- 0
      adjustment_notes <- c(adjustment_notes, "instrument auto-adjusted to 0 (ni has no variation); AR disabled")
    }
  }

  # AR intervals are expensive and sometimes unavailable; disable automatically
  # for large samples to avoid memory/time issues.
  if (maive_AR == 1 && length(yi) > 5000) {
    maive_AR <- 0
    adjustment_notes <- c(adjustment_notes, "AR disabled automatically for n > 5000")
  }

  # Validate settings ranges
  if (!maive_method %in% 1:4) {
    stop("Invalid method parameter: must be 1 (PET), 2 (PEESE), 3 (PET-PEESE), or 4 (EK)",
      call. = FALSE
    )
  }

  if (!maive_weight %in% 0:2) {
    stop("Invalid weight parameter: must be 0 (none), 1 (inverse-variance), or 2 (MAIVE-adjusted)",
      call. = FALSE
    )
  }

  if (!maive_instrument %in% 0:1) {
    stop("Invalid instrument parameter: must be 0 (no IV) or 1 (use IV)",
      call. = FALSE
    )
  }


  # ============================================================================
  # SECTION 4: CALL MAIVE OR WAIVE
  # ============================================================================

  # Prepare detailed error context for debugging
  error_context <- list(
    n_studies = length(yi),
    has_study_id = !is.null(study_id),
    method = maive_method,
    weight = maive_weight,
    instrument = maive_instrument,
    studylevel = maive_studylevel,
    SE = maive_SE,
    AR = maive_AR,
    first_stage = maive_first_stage,
    use_waive = use_waive
  )

  # Call the appropriate MAIVE function
  maive_result <- tryCatch(
    {
      if (use_waive) {
        # Call WAIVE (weighted adjusted IV estimator)
        MAIVE::waive(
          dat = maive_data,
          method = maive_method,
          weight = maive_weight,
          instrument = maive_instrument,
          studylevel = maive_studylevel,
          SE = maive_SE,
          AR = maive_AR,
          first_stage = maive_first_stage
        )
      } else {
        # Call standard MAIVE
        MAIVE::maive(
          dat = maive_data,
          method = maive_method,
          weight = maive_weight,
          instrument = maive_instrument,
          studylevel = maive_studylevel,
          SE = maive_SE,
          AR = maive_AR,
          first_stage = maive_first_stage
        )
      }
    },
    error = function(e) {
      # Return a standardized non-converged row instead of stopping; this keeps
      # run_method(..., silent=FALSE) from printing noisy errors while still
      # recording what happened.
      error_msg <- paste0(
        "MAIVE execution failed: ", conditionMessage(e),
        "\nContext: n=", error_context$n_studies,
        ", method=", error_context$method,
        ", weight=", error_context$weight,
        ", instrument=", error_context$instrument,
        ", studylevel=", error_context$studylevel,
        ", SE=", error_context$SE,
        ", AR=", error_context$AR,
        ", first_stage=", error_context$first_stage,
        ", use_waive=", error_context$use_waive
      )

      if (length(adjustment_notes) > 0) {
        error_msg <- paste0("Auto-adjustments: ", paste(adjustment_notes, collapse = "; "), "\n", error_msg)
      }

      return(
        create_empty_result(
          method_name = method_name,
          note = error_msg,
          extra_columns = method_extra_columns.MAIVE(method_name)
        )
      )
    }
  )

  # If the tryCatch returned an empty-result data frame, pass it through.
  if (is.data.frame(maive_result) && isTRUE(maive_result$convergence[1] == FALSE)) {
    return(maive_result)
  }


  # ============================================================================
  # SECTION 5: EXTRACT RESULTS FROM MAIVE OUTPUT
  # ============================================================================

  # Extract main estimate and standard error
  estimate <- maive_result$beta
  se <- maive_result$SE

  # Validate extracted values
  if (is.null(estimate) || is.na(estimate)) {
    stop("MAIVE returned NULL or NA estimate", call. = FALSE)
  }

  if (is.null(se) || is.na(se) || se <= 0) {
    stop("MAIVE returned invalid standard error (NULL, NA, or non-positive)", call. = FALSE)
  }


  # ============================================================================
  # SECTION 6: EXTRACT CONFIDENCE INTERVALS
  # ============================================================================

  # Anderson-Rubin CI is preferred when available (robust to weak instruments)
  # MAIVE may return AR_CI as:
  #   - Numeric vector c(lower, upper)
  #   - String "NA"
  #   - Named vector c(lower=NA, upper=NA)
  #   - NULL

  ar_ci <- maive_result$AR_CI

  # Robust checking for valid AR CI
  is_valid_ar_ci <- !is.null(ar_ci) &&
    is.numeric(ar_ci) &&
    length(ar_ci) == 2 &&
    !anyNA(ar_ci) &&
    is.finite(ar_ci[1]) &&
    is.finite(ar_ci[2])

  if (is_valid_ar_ci) {
    # Use Anderson-Rubin CI
    ci_lower <- ar_ci[1]
    ci_upper <- ar_ci[2]
    used_ar_ci <- TRUE
    ar_ci_available <- TRUE
  } else {
    # Fallback to standard Wald CI
    ci_lower <- estimate - 1.96 * se
    ci_upper <- estimate + 1.96 * se
    used_ar_ci <- FALSE

    # Check if AR CI was attempted but failed
    ar_ci_available <- !is.null(ar_ci) &&
      (is.character(ar_ci) ||
        (is.numeric(ar_ci) && anyNA(ar_ci)))
  }


  # ============================================================================
  # SECTION 7: COMPUTE P-VALUE
  # ============================================================================

  # Two-tailed p-value using Wald test
  z_stat <- estimate / se
  p_value <- 2 * (1 - pnorm(abs(z_stat)))


  # ============================================================================
  # SECTION 8: EXTRACT EXTRA DIAGNOSTICS
  # ============================================================================

  # First-stage F-test for instrument strength
  # MAIVE returns "NA" string when instrument=0
  f_test <- maive_result[["F-test"]]
  if (!is.null(f_test) && is.character(f_test) && f_test == "NA") {
    f_test <- NA_real_
  } else if (!is.null(f_test)) {
    f_test <- as.numeric(f_test)
  } else {
    f_test <- NA_real_
  }

  # Hausman test statistic (IV vs OLS comparison)
  hausman_stat <- maive_result$Hausman
  if (is.null(hausman_stat)) {
    hausman_stat <- NA_real_
  }

  # Publication bias p-value
  pbias_pval <- maive_result$pbias_pval %||% maive_result[["pub bias p-value"]]
  if (is.null(pbias_pval)) {
    pbias_pval <- NA_real_
  }


  # ============================================================================
  # SECTION 9: PREPARE NOTE FIELD
  # ============================================================================

  # Combine any adjustment notes
  note_parts <- adjustment_notes

  if (used_ar_ci) {
    note_parts <- c(note_parts, "Using Anderson-Rubin CI")
  } else if (ar_ci_available) {
    note_parts <- c(note_parts, "AR CI attempted but unavailable; using Wald CI")
  }

  note <- if (length(note_parts) > 0) {
    paste(note_parts, collapse = "; ")
  } else {
    NA_character_
  }


  # ============================================================================
  # SECTION 10: RETURN STANDARDIZED OUTPUT
  # ============================================================================

  result <- data.frame(
    # Required columns (9 total)
    method = method_name,
    estimate = estimate,
    standard_error = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    BF = NA_real_, # MAIVE does not compute Bayes factors
    convergence = TRUE,
    note = note,

    # Extra columns (5 total - declared in method_extra_columns.MAIVE)
    first_stage_f = f_test,
    hausman_stat = hausman_stat,
    bias_p_value = pbias_pval,
    used_ar_ci = used_ar_ci,
    ar_ci_available = ar_ci_available,
    stringsAsFactors = FALSE
  )

  return(result)
}


#' @title MAIVE Method Settings
#'
#' @author Petr Cala \email{cala.p@@seznam.cz}
#'
#' @description
#' Defines available configurations for the MAIVE method. Each configuration
#' specifies different combinations of estimators, weighting schemes, and
#' methodological choices.
#'
#' @param method_name Method identifier (automatically passed)
#'
#' @return Named list of settings configurations. Each configuration is a list
#'   with the following parameters:
#'   \describe{
#'     \item{method}{Integer: 1=PET, 2=PEESE, 3=PET-PEESE, 4=EK}
#'     \item{weight}{Integer: 0=none, 1=inverse-variance, 2=MAIVE-adjusted}
#'     \item{instrument}{Integer: 0=no IV, 1=use IV}
#'     \item{studylevel}{Integer: 0=none, 1=fixed effects, 2=cluster, 3=both}
#'     \item{SE}{Integer: 0=CR0, 1=CR1, 2=CR2, 3=wild bootstrap}
#'     \item{AR}{Integer: 0=no, 1=yes - compute Anderson-Rubin CI}
#'     \item{first_stage}{Integer: 0=levels, 1=log - first-stage specification}
#'     \item{use_waive}{Logical: FALSE=maive(), TRUE=waive()}
#'   }
#'
#' @details
#' Available configurations:
#' \describe{
#'   \item{default}{PET-PEESE with IV, unweighted, cluster SE, levels first-stage}
#'   \item{PET}{Linear PET with IV, unweighted}
#'   \item{PEESE}{Quadratic PEESE with IV, unweighted}
#'   \item{EK}{Endogenous Kink model with IV}
#'   \item{weighted}{PET-PEESE with MAIVE-adjusted inverse-variance weighting}
#'   \item{WAIVE}{Robust variant with outlier downweighting}
#'   \item{log_first_stage}{PET-PEESE with log-linear first stage (for heteroskedasticity)}
#'   \item{no_IV}{Standard PET-PEESE without instrumentation (baseline comparison)}
#' }
#'
#' All configurations use SE=0 (CR0/Huber-White) for computational speed
#' in benchmarking. Settings are immutable once published for reproducibility.
#'
#' @export
method_settings.MAIVE <- function(method_name) {
  settings <- list(

    # -------------------------------------------------------------------------
    # GROUP A: CORE MAIVE METHODS (IV with different estimators)
    # -------------------------------------------------------------------------
    "default" = list(
      method = 3, # PET-PEESE conditional selection
      weight = 0, # No weighting (unweighted)
      instrument = 1, # Use variance instrumentation
      studylevel = 2, # Cluster-robust standard errors
      SE = 0, # CR0 (Huber-White) - fastest
      AR = 1, # Compute Anderson-Rubin CI
      first_stage = 0, # Levels specification (linear)
      use_waive = FALSE # Use maive() not waive()
    ),
    "PET" = list(
      method = 1, # Linear PET only
      weight = 0,
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 1,
      first_stage = 0,
      use_waive = FALSE
    ),
    "PEESE" = list(
      method = 2, # Quadratic PEESE only
      weight = 0,
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 1,
      first_stage = 0,
      use_waive = FALSE
    ),
    "EK" = list(
      method = 4, # Endogenous Kink model
      weight = 0,
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 0, # AR not available for EK
      first_stage = 0,
      use_waive = FALSE
    ),

    # -------------------------------------------------------------------------
    # GROUP B: WEIGHTED MAIVE VARIANTS
    # -------------------------------------------------------------------------

    "weighted" = list(
      method = 3, # PET-PEESE
      weight = 2, # MAIVE-adjusted inverse-variance weighting
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 1, # AR available for weight=2
      first_stage = 0,
      use_waive = FALSE
    ),
    "WAIVE" = list(
      method = 3, # PET-PEESE
      weight = 0,
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 0, # Typically don't use AR with WAIVE
      first_stage = 0,
      use_waive = TRUE # Use waive() function
    ),

    # -------------------------------------------------------------------------
    # GROUP C: ALTERNATIVE SPECIFICATIONS
    # -------------------------------------------------------------------------

    "log_first_stage" = list(
      method = 3, # PET-PEESE
      weight = 0,
      instrument = 1,
      studylevel = 2,
      SE = 0,
      AR = 1,
      first_stage = 1, # Log-linear first stage with smearing
      use_waive = FALSE
    ),
    "no_IV" = list(
      method = 3, # PET-PEESE
      weight = 0,
      instrument = 0, # No instrumentation (standard estimator)
      studylevel = 2,
      SE = 0,
      AR = 0, # AR only for IV estimators
      first_stage = 0,
      use_waive = FALSE
    )
  )

  return(settings)
}


#' @title MAIVE Extra Output Columns
#'
#' @author Petr Cala \email{cala.p@@seznam.cz}
#'
#' @description
#' Declares additional output columns beyond the standard 9 required columns
#' that are returned by the MAIVE method.
#'
#' @param method_name Method identifier (automatically passed)
#'
#' @return Character vector of extra column names
#'
#' @details
#' The MAIVE method returns 5 additional diagnostic columns:
#' \describe{
#'   \item{first_stage_f}{First-stage F-statistic testing instrument strength.
#'     High values (>10) indicate strong instruments. NA when instrument=0.}
#'   \item{hausman_stat}{Hausman test statistic comparing IV vs OLS estimates.
#'     Large values suggest endogeneity (need for IV). Based on difference-in-estimators.}
#'   \item{bias_p_value}{P-value for publication bias test using instrumented
#'     precision-effect relationship (FAT test). Low values indicate bias.}
#'   \item{used_ar_ci}{Logical indicator of whether Anderson-Rubin confidence
#'     interval was used (TRUE) or standard Wald CI (FALSE).}
#'   \item{ar_ci_available}{Logical indicator of whether AR CI computation
#'     was attempted and potentially available.}
#' }
#'
#' @export
method_extra_columns.MAIVE <- function(method_name) {
  c(
    "first_stage_f", # First-stage F-test for instrument strength
    "hausman_stat", # Hausman test: IV vs OLS comparison
    "bias_p_value", # Publication bias test p-value
    "used_ar_ci", # Logical: whether Anderson-Rubin CI was used
    "ar_ci_available" # Logical: whether AR CI was computed/available
  )
}
