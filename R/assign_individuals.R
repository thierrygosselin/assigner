#' Assign individuals to reference strata
#'
#' @description
#' Estimate allele frequencies in reference strata and calculate the genotype
#' log-likelihood of every individual under every candidate stratum. The native
#' likelihood engine is implemented entirely in R and returns the complete
#' individual-by-stratum score table required for assignment diagnostics and
#' Dlr calculations.
#'
#' @param data A GDS file or object, or tidy genomic data accepted by
#'   [genometranslator::read_genome()]. Genotypes are converted internally to
#'   diploid alternate-allele dosage (`0`, `1`, `2`, or `NA`).
#' @param strata Optional strata file or object accepted by
#'   [genometranslator::read_strata()]. When `NULL`, `STRATA` must be present in
#'   `data`. Default: \code{strata = NULL}.
#' @param engine Assignment engine: native population-genetic likelihood,
#'   elastic-net multinomial regression, XGBoost, or experimental TabPFN.
#'   Default: \code{engine = "likelihood"}.
#' @param genotype.method Genotype information used by the native likelihood
#'   engine: called genotypes (`"GT"`), log10-scaled or normalized genotype
#'   likelihoods (`"GL"`), or Phred-scaled genotype likelihoods (`"PL"`).
#'   GL requires `GL_HOM_REF`, `GL_HET`, and `GL_HOM_ALT`; PL uses the analogous
#'   `PL_` columns. Default: \code{genotype.method = "GT"}.
#' @param leave.one.out Logical. When assigning a reference individual, exclude
#'   that individual from allele-frequency estimation in its home stratum.
#'   Default: \code{leave.one.out = TRUE}.
#' @param frequency.floor Minimum allele frequency used to prevent zero
#'   likelihoods. With `NULL`, the floor is estimated independently for every
#'   stratum and marker as `1 / (2 * N + 1)`, where `N` is the number of
#'   genotyped reference individuals. Default:
#'   \code{frequency.floor = NULL}.
#' @param folds Number of stratified outer validation folds used by machine-
#'   learning engines. Default: \code{folds = 5}.
#' @param random.seed Integer controlling fold creation and stochastic engines.
#'   Default: \code{random.seed = NULL}.
#' @param effective.n.tolerance Relative difference in effective reference
#'   sample size above which GL/PL analyses issue an imbalance warning.
#'   Default: \code{effective.n.tolerance = 0.20}.
#' @param em.tolerance Convergence tolerance for GL/PL allele-frequency
#'   estimation. Default: \code{em.tolerance = 1e-7}.
#' @param em.max.iterations Maximum EM iterations for each stratum-marker
#'   allele-frequency estimate. Default: \code{em.max.iterations = 200}.
#' @param unsampled.source Logical. For GL/PL data, calculate a coverage-aware,
#'   reference-standardized score for detecting individuals that may originate
#'   from an unsampled source. Requires `ALLELE_REF_DEPTH` and
#'   `ALLELE_ALT_DEPTH`. Default: \code{unsampled.source = FALSE}.
#' @param unsampled.cutoff Lower standardized-score threshold used to flag a
#'   possible unsampled source. This is an exploratory value that must be
#'   calibrated for the dataset. Default: \code{unsampled.cutoff = -3}.
#' @param sequencing.error Per-read error probability used by the
#'   unsampled-source expectation model. Default:
#'   \code{sequencing.error = 0.01}.
#' @param evaluation.data Optional GL/PL dataset used only to evaluate focal
#'   individuals while reference allele frequencies remain estimated from
#'   `data`. It must contain the same individual-marker combinations. This
#'   supports coverage-matched or deliberately downsampled leave-one-out
#'   evaluation. Default: \code{evaluation.data = NULL}.
#' @param marker.blocks Optional marker-block definition used to summarize how
#'   assignment likelihood accumulates across the genome. Supply a named vector,
#'   a table containing `MARKERS` and `MARKER_BLOCK`, or an integer number of
#'   interleaved blocks. Default: \code{marker.blocks = NULL}.
#' @param reporting.unit Optional metadata column defining broader reporting
#'   units containing one or more `STRATA`. A second complete assignment result
#'   is returned without replacing collection-level results. Default:
#'   \code{reporting.unit = NULL}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{assignment}{One row per individual with its current stratum,
#'     inferred stratum, home likelihood, maximum likelihood, likelihood-ratio
#'     statistic, and second-best stratum.}
#'   \item{likelihoods}{The complete table containing the log-likelihood of
#'     every individual under every candidate stratum.}
#'   \item{allele.frequencies}{Reference-stratum allele-frequency estimates
#'     before individual leave-one-out adjustment.}
#'   \item{effective.sample.size}{For GL/PL, mean Fisher-information effective
#'     sample size by reference stratum, calculated across loci with estimated
#'     allele frequency between 0.05 and 0.95.}
#'   \item{effective.sample.size.markers}{For GL/PL, effective sample size for
#'     every reference-stratum and marker combination.}
#'   \item{individual.effective.sample.size}{For GL/PL, the mean effective
#'     information contributed by each reference individual across markers.}
#'   \item{marker.block.likelihoods}{When requested, individual-by-candidate
#'     likelihoods calculated separately for every marker block.}
#'   \item{reporting.unit}{When requested, a nested assignment result for the
#'     broader reporting-unit level.}
#'   \item{unsampled.calibration}{When requested, reference-stratum mean and
#'     standard deviation used to standardize the unsampled-source score.}
#' }
#'
#' @details
#' For genotype dosage `0`, `1`, and `2`, the Hardy-Weinberg genotype
#' probabilities are `(1-p)^2`, `2p(1-p)`, and `p^2`, respectively, where `p`
#' is the alternate-allele frequency in a candidate stratum. Probabilities are
#' accumulated on the log scale. Missing genotypes do not contribute to the
#' score. For called genotypes, `likelihoods$Z_LIKELIHOOD` standardizes each
#' candidate likelihood using the expected mean and variance from exactly the
#' markers observed in that individual. This makes individuals with different
#' missing-data patterns more comparable without imputing genotypes.
#'
#' The complete `likelihoods` table is intentionally retained. The best and
#' second-best assignments alone are insufficient for Dlr calculations and can
#' conceal assignment ambiguity among additional candidate strata.
#'
#' This is a classification method: candidate source strata must be defined and
#' represented by reference samples. The native engine assigns an individual to
#' the stratum with the greatest multilocus genotype likelihood. It does not
#' perform a Monte Carlo exclusion test, estimate migration rates, or establish
#' that the true source population was sampled. Those are separate inferential
#' questions (Manel et al. 2005; Paetkau et al. 2004).
#'
#' Reference individuals should normally be evaluated with
#' `leave.one.out = TRUE`, because using the focal genotype to estimate its home
#' allele frequencies biases assignment in favour of that stratum. Alleles that
#' are observed in an individual but absent from a finite reference sample would
#' otherwise produce a zero likelihood. `frequency.floor` controls the
#' correction. The adaptive default reflects the number of called reference
#' individuals at each marker; set `frequency.floor = 0.005` to use the fixed
#' baseline described by Paetkau et al. (2004) and used by GenoDive.
#'
#' Assignment reliability depends on representative reference sampling,
#' differentiation among candidate strata, marker information, and reference
#' sample size. Very small reference strata can yield unstable allele-frequency
#' estimates. A high maximum likelihood is comparative evidence among the
#' supplied candidates, not proof that the assigned stratum is the true source.
#'
#' @section GT, GL, and PL assignment features:
#' \itemize{
#'   \item Choose the genotype representation explicitly with
#'     `genotype.method = "GT"`, `"GL"`, or `"PL"`; the function does not guess.
#'   \item GL/PL reference allele frequencies are estimated by expectation-
#'     maximization rather than from hard genotype calls.
#'   \item GL/PL assignment integrates over `0/0`, `0/1`, and `1/1`, retaining
#'     genotype uncertainty for the focal individual.
#'   \item Leave-one-out recalculates the focal reference individual's home-
#'     stratum allele frequencies.
#'   \item Observed Fisher information is summarized as effective reference
#'     sample size, revealing imbalances caused by sample number, coverage, and
#'     genotype uncertainty.
#'   \item With the adaptive default, GL/PL frequency floors use marker-specific
#'     effective reference size when available, falling back to the number of
#'     usable reference individuals.
#'   \item `evaluation.data` can hold lower-coverage or deliberately downsampled
#'     GL/PL values for honest coverage-matched evaluation while the original
#'     reference data continue to estimate allele frequencies.
#'   \item `unsampled.source = TRUE` adds a coverage-aware, reference-standardized
#'     diagnostic for individuals that may not originate from any sampled source.
#'   \item The complete individual-by-stratum likelihood table is retained for
#'     ambiguity assessment, Dlr, and downstream analyses.
#' }
#'
#' @section Relationship to gsi_sim and rubias:
#' Eric C. Anderson's
#' \href{https://github.com/eriqande/gsi_sim}{gsi_sim} supports genetic stock
#' identification simulations and is used by \code{\link{evaluate_assignment}}
#' to evaluate assignment performance across marker panels and resampling
#' designs. Eric C. Anderson also developed assigner's original `gsi_sim`
#' wrapper.
#'
#' \href{https://github.com/eriqande/rubias}{rubias}, developed by Eric C.
#' Anderson and Ben Moran, implements Bayesian inference for the conditional
#' genetic stock identification model. It is designed principally to estimate
#' reporting-unit or population proportions in a mixture while also producing
#' individual posterior assignments. In that joint analysis, estimated mixture
#' proportions contribute information to individual membership.
#'
#' By contrast, `assign_individuals()` evaluates each individual against the
#' supplied reference strata and retains its full source-likelihood table. It
#' does not estimate mixture proportions. Use rubias when mixed-stock composition
#' is the principal estimand; use `assign_individuals()` for direct individual
#' assignment, GL/PL uncertainty, Dlr, and unsampled-source exploration.
#'
#' For GL/PL, allele frequencies are estimated from genotype likelihoods by an
#' expectation-maximization algorithm. The genotype likelihood of the focal
#' sample is marginalized over `0/0`, `0/1`, and `1/1`, rather than replacing
#' uncertainty with a hard genotype call. Effective sample size follows the
#' observed-Fisher-information formulation of DeSaix et al. (2024).
#'
#' With `unsampled.source = TRUE`, the function calculates a coverage-aware raw
#' score by comparing the observed log probability with the expectation under a
#' binomial read-error model, then standardizes it using leave-one-out scores of
#' reference individuals from the same source. This is an exploratory adaptation
#' of DeSaix et al. (2024). The cutoff is dataset-specific: inspect reference
#' score distributions and known controls rather than treating `-3` as universal.
#' Use [compare_unsampled_scores()] to compare this independently implemented
#' diagnostic with scores generated by WGSassign. Agreement should be evaluated
#' on simulated data, known-origin controls, and realistic coverage profiles
#' before selecting a cutoff.
#'
#' @section Marker blocks and reporting units:
#' Linked loci do not provide independent evidence. `marker.blocks` retains the
#' full-data assignment while also returning likelihood contributions for
#' biologically defined linkage groups, chromosomes, or reproducible marker
#' partitions. These are sensitivity diagnostics, not independent replicates.
#'
#' `reporting.unit` evaluates a broader classification in addition to the
#' original `STRATA` analysis. Collections and reporting units answer different
#' questions; the reporting-unit result is nested under `result$reporting.unit`
#' so fine-scale ambiguity is not discarded.
#'
#' @section Missing data:
#' Missingness can encode sequencing batch, library, coverage, or processing
#' differences and can therefore create misleading assignment accuracy. The
#' likelihood engine ignores missing focal genotypes, estimates each
#' stratum-marker frequency from called reference genotypes only, and compares
#' candidate strata using the same called-marker set for an individual.
#' Elastic net uses training-fold mean imputation. XGBoost uses its native
#' missing-value routing. TabPFN uses its native preprocessing, including a
#' missingness indicator. All preprocessing and model fitting occur within the
#' training portion of each outer fold. Inspect the returned missingness table
#' and test whether missingness is associated with strata or processing batch.
#'
#' @section Data quality and assignment assumptions:
#' Assignment can produce a confident-looking result from technical structure
#' rather than population structure. Before interpreting the output, examine:
#' \itemize{
#'   \item \strong{Candidate sources and reference sampling.} The likelihood
#'     analysis assumes that candidate source populations are defined in advance,
#'     sampled representatively, and sufficiently differentiated. If the true
#'     source is absent, the function still reports the best of the supplied
#'     candidates. Unequal or small reference samples make allele-frequency
#'     estimates unevenly precise.
#'   \item \strong{Hardy-Weinberg assumptions.} The native diploid likelihood
#'     uses Hardy-Weinberg genotype probabilities within each stratum. Strong
#'     departures can reflect genotyping artefacts, inbreeding, family structure,
#'     admixture, or incorrectly pooled populations. Investigate the cause;
#'     Hardy-Weinberg departure is not by itself a reason to discard a marker.
#'   \item \strong{Linkage disequilibrium.} Multilocus likelihoods treat marker
#'     contributions as independent. With linked markers the result is a
#'     composite-likelihood approximation: assignment rankings can remain useful,
#'     but likelihood differences generally overstate independent information
#'     and should not be treated as calibrated uncertainty (DeSaix et al. 2024).
#'     Pruning is therefore not always required for classification, but sensitivity
#'     across pruned panels or biologically justified linkage blocks is important.
#'   \item \strong{Genotyping error.} Allelic dropout, paralogs, strand or allele
#'     coding errors, low-depth calls, and batch-specific calling differences can
#'     resemble population differentiation. Replicates, depth and balance
#'     diagnostics, reproducibility checks, and consistent joint genotyping are
#'     especially important.
#'   \item \strong{Missingness and batch effects.} Population, extraction batch,
#'     library, sequencing lane, coverage, and processing history must not be
#'     confounded. Compare missingness and other QC metrics across both strata and
#'     technical batches. Random cross-validation does not remove a confounded
#'     batch effect; when possible, validate on an independent batch.
#'   \item \strong{Relatedness, duplicates, and family structure.} Close relatives
#'     split across training and evaluation data can inflate apparent accuracy.
#'     Remove unintended duplicates and construct holdouts by family, collection,
#'     site, or batch when those groupings could leak information.
#'   \item \strong{Marker frequency and selection.} Very rare alleles and low
#'     minor-allele counts are unstable under leave-one-out estimation. Marker
#'     selection performed on the complete dataset leaks information into the
#'     evaluation; repeat selection inside each training fold. Anderson (2010)
#'     calls the resulting upward bias in predicted accuracy \emph{high-grading
#'     bias}. Ordinary leave-one-out assignment does not remove it because the
#'     focal individual has already influenced which loci were selected. When a
#'     panel is chosen for its apparent differentiation among the reference
#'     samples, reserve untouched holdout individuals or use a complete
#'     Training-Holdout-Leave-one-out design in which marker ranking occurs only
#'     in the training data.
#' }
#'
#' This function performs assignment, not general genomic filtering. Prepare
#' and quality-control the input with
#' \href{https://thierrygosselin.github.io/genometranslator/}{genometranslator}
#' and \href{https://thierrygosselin.github.io/radr/}{radr} first.
#'
#' @references Manel S, Gaggiotti OE, Waples RS (2005). Assignment methods:
#'   matching biological questions with appropriate techniques. Trends in
#'   Ecology & Evolution, 20(3), 136-142.
#'   \doi{10.1016/j.tree.2004.12.004}.
#' @references Paetkau D, Calvert W, Stirling I, Strobeck C (1995).
#'   Microsatellite analysis of population structure in Canadian polar bears.
#'   Molecular Ecology, 4, 347-354.
#' @references Paetkau D, Slade R, Burden M, Estoup A (2004). Genetic assignment
#'   methods for the direct, real-time estimation of migration rate: a
#'   simulation-based exploration of accuracy and power. Molecular Ecology, 13,
#'   55-65. \doi{10.1046/j.1365-294X.2003.02008.x}.
#' @references Cornuet J-M, Piry S, Luikart G, Estoup A, Solignac M (1999).
#'   New methods employing multilocus genotypes to select or exclude populations
#'   as origins of individuals. Genetics, 153, 1989-2000.
#' @references DeSaix MG, Rodriguez MD, Ruegg KC, Anderson EC (2024). Population
#'   assignment from genotype likelihoods for low-coverage whole-genome
#'   sequencing data. Methods in Ecology and Evolution, 15, 493-510.
#'   \doi{10.1111/2041-210X.14286}.
#' @references Moran BM, Anderson EC (2019). Bayesian inference from the
#'   conditional genetic stock identification model. Canadian Journal of
#'   Fisheries and Aquatic Sciences, 76(4), 551-560.
#'   \doi{10.1139/cjfas-2018-0016}.
#' @references Anderson EC (2010). Assessing the power of informative subsets
#'   of loci for population assignment: standard methods are upwardly biased.
#'   Molecular Ecology Resources, 10(4), 701-710.
#'   \doi{10.1111/j.1755-0998.2010.02846.x}.
#'
#' @examples
#' \dontrun{
#' result <- assigner::assign_individuals(
#'   data = genome,
#'   leave.one.out = TRUE
#' )
#'
#' result$assignment
#' result$likelihoods
#'
#' gl_result <- assigner::assign_individuals(
#'   data = low_coverage_genome,
#'   genotype.method = "GL",
#'   unsampled.source = TRUE
#' )
#' }
#'
#' @export
assign_individuals <- function(
  data,
  strata = NULL,
  engine = c("likelihood", "elastic_net", "xgboost", "tabpfn"),
  genotype.method = c("GT", "GL", "PL"),
  leave.one.out = TRUE,
  frequency.floor = NULL,
  folds = 5L,
  random.seed = NULL,
  effective.n.tolerance = 0.20,
  em.tolerance = 1e-7,
  em.max.iterations = 200L,
  unsampled.source = FALSE,
  unsampled.cutoff = -3,
  sequencing.error = 0.01,
  evaluation.data = NULL,
  marker.blocks = NULL,
  reporting.unit = NULL,
  verbose = TRUE
) {
  .start <- tgbase::startup(
    f.name = "assign_individuals", package = "assigner", verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (missing(data)) rlang::abort("`data` is required.")
  engine <- match.arg(engine)
  genotype.method <- match.arg(genotype.method)
  if (engine != "likelihood") {
    if (genotype.method != "GT") {
      rlang::abort("`genotype.method = \"GL\"` and `\"PL\"` require `engine = \"likelihood\"`.")
    }
    return(assign_individuals_ml(
      data = data, strata = strata, engine = engine, folds = folds,
      random.seed = random.seed, verbose = verbose
    ))
  }
  if (!is.logical(leave.one.out) || length(leave.one.out) != 1L || is.na(leave.one.out)) {
    rlang::abort("`leave.one.out` must be TRUE or FALSE.")
  }
  if (!is.null(frequency.floor) &&
      (!is.numeric(frequency.floor) || length(frequency.floor) != 1L ||
       is.na(frequency.floor) || frequency.floor <= 0 || frequency.floor >= 0.5)) {
    rlang::abort("`frequency.floor` must be NULL or one number between 0 and 0.5.")
  }
  if (!is.numeric(effective.n.tolerance) || length(effective.n.tolerance) != 1L ||
      is.na(effective.n.tolerance) || effective.n.tolerance < 0 ||
      effective.n.tolerance >= 1) {
    rlang::abort("`effective.n.tolerance` must be one number from 0 up to, but not including, 1.")
  }
  if (!is.numeric(em.tolerance) || length(em.tolerance) != 1L ||
      is.na(em.tolerance) || em.tolerance <= 0) {
    rlang::abort("`em.tolerance` must be one positive number.")
  }
  if (!is.numeric(em.max.iterations) || length(em.max.iterations) != 1L ||
      is.na(em.max.iterations) || em.max.iterations < 1) {
    rlang::abort("`em.max.iterations` must be one positive integer.")
  }
  em.max.iterations <- as.integer(em.max.iterations)
  if (!is.logical(unsampled.source) || length(unsampled.source) != 1L ||
      is.na(unsampled.source)) {
    rlang::abort("`unsampled.source` must be TRUE or FALSE.")
  }
  if (!is.numeric(unsampled.cutoff) || length(unsampled.cutoff) != 1L ||
      is.na(unsampled.cutoff)) {
    rlang::abort("`unsampled.cutoff` must be one number.")
  }
  if (!is.numeric(sequencing.error) || length(sequencing.error) != 1L ||
      is.na(sequencing.error) || sequencing.error <= 0 || sequencing.error >= 0.5) {
    rlang::abort("`sequencing.error` must be one number between 0 and 0.5.")
  }

  genome <- genometranslator::read_genome(data = data, import.metadata = TRUE)
  required <- c("INDIVIDUALS", "MARKERS")
  missing.columns <- setdiff(required, names(genome))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The genomic data is missing required columns: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }

  if (genotype.method != "GT") {
    result <- assign_individuals_gl(
      genome = genome, strata = strata, genotype.method = genotype.method,
      leave.one.out = leave.one.out, frequency.floor = frequency.floor,
      effective.n.tolerance = effective.n.tolerance,
      em.tolerance = em.tolerance, em.max.iterations = em.max.iterations,
      unsampled.source = unsampled.source,
      unsampled.cutoff = unsampled.cutoff,
      sequencing.error = sequencing.error,
      evaluation.data = evaluation.data,
      marker.blocks = marker.blocks,
      verbose = verbose
    )
    if (!is.null(reporting.unit)) {
      metadata <- if (is.null(strata)) genome else
        genometranslator::read_strata(strata = strata)$strata
      if (!reporting.unit %in% names(metadata)) {
        rlang::abort(paste0("Reporting-unit column `", reporting.unit, "` was not found."))
      }
      reporting.strata <- dplyr::transmute(
        metadata, INDIVIDUALS,
        STRATA = as.character(.data[[reporting.unit]])
      ) |>
        dplyr::distinct()
      result$reporting.unit <- assign_individuals(
        data = genome, strata = reporting.strata, genotype.method = genotype.method,
        leave.one.out = leave.one.out, frequency.floor = frequency.floor,
        effective.n.tolerance = effective.n.tolerance,
        em.tolerance = em.tolerance, em.max.iterations = em.max.iterations,
        unsampled.source = FALSE, evaluation.data = evaluation.data,
        marker.blocks = marker.blocks, reporting.unit = NULL, verbose = FALSE
      )
    }
    return(result)
  }

  if (!is.null(evaluation.data)) {
    rlang::abort("`evaluation.data` is currently available for GL/PL assignment only.")
  }

  if ("ALT_DOSAGE" %in% names(genome)) {
    genome <- dplyr::mutate(genome, GT = as.integer(ALT_DOSAGE))
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | as.character(genome$GT) %in% c("0", "1", "2"))) {
    genome <- dplyr::mutate(genome, GT = as.integer(GT))
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | grepl("^[0-9]{6}$", as.character(genome$GT)))) {
    gt.character <- as.character(genome$GT)
    genome$A1 <- substr(gt.character, 1L, 3L)
    genome$A2 <- substr(gt.character, 4L, 6L)
    genome$A1[gt.character == "000000" | is.na(gt.character)] <- NA_character_
    genome$A2[gt.character == "000000" | is.na(gt.character)] <- NA_character_
    marker.alleles <- genome |>
      dplyr::select(MARKERS, A1, A2) |>
      tidyr::pivot_longer(c(A1, A2), values_to = "ALLELE") |>
      dplyr::filter(!is.na(ALLELE)) |>
      dplyr::group_by(MARKERS) |>
      dplyr::summarise(
        REF_ALLELE = sort(unique(ALLELE))[1L],
        N_ALLELES = dplyr::n_distinct(ALLELE),
        .groups = "drop"
      )
    if (any(marker.alleles$N_ALLELES > 2L)) {
      rlang::abort("The native likelihood engine currently requires biallelic markers.")
    }
    genome <- genome |>
      dplyr::left_join(marker.alleles, by = "MARKERS") |>
      dplyr::mutate(
        GT = dplyr::if_else(
          is.na(A1) | is.na(A2), NA_integer_,
          as.integer(A1 != REF_ALLELE) + as.integer(A2 != REF_ALLELE)
        )
      ) |>
      dplyr::select(-A1, -A2, -REF_ALLELE, -N_ALLELES)
  } else {
    rlang::abort(paste0(
      "The native likelihood engine requires `ALT_DOSAGE`, numeric dosage in ",
      "`GT`, or six-digit biallelic genotypes in `GT`."
    ))
  }

  if (!is.null(strata)) {
    strata.data <- genometranslator::read_strata(strata = strata) $strata
    genome <- dplyr::select(genome, -tidyselect::any_of("STRATA")) |>
      dplyr::left_join(
        dplyr::select(strata.data, INDIVIDUALS, STRATA),
        by = "INDIVIDUALS"
      )
  }
  if (!"STRATA" %in% names(genome)) {
    rlang::abort("Supply `strata`, or include a `STRATA` column in `data`.")
  }
  if (any(!is.na(genome$GT) & !genome$GT %in% 0:2)) {
    rlang::abort("`GT` must contain diploid alternate-allele dosages: 0, 1, 2, or NA.")
  }

  individuals <- dplyr::distinct(genome, INDIVIDUALS, STRATA)
  duplicated.ids <- individuals |>
    dplyr::count(INDIVIDUALS, name = "N") |>
    dplyr::filter(N > 1L)
  if (nrow(duplicated.ids)) {
    rlang::abort("Each individual must belong to at most one reference stratum.")
  }

  reference <- dplyr::filter(genome, !is.na(STRATA), nzchar(as.character(STRATA)))
  if (!nrow(reference)) rlang::abort("No reference individuals have a valid `STRATA` value.")

  allele.frequencies <- reference |>
    dplyr::group_by(STRATA, MARKERS) |>
    dplyr::summarise(
      ALT_COUNT = sum(GT, na.rm = TRUE),
      N_GENOTYPED = sum(!is.na(GT)),
      ALT_FREQ = dplyr::if_else(N_GENOTYPED > 0L, ALT_COUNT / (2 * N_GENOTYPED), NA_real_),
      .groups = "drop"
    )

  candidate.strata <- sort(unique(as.character(allele.frequencies$STRATA)))
  common.frequency.markers <- allele.frequencies |>
    dplyr::group_by(MARKERS) |>
    dplyr::summarise(
      N_CANDIDATES = sum(N_GENOTYPED > 0L),
      .groups = "drop"
    ) |>
    dplyr::filter(N_CANDIDATES == length(candidate.strata)) |>
    dplyr::pull(MARKERS)

  genotype.data <- dplyr::select(
    genome, INDIVIDUALS, CURRENT = STRATA, MARKERS, GT
  ) |>
    dplyr::filter(MARKERS %in% common.frequency.markers)
  marker.blocks <- normalize_marker_blocks(marker.blocks, common.frequency.markers)

  if (leave.one.out) {
    home.counts <- allele.frequencies |>
      dplyr::select(CURRENT = STRATA, MARKERS, HOME_N_GENOTYPED = N_GENOTYPED)
    genotype.data <- genotype.data |>
      dplyr::left_join(home.counts, by = c("CURRENT", "MARKERS")) |>
      dplyr::filter(
        is.na(CURRENT) | is.na(GT) | is.na(HOME_N_GENOTYPED) |
          HOME_N_GENOTYPED > 1L
      ) |>
      dplyr::select(-HOME_N_GENOTYPED)
  }

  score_candidate <- function(candidate) {
    frequencies <- dplyr::filter(allele.frequencies, as.character(STRATA) == candidate) |>
      dplyr::select(MARKERS, ALT_COUNT, N_GENOTYPED, ALT_FREQ)

    scored <- dplyr::left_join(genotype.data, frequencies, by = "MARKERS")
    use.loo <- leave.one.out & !is.na(scored$CURRENT) &
      as.character(scored$CURRENT) == candidate & !is.na(scored$GT)
    alt.count <- scored$ALT_COUNT - ifelse(use.loo, scored$GT, 0)
    n.genotyped <- scored$N_GENOTYPED - ifelse(use.loo, 1L, 0L)
    p <- ifelse(n.genotyped > 0L, alt.count / (2 * n.genotyped), NA_real_)
    floor.value <- if (is.null(frequency.floor)) {
      ifelse(n.genotyped > 0L, 1 / (2 * n.genotyped + 1), NA_real_)
    } else {
      rep(frequency.floor, length(p))
    }
    p <- pmin(1 - floor.value, pmax(floor.value, p))

    log.probability <- rep(NA_real_, length(p))
    valid <- !is.na(scored$GT) & !is.na(p)
    log.probability[valid & scored$GT == 0] <- 2 * log1p(-p[valid & scored$GT == 0])
    log.probability[valid & scored$GT == 1] <- log(2) + log(p[valid & scored$GT == 1]) +
      log1p(-p[valid & scored$GT == 1])
    log.probability[valid & scored$GT == 2] <- 2 * log(p[valid & scored$GT == 2])

    scored$LOG_PROBABILITY <- log.probability
    q0 <- (1 - p)^2
    q1 <- 2 * p * (1 - p)
    q2 <- p^2
    locus.mean <- q0 * log(q0) + q1 * log(q1) + q2 * log(q2)
    locus.var <- q0 * (log(q0) - locus.mean)^2 +
      q1 * (log(q1) - locus.mean)^2 + q2 * (log(q2) - locus.mean)^2
    locus.mean[!valid] <- NA_real_
    locus.var[!valid] <- NA_real_
    scored$EXPECTED_MEAN <- locus.mean
    scored$EXPECTED_VAR <- locus.var
    summary <- scored |>
      dplyr::group_by(INDIVIDUALS, CURRENT) |>
      dplyr::summarise(
        LOG_LIKELIHOOD = if (all(is.na(LOG_PROBABILITY))) NA_real_ else sum(LOG_PROBABILITY, na.rm = TRUE),
        N_MARKERS = sum(!is.na(LOG_PROBABILITY)),
        EXPECTED_MEAN = sum(EXPECTED_MEAN, na.rm = TRUE),
        EXPECTED_VAR = sum(EXPECTED_VAR, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        CANDIDATE = candidate,
        Z_LIKELIHOOD = dplyr::if_else(
          EXPECTED_VAR > 0,
          (LOG_LIKELIHOOD - EXPECTED_MEAN) / sqrt(EXPECTED_VAR), NA_real_
        )
      )
    per.marker <- if (is.null(marker.blocks)) NULL else
      dplyr::transmute(
        scored, INDIVIDUALS, CURRENT, MARKERS,
        CANDIDATE = candidate, LOG_PROBABILITY
      )
    list(summary = summary, per.marker = per.marker)
  }

  if (isTRUE(verbose)) {
    message("Calculating likelihoods for ", length(candidate.strata), " candidate strata")
  }
  scored.candidates <- purrr::map(candidate.strata, score_candidate)
  likelihoods <- purrr::map(scored.candidates, "summary") |>
    purrr::list_rbind() |>
    dplyr::select(
      INDIVIDUALS, CURRENT, CANDIDATE, LOG_LIKELIHOOD, N_MARKERS,
      EXPECTED_MEAN, EXPECTED_VAR, Z_LIKELIHOOD
    ) |>
    dplyr::arrange(INDIVIDUALS, dplyr::desc(LOG_LIKELIHOOD), CANDIDATE)
  marker.block.likelihoods <- if (is.null(marker.blocks)) NULL else {
    purrr::map(scored.candidates, "per.marker") |>
      purrr::list_rbind() |>
      summarize_block_likelihoods(marker.blocks)
  }

  ranked <- likelihoods |>
    dplyr::group_by(INDIVIDUALS) |>
    dplyr::arrange(dplyr::desc(LOG_LIKELIHOOD), CANDIDATE, .by_group = TRUE) |>
    dplyr::mutate(RANK = dplyr::row_number()) |>
    dplyr::ungroup()

  best <- ranked |>
    dplyr::filter(RANK == 1L) |>
    dplyr::transmute(
      INDIVIDUALS, CURRENT, INFERRED = CANDIDATE,
      LOG_LIK_MAX = LOG_LIKELIHOOD, N_MARKERS
    )
  second <- ranked |>
    dplyr::filter(RANK == 2L) |>
    dplyr::transmute(
      INDIVIDUALS, SECOND_BEST = CANDIDATE,
      LOG_LIK_SECOND = LOG_LIKELIHOOD
    )
  home <- likelihoods |>
    dplyr::filter(!is.na(CURRENT), as.character(CURRENT) == CANDIDATE) |>
    dplyr::transmute(INDIVIDUALS, LOG_LIK_HOME = LOG_LIKELIHOOD)

  assignment <- best |>
    dplyr::left_join(second, by = "INDIVIDUALS") |>
    dplyr::left_join(home, by = "INDIVIDUALS") |>
    dplyr::mutate(
      LOG_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_SECOND,
      HOME_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_HOME,
      CORRECT = dplyr::if_else(
        is.na(CURRENT), NA, as.character(CURRENT) == INFERRED
      )
    ) |>
    dplyr::select(
      INDIVIDUALS, CURRENT, INFERRED, CORRECT,
      LOG_LIK_HOME, LOG_LIK_MAX, HOME_LIK_RATIO,
      SECOND_BEST, LOG_LIK_SECOND, LOG_LIK_RATIO, N_MARKERS
    ) |>
    dplyr::arrange(CURRENT, INDIVIDUALS)

  result <- list(
    assignment = assignment,
    likelihoods = likelihoods,
    allele.frequencies = allele.frequencies,
    marker.block.likelihoods = marker.block.likelihoods
  )
  if (!is.null(reporting.unit)) {
    metadata <- if (is.null(strata)) genome else
      genometranslator::read_strata(strata = strata)$strata
    if (!reporting.unit %in% names(metadata)) {
      rlang::abort(paste0("Reporting-unit column `", reporting.unit, "` was not found."))
    }
    reporting.strata <- dplyr::transmute(
      metadata, INDIVIDUALS, STRATA = as.character(.data[[reporting.unit]])
    ) |>
      dplyr::distinct()
    result$reporting.unit <- assign_individuals(
      data = genome, strata = reporting.strata, genotype.method = "GT",
      leave.one.out = leave.one.out, frequency.floor = frequency.floor,
      marker.blocks = marker.blocks, reporting.unit = NULL, verbose = FALSE
    )
  }
  result
}
