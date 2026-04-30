# ==============================================================================
# MASTER THESIS: THE DISTRIBUTIONAL CONSEQUENCES OF HIGHER EDUCATION EXPANSION
# SCRIPT: thesis.R
# ==============================================================================
# Author:      Baha Khammari
# Date:        2026-04


rm(list = ls())
gc()

# ==============================================================================
# SECTION 0: SETUP -- PARAMETERS, LIBRARIES, DIRECTORIES, THEME
# ==============================================================================

params <- list(
  age_min        = 22L,
  age_max        = 35L,
  lag_years      = 4L,
  q_low          = 0.20,
  q_mid          = 0.50,
  q_high         = 0.75,
  set_seed       = 7364599L,
  years          = 2005L:2019L,
  central_states = c("SP", "RJ", "MG", "ES", "RS", "SC", "PR")
)

set.seed(params$set_seed)

library(basedosdados)
library(dplyr)
library(tidyr)
library(ggplot2)
library(did)
library(readr)
library(xtable)
library(scales)
library(broom)

for (pkg in c("fixest", "sf", "geobr")) {
  tryCatch(library(pkg, character.only = TRUE),
           error = function(e) message(sprintf("  WARNING: package '%s' not available.", pkg)))
}

tryCatch(library(HonestDiD),
         error = function(e) {
           message("  WARNING: HonestDiD not installed. Attempting GitHub install...")
           tryCatch({
             if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
             devtools::install_github("asheshrambachan/HonestDiD")
             library(HonestDiD)
           }, error = function(e2) message("  WARNING: HonestDiD install failed -- sensitivity section skipped."))
         })

basedosdados::set_billing_id("microdadosbrazil")

dirs <- c("Figures_Final", "Tables_Final", "Documentation_Final")
sapply(dirs, function(x) if (!dir.exists(x)) dir.create(x))

theme_set(theme_classic() + theme(
  plot.title         = element_text(face = "bold", size = 14, color = "black"),
  plot.subtitle      = element_text(size = 11, color = "grey40"),
  axis.text          = element_text(size = 10, color = "black"),
  legend.position    = "bottom",
  panel.grid.major.y = element_line(color = "grey95")
))

# Dose-bin color palette -- used by all three bin histograms and violin plot
DOSE_COLORS <- c(Low = "#56B4E9", Mid = "#E69F00", High = "#D55E00")

# 95% normal-approximation critical value (two-sided)
Z_95 <- qnorm(0.975)

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

safe_err_msg <- function(e) gsub("\\{", "{{", gsub("\\}", "}}", conditionMessage(e)))

fixest_pval <- function(fit) coeftable(fit)[, "Pr(>|t|)"]

add_micro_names <- function(df, geo) {
  cols_needed <- c("nome_microrregiao", "sigla_uf")
  cols_missing <- setdiff(cols_needed, names(df))
  if (length(cols_missing) == 0) return(df)
  lkp <- geo %>%
    distinct(id_microrregiao, .keep_all = TRUE) %>%
    select(id_microrregiao, all_of(cols_missing))
  df %>% left_join(lkp, by = "id_microrregiao")
}

# Wald chi-squared pre-trend joint test on ATT_ES(e < 0).
# FIX v5: exhaustive VCov field search; falls back to diagonal SEs if not found.
pretrend_ftest <- function(att_gt_obj) {
  es      <- did::aggte(att_gt_obj, type = "dynamic", na.rm = TRUE)
  pre_idx <- which(es$egt < 0)
  if (length(pre_idx) == 0) {
    warning("No pre-treatment periods available for Wald test.")
    return(list(statistic = NA_real_, df = 0L, p.value = NA_real_,
                used_full_vcov = FALSE))
  }
  ests <- es$att.egt[pre_idx]

  # v5 FIX: search all fields for a matrix of the right dimension
  vcv <- NULL
  used_full_vcov <- FALSE
  for (nm in names(es)) {
    mat <- es[[nm]]
    if (is.matrix(mat) && nrow(mat) == length(es$att.egt) &&
        ncol(mat) == length(es$att.egt)) {
      vcv <- mat[pre_idx, pre_idx, drop = FALSE]
      used_full_vcov <- TRUE
      break
    }
  }
  if (is.null(vcv) || any(is.na(vcv))) {
    ses <- es$se.egt[pre_idx]
    W   <- sum((ests / ses)^2, na.rm = TRUE)
    df  <- sum(!is.na(ests / ses))
  } else {
    W  <- tryCatch(as.numeric(t(ests) %*% solve(vcv) %*% ests),
                   error = function(e) { warning("VCov inversion failed."); NA_real_ })
    df <- length(ests)
  }
  list(statistic = W, df = df,
       p.value = pchisq(W, df = df, lower.tail = FALSE),
       used_full_vcov = used_full_vcov)
}

#
# run_cs -- wrapper on did::att_gt that:
#   1) aggregates demographic cells to the microregion x year x cohort panel;
#   2) left-joins pre-treatment covariates from df_covars_dr (if supplied);
#   3) calls did::att_gt with `est_method = "dr"` so that covariate-adjusted
#      specifications use the Sant'Anna & Zhao (2020) doubly-robust estimator
#      (propensity-score reweighting x outcome-regression adjustment).  When
#      `xformla = ~1` the DR estimator degenerates to the unadjusted comparison
#      of not-yet-treated controls, which is identical to est_method = "reg"
#      with no covariates -- see the `did` package documentation.
#
run_cs <- function(d, outcome_var, label = "", agg = c("wmean", "logsum"),
                   xformla = ~1, covars_df = NULL,
                   est_method = c("dr", "reg", "ipw")) {
  agg        <- match.arg(agg)
  est_method <- match.arg(est_method)

  d_agg <- if (agg == "wmean") {
    d %>%
      group_by(id_num, ano, g) %>%
      summarise(y = weighted.mean(.data[[outcome_var]], w = n_vinculos, na.rm = TRUE),
                .groups = "drop")
  } else {
    d %>%
      group_by(id_num, ano, g) %>%
      summarise(y = log(sum(.data[[outcome_var]], na.rm = TRUE)),
                .groups = "drop")
  }

  if (!is.null(covars_df)) {
    d_agg <- d_agg %>% left_join(covars_df, by = "id_num")
  }

  n_treated <- d_agg %>% filter(g > 0) %>% pull(id_num) %>% n_distinct()
  n_notyet  <- d_agg %>% filter(g == 0) %>% pull(id_num) %>% n_distinct()
  n_obs     <- nrow(d_agg)

  if (n_treated < 20) {
    warning(sprintf("[%s] Only %d treated units (need >= 20). Skipping.", label, n_treated))
    return(NULL)
  }
  if (n_notyet < 20) {
    warning(sprintf("[%s] Only %d not-yet-treated units (need >= 20). Skipping.", label, n_notyet))
    return(NULL)
  }
  message(sprintf("   [%s] %d treated | %d not-yet-treated | %d obs | est_method = %s",
                  label, n_treated, n_notyet, n_obs, est_method))

  att <- did::att_gt(
    yname                  = "y",
    tname                  = "ano",
    idname                 = "id_num",
    gname                  = "g",
    xformla                = xformla,
    data                   = d_agg,
    control_group          = "notyettreated",
    allow_unbalanced_panel = TRUE,
    est_method             = est_method,
    print_details          = FALSE
  )

  es         <- did::aggte(att, type = "dynamic", na.rm = TRUE)
  es_tidy    <- broom::tidy(es)
  att_simple <- did::aggte(att, type = "simple", na.rm = TRUE)
  pretrend   <- pretrend_ftest(att)

  message(sprintf("   [%s] ATT = %.4f (SE = %.4f) | Pre-trend Wald p = %.4f",
                  label, att_simple$overall.att, att_simple$overall.se,
                  pretrend$p.value))

  list(att_gt = att, es = es, tidy = es_tidy, simple = att_simple,
       pretrend = pretrend, label = label,
       n_obs = n_obs, n_treated = n_treated, n_notyet = n_notyet)
}

build_att_row <- function(res) {
  s   <- res$simple
  t_val <- s$overall.att / s$overall.se
  p_val <- 2 * pnorm(-abs(t_val))
  tibble(
    Model      = res$label,
    ATT        = round(s$overall.att, 4),
    SE         = round(s$overall.se,  4),
    t_val      = round(t_val, 3),
    p_val      = round(p_val, 4),
    CI_low     = round(s$overall.att - 1.96 * s$overall.se, 4),
    CI_high    = round(s$overall.att + 1.96 * s$overall.se, 4),
    Pretrend_p = round(res$pretrend$p.value, 4),
    N_obs      = res$n_obs,
    N_treated  = res$n_treated,
    N_notyet   = res$n_notyet
  )
}

plot_event_study <- function(df1, df2, l1, l2, c1, c2, title_txt,
                             ylab = "ATT on log(wages in min. wage units)") {
  ggplot() +
    geom_hline(yintercept = 0, color = "black",  linetype = "solid",  linewidth = 0.4) +
    geom_vline(xintercept = -1, color = "grey40", linetype = "dashed", linewidth = 0.5) +
    geom_ribbon(data = df2, aes(x = event.time, ymin = conf.low, ymax = conf.high),
                alpha = 0.12, fill = c2) +
    geom_line(data = df2, aes(x = event.time, y = estimate, color = l2), linewidth = 1.1) +
    geom_point(data = df2, aes(x = event.time, y = estimate, color = l2), size = 1.8) +
    geom_ribbon(data = df1, aes(x = event.time, ymin = conf.low, ymax = conf.high),
                alpha = 0.20, fill = c1) +
    geom_line(data = df1, aes(x = event.time, y = estimate, color = l1), linewidth = 1.1) +
    geom_point(data = df1, aes(x = event.time, y = estimate, color = l1), size = 1.8) +
    scale_color_manual(name = "", values = setNames(c(c1, c2), c(l1, l2))) +
    scale_x_continuous(breaks = seq(-12, 10, 2), limits = c(NA, 10)) +
    labs(title = title_txt,
         subtitle = "95% confidence bands | CS (2021), not-yet-treated comparison",
         x = "Event time relative to first treatment year", y = ylab)
}

run_cs_preagg <- function(d, label) {
  n_treated <- d %>% filter(g > 0) %>% pull(id_num) %>% n_distinct()
  n_notyet  <- d %>% filter(g == 0) %>% pull(id_num) %>% n_distinct()

  if (n_treated < 20 || n_notyet < 20) {
    warning(sprintf("[%s] Insufficient units: %d treated, %d not-yet-treated",
                    label, n_treated, n_notyet))
    return(NULL)
  }

  message(sprintf("   [%s] %d treated | %d not-yet-treated | %d obs",
                  label, n_treated, n_notyet, nrow(d)))

  att <- did::att_gt(
    yname  = "y", tname = "ano", idname = "id_num", gname = "g",
    data   = d,
    control_group = "notyettreated",
    allow_unbalanced_panel = TRUE,
    print_details = FALSE
  )

  es         <- did::aggte(att, type = "dynamic", na.rm = TRUE)
  es_tidy    <- broom::tidy(es)
  att_simple <- did::aggte(att, type = "simple", na.rm = TRUE)
  pretrend   <- pretrend_ftest(att)

  message(sprintf("   [%s] ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                  label, att_simple$overall.att, att_simple$overall.se,
                  pretrend$p.value))

  list(att_gt = att, es = es, tidy = es_tidy, simple = att_simple,
       pretrend = pretrend, label = label, n_obs = nrow(d),
       n_treated = n_treated, n_notyet = n_notyet)
}

# v5 FIX 2: HonestDiD helper with exhaustive VCov search + diagonal fallback
run_honest_v5 <- function(es_obj, label) {
  betahat     <- es_obj$att.egt
  event_times <- es_obj$egt
  post_start  <- min(which(event_times >= 0))
  n_pre       <- post_start - 1L
  n_post      <- length(betahat) - post_start + 1L

  # Guard against interactive debuggers AND warn=2 promotion BEFORE any
  # warnings are issued. A user's .Rprofile may set options(warn = 2) and
  # options(error = recover), which would otherwise trap on the internal
  # lpSolveAPI error path and hang the script.
  old_err  <- getOption("error")
  old_warn <- getOption("warn")
  on.exit({ options(error = old_err); options(warn = old_warn) }, add = TRUE)
  options(error = NULL, warn = 0)

  if (n_pre < 1L || n_post < 1L) {
    message(sprintf("   [%s] Not enough pre/post periods for HonestDiD -- skipped.", label))
    return(list(label = label, results = NULL, breakdown = NA_real_,
                used_full_vcov = FALSE))
  }

  # Exhaustive search for a square matrix of dimension == length(betahat)
  sigma          <- NULL
  used_full_vcov <- FALSE
  for (nm in names(es_obj)) {
    mat <- es_obj[[nm]]
    if (is.matrix(mat) &&
        nrow(mat) == length(betahat) &&
        ncol(mat) == length(betahat) &&
        !any(is.na(mat))) {
      sigma          <- mat
      used_full_vcov <- TRUE
      break
    }
  }

  # Fallback: diagonal matrix from pointwise SEs. Use message() (not warning())
  # so the options(warn = 2) user setting cannot promote it to an error.
  if (is.null(sigma)) {
    message(sprintf("   [%s] Full VCov not found -- using diag(se^2) fallback.", label))
    sigma <- diag(es_obj$se.egt^2)
  }

  tryCatch({
    results <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat        = betahat,
      sigma          = sigma,
      numPrePeriods  = n_pre,
      numPostPeriods = n_post,
      Mbarvec        = seq(0, 2, by = 0.25)
    )
    # Breakdown: largest Mbar for which the lower CI bound still > 0
    pos_rows <- !is.na(results$lb) & results$lb > 0
    bd <- if (any(pos_rows)) max(results$Mbar[pos_rows]) else NA_real_
    message(sprintf("   [%s] HonestDiD breakdown Mbar = %s (full VCov: %s)",
                    label,
                    ifelse(is.na(bd), "NA", sprintf("%.2f", bd)),
                    used_full_vcov))
    list(label = label, results = results, breakdown = bd,
         used_full_vcov = used_full_vcov)
  }, error = function(e) {
    message(sprintf("   [%s] HonestDiD failed: %s", label, safe_err_msg(e)))
    list(label = label, results = NULL, breakdown = NA_real_,
         used_full_vcov = used_full_vcov)
  })
}

# ==============================================================================
# SECTION 2: DOCUMENTATION & PROVENANCE
# ==============================================================================
message(" [1/12] Fetching Official Dictionaries...")

tryCatch({
  dict_prouni <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_mec_prouni.dicionario`
     WHERE id_tabela = 'microdados'")
  write_csv(dict_prouni, "Documentation_Final/appendix_prouni_codes.csv")
}, error = function(e) message("  Warning: ProUni dict skipped."))

tryCatch({
  dict_rais <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_me_rais.dicionario`
     WHERE id_tabela = 'microdados_vinculos'")
  write_csv(dict_rais, "Documentation_Final/appendix_rais_codes.csv")
}, error = function(e) message("  Warning: RAIS dict skipped."))

tryCatch({
  dict_census <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_ibge_censo_demografico.dicionario`
     WHERE id_tabela = 'microdados_pessoa_2010'")
  write_csv(dict_census, "Documentation_Final/appendix_census_codes.csv")
}, error = function(e) message("  Warning: Census dict skipped."))

# ==============================================================================
# SECTION 3: DATA EXTRACTION FROM BIGQUERY  (v5: RDS caching on every call)
# ==============================================================================
# RDS cache helper -- loads from disk if available, otherwise queries BigQuery.
load_or_fetch <- function(cache_path, query_fn, label) {
  if (file.exists(cache_path)) {
    message(sprintf("   -> Cache HIT: %s", cache_path))
    return(readRDS(cache_path))
  }
  message(sprintf("   -> Cache MISS -- querying BigQuery: %s", label))
  result <- query_fn()
  saveRDS(result, cache_path)
  result
}

message(" [2/12] Extracting Data from BigQuery (with RDS caching)...")

# 3A. Geography crosswalk
df_geo <- load_or_fetch(
  "Documentation_Final/cache_geo.rds",
  function() basedosdados::read_sql("
    SELECT
      id_municipio,
      id_microrregiao,
      nome_microrregiao,
      sigla_uf,
      nome_regiao AS regiao
    FROM `basedosdados.br_bd_diretorios_brasil.municipio`
  "),
  "Geography crosswalk"
)

# 3B. Census 2010 -- youth population denominator + covariates
df_census <- load_or_fetch(
  "Documentation_Final/cache_census.rds",
  function() basedosdados::read_sql("
    SELECT
      id_municipio,
      v0601 AS sexo,
      v0606 AS raca_cor,
      SUM(peso_amostral) AS pop_18_24,
      SUM(CASE WHEN v0606 IN ('2', '4')
               THEN peso_amostral ELSE 0 END) AS pop_nonwhite,
      SUM(CASE WHEN v6400 = '3'
                AND (v0633 IS NULL OR CAST(v0633 AS INT64) < 11)
               THEN peso_amostral ELSE 0 END) AS pop_18_24_eligible,
      AVG(CASE WHEN CAST(v6511 AS FLOAT64) = 0 OR v6511 IS NULL
               THEN 1.0 ELSE 0.0 END) AS share_no_formal_income,
      AVG(v6532) AS avg_hh_income_pc_mw,
      AVG(CASE WHEN v0627 = '1' THEN 1.0 ELSE 0.0 END) AS literacy_rate,
      AVG(CASE WHEN v0628 = '1' THEN 1.0 ELSE 0.0 END) AS share_in_school,
      AVG(CASE WHEN v0616 IN ('3','4') THEN 1.0 ELSE 0.0 END) AS share_mobility_disability
    FROM `basedosdados.br_ibge_censo_demografico.microdados_pessoa_2010`
    WHERE
      v6036 BETWEEN 18 AND 24
      AND v0601 NOT IN ('9')
      AND v0606 NOT IN ('9')
    GROUP BY id_municipio, v0601, v0606
  "),
  "Census 2010"
)

# 3C. ProUni scholarships -- batched download (now also keeps id_ies for HEI linkage).
# v5.1 schema change: if the cached RDS predates the id_ies addition, invalidate it.
if (file.exists("Documentation_Final/cache_prouni.rds")) {
  .tmp_prouni_cache <- tryCatch(readRDS("Documentation_Final/cache_prouni.rds"),
                                error = function(e) NULL)
  if (!is.null(.tmp_prouni_cache) && !"id_ies" %in% names(.tmp_prouni_cache)) {
    message("   -> ProUni cache predates id_ies column; invalidating for refresh.")
    file.rename("Documentation_Final/cache_prouni.rds",
                "Documentation_Final/cache_prouni_legacy_no_ies.rds")
  }
  rm(.tmp_prouni_cache)
}
message("   -> ProUni (batched 2005-2009 / 2010-2014 / 2015-2019)...")
df_prouni_raw <- load_or_fetch(
  "Documentation_Final/cache_prouni.rds",
  function() {
    q_base  <- "
      SELECT ano, id_municipio, id_ies, raca_cor, sexo, tipo_bolsa,
             modalidade_ensino, nivel_bolsa,
             COUNT(*) AS n_bolsas
      FROM `basedosdados.br_mec_prouni.microdados`
    "
    q_group <- "GROUP BY ano, id_municipio, id_ies, raca_cor, sexo, tipo_bolsa,
                         modalidade_ensino, nivel_bolsa"
    safe_pull <- function(where_clause) {
      tryCatch(basedosdados::read_sql(paste(q_base, where_clause, q_group)),
               error = function(e) {
                 message(sprintf("   ProUni SQL w/ id_ies failed (%s) -- falling back to legacy schema",
                                 e$message))
                 q_legacy_base  <- "
                   SELECT ano, id_municipio, raca_cor, sexo, tipo_bolsa,
                          COUNT(*) AS n_bolsas
                   FROM `basedosdados.br_mec_prouni.microdados`
                 "
                 q_legacy_group <- "GROUP BY ano, id_municipio, raca_cor, sexo, tipo_bolsa"
                 out <- basedosdados::read_sql(paste(q_legacy_base, where_clause, q_legacy_group))
                 out$id_ies             <- NA_character_
                 out$modalidade_ensino  <- NA_character_
                 out$nivel_bolsa        <- NA_character_
                 out
               })
    }
    p1 <- safe_pull("WHERE ano BETWEEN 2005 AND 2009")
    p2 <- safe_pull("WHERE ano BETWEEN 2010 AND 2014")
    p3 <- safe_pull("WHERE ano BETWEEN 2015 AND 2019")
    bind_rows(p1, p2, p3)
  },
  "ProUni scholarships"
)

# 3D. RAIS outcomes -- year-by-year loop
message("   -> RAIS main (year-by-year loop)...")
df_rais <- load_or_fetch(
  "Documentation_Final/cache_rais_main.rds",
  function() {
    rais_list <- list()
    for (y in params$years) {
      q <- paste0("
        SELECT
          ano, id_municipio, raca_cor, sexo,
          COUNT(*) AS n_vinculos,
          AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1))) AS log_wage_sm
        FROM `basedosdados.br_me_rais.microdados_vinculos`
        WHERE
          ano = ", y, "
          AND idade BETWEEN ", params$age_min, " AND ", params$age_max, "
          AND valor_remuneracao_media_sm > 0
          AND tipo_vinculo = '10'
          AND vinculo_ativo_3112 = '1'
        GROUP BY ano, id_municipio, raca_cor, sexo
      ")
      tryCatch({
        rais_list[[as.character(y)]] <- basedosdados::read_sql(q)
        message(sprintf("     Year %d OK", y))
      }, error = function(e) {
        message(sprintf("     Error fetching year %d: %s", y, safe_err_msg(e)))
      })
      gc()
    }
    missing_years <- setdiff(as.character(params$years), names(rais_list))
    if (length(missing_years) > 0) {
      stop(sprintf("RAIS download incomplete -- missing years: %s",
                   paste(missing_years, collapse = ", ")))
    }
    bind_rows(rais_list)
  },
  "RAIS main"
)
gc()

# 3E. RAIS age-cohort split (Duflo-style)
message(" [NEW] RAIS age-cohort split (young 22-28 / old 29-35)...")
df_rais_cohort <- load_or_fetch(
  "Documentation_Final/cache_rais_cohort.rds",
  function() {
    cohort_list <- list()
    for (y in params$years) {
      q <- paste0("
        SELECT
          ano, id_municipio, raca_cor, sexo,
          CASE
            WHEN idade BETWEEN 22 AND 28 THEN 'young'
            WHEN idade BETWEEN 29 AND 35 THEN 'old'
          END AS age_group,
          COUNT(*) AS n_vinculos,
          AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1))) AS log_wage_sm
        FROM `basedosdados.br_me_rais.microdados_vinculos`
        WHERE
          ano = ", y, "
          AND idade BETWEEN 22 AND 35
          AND valor_remuneracao_media_sm > 0
          AND tipo_vinculo = '10'
          AND vinculo_ativo_3112 = '1'
        GROUP BY ano, id_municipio, raca_cor, sexo, age_group
      ")
      tryCatch({
        cohort_list[[as.character(y)]] <- basedosdados::read_sql(q)
        message(sprintf("     Cohort year %d OK", y))
      }, error = function(e)
        message(sprintf("     Error cohort year %d: %s", y, safe_err_msg(e))))
      gc()
    }
    bind_rows(cohort_list)
  },
  "RAIS cohort"
)
message(sprintf("   Cohort RAIS: %d rows", nrow(df_rais_cohort)))

# 3F. RAIS with education level (first stage for LATE)
message("\n SECTION 3F: RAIS with education level...")
df_rais_educ <- load_or_fetch(
  "Documentation_Final/cache_rais_educ.rds",
  function() {
    educ_list <- list()
    for (y in params$years) {
      q <- paste0("
        SELECT
          ano, id_municipio, raca_cor,
          CASE
            WHEN idade BETWEEN 22 AND 28 THEN 'young'
            WHEN idade BETWEEN 29 AND 35 THEN 'old'
          END AS age_group,
          COUNT(*) AS n_vinculos,
          SUM(CASE WHEN CAST(grau_instrucao_apos_2005 AS INT64) >= 9
                   THEN 1 ELSE 0 END) AS n_college
        FROM `basedosdados.br_me_rais.microdados_vinculos`
        WHERE
          ano = ", y, "
          AND idade BETWEEN 22 AND 35
          AND valor_remuneracao_media_sm > 0
          AND tipo_vinculo = '10'
          AND vinculo_ativo_3112 = '1'
        GROUP BY ano, id_municipio, raca_cor, age_group
      ")
      tryCatch({
        educ_list[[as.character(y)]] <- basedosdados::read_sql(q)
        message(sprintf("     Educ year %d OK", y))
      }, error = function(e)
        message(sprintf("     ERROR educ year %d: %s", y, safe_err_msg(e))))
      gc()
    }
    bind_rows(educ_list)
  },
  "RAIS education"
)
message(sprintf("   Education RAIS: %d rows", nrow(df_rais_educ)))

# 3G. RAIS birth-cohort design (exposed 1987-93 vs control 1975-82)
message("\n SECTION 3G: RAIS birth-cohort data...")
df_rais_bc <- load_or_fetch(
  "Documentation_Final/cache_rais_bc.rds",
  function() {
    tryCatch({
      q_bc <- paste0("
        SELECT
          ano, id_municipio, raca_cor, sexo,
          CASE
            WHEN (ano - idade) BETWEEN 1987 AND 1993 THEN 'exposed'
            WHEN (ano - idade) BETWEEN 1975 AND 1982 THEN 'control'
          END AS birth_cohort,
          COUNT(*) AS n_vinculos,
          AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1))) AS log_wage_sm,
          SUM(CASE WHEN CAST(grau_instrucao_apos_2005 AS INT64) >= 9
                   THEN 1 ELSE 0 END) AS n_college
        FROM `basedosdados.br_me_rais.microdados_vinculos`
        WHERE
          ano BETWEEN ", min(params$years), " AND ", max(params$years), "
          AND idade BETWEEN 18 AND 51
          AND valor_remuneracao_media_sm > 0
          AND tipo_vinculo = '10'
          AND vinculo_ativo_3112 = '1'
          AND (
            (ano - idade) BETWEEN 1975 AND 1982
            OR (ano - idade) BETWEEN 1987 AND 1993
          )
        GROUP BY ano, id_municipio, raca_cor, sexo, birth_cohort
      ")
      basedosdados::read_sql(q_bc)
    }, error = function(e) {
      message(sprintf("   ERROR birth-cohort query: %s", safe_err_msg(e)))
      NULL
    })
  },
  "RAIS birth-cohort"
)
if (!is.null(df_rais_bc))
  message(sprintf("   Birth-cohort RAIS: %d rows", nrow(df_rais_bc)))
gc()

# ---- DESCRIPTIVE CHECKS: loaded datasets ----------------------------------------
message("\n=== DATASET DESCRIPTIVE CHECKS ===")

message(sprintf(
  "  df_geo:    %d rows | %d municipalities | %d microregions | %d UFs",
  nrow(df_geo), n_distinct(df_geo$id_municipio),
  n_distinct(df_geo$id_microrregiao), n_distinct(df_geo$sigla_uf)
))

message(sprintf(
  "  df_census: %d rows | pop_18_24 range [%.0f, %.0f] | avg_hh_income range [%.2f, %.2f]",
  nrow(df_census),
  min(df_census$pop_18_24,          na.rm = TRUE),
  max(df_census$pop_18_24,          na.rm = TRUE),
  min(df_census$avg_hh_income_pc_mw, na.rm = TRUE),
  max(df_census$avg_hh_income_pc_mw, na.rm = TRUE)
))
message(sprintf(
  "             NA counts -- pop_18_24: %d | avg_hh_income_pc_mw: %d | literacy_rate: %d",
  sum(is.na(df_census$pop_18_24)),
  sum(is.na(df_census$avg_hh_income_pc_mw)),
  sum(is.na(df_census$literacy_rate))
))

prouni_total <- sum(df_prouni_raw$n_bolsas, na.rm = TRUE)
message(sprintf(
  "  df_prouni: %d rows | years %d-%d | %d municipalities | %.0f total scholarships",
  as.integer(nrow(df_prouni_raw)),
  as.integer(min(df_prouni_raw$ano, na.rm = TRUE)),
  as.integer(max(df_prouni_raw$ano, na.rm = TRUE)),
  as.integer(n_distinct(df_prouni_raw$id_municipio)),
  prouni_total
))

rais_stats <- df_rais %>%
  summarise(min_ano = min(ano), max_ano = max(ano),
            min_lw = min(log_wage_sm, na.rm = TRUE),
            max_lw = max(log_wage_sm, na.rm = TRUE),
            na_lw  = sum(is.na(log_wage_sm)))
message(sprintf(
  "  df_rais:   %d rows | years %d-%d | log_wage range [%.2f, %.2f] | NA wages: %d",
  as.integer(nrow(df_rais)),
  as.integer(rais_stats$min_ano), as.integer(rais_stats$max_ano),
  rais_stats$min_lw, rais_stats$max_lw,
  as.integer(rais_stats$na_lw)
))
rm(rais_stats)
message(sprintf(
  "  df_rais_cohort: %d rows | young/old split: %s",
  as.integer(nrow(df_rais_cohort)),
  paste(table(df_rais_cohort$age_group), collapse = " / ")
))
message("=== END DATASET CHECKS ===\n")

# ==============================================================================
# SECTION 4: TREATMENT CONSTRUCTION
# ==============================================================================
# D_{r,t} = 1000 * cumulative scholarships_{r, up to t-4} / youth_pop_{r,2010}
# Dose bins from 2019 cross-section of D_{r,t}:
#   Low = (tau_20, tau_50] | Mid = (tau_50, tau_75] (NEW v5) | High > tau_75
# Treatment timing g = first year D_{r,t} > tau_20; g=0 = not-yet-treated
# ==============================================================================
message(" [3/12] Defining Treatment Dose & Timing...")

# 4.1 Aggregate ProUni to microregion x year
# NOTE on id_municipio: in br_mec_prouni.microdados, id_municipio is the
# BENEFICIARY'S home municipality (where the student lives), NOT the HEI
# municipality.  This gives a demand-side dose: scholarships received by
# residents of each microregion.  It is the appropriate measure for studying
# local labour-market outcomes because (a) RAIS workplace municipality
# correlates strongly with residence municipality, and (b) the intent-to-treat
# interprets as "local youth exposed to ProUni".  An HEI-location supply-side
# alternative would require joining id_ies to the IES directory.
df_prouni_agg <- df_prouni_raw %>%
  mutate(ano = as.integer(ano)) %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao, ano) %>%
  summarise(new_schol = sum(n_bolsas, na.rm = TRUE), .groups = "drop")

# 4.2 Youth population denominator at microregion level
df_pop_micro <- df_census %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao) %>%
  summarise(pop = sum(pop_18_24, na.rm = TRUE), .groups = "drop")

# 4.3 Build balanced panel and compute lagged density D_rt
df_panel <- expand_grid(
  id_microrregiao = unique(df_pop_micro$id_microrregiao),
  ano             = params$years
) %>%
  left_join(df_prouni_agg, by = c("id_microrregiao", "ano")) %>%
  mutate(new_schol = replace_na(new_schol, 0L)) %>%
  left_join(df_pop_micro, by = "id_microrregiao") %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(
    cum_schol = cumsum(new_schol),
    D_rt = (dplyr::lag(cum_schol, params$lag_years, default = 0) / pop) * 1000
  ) %>%
  ungroup() %>%
  filter(pop > 0)

# 4.4 Compute dose thresholds from 2019 cross-section of D_rt
df_2019 <- df_panel %>% filter(ano == max(params$years))

tau_20 <- quantile(df_2019$D_rt, params$q_low,  na.rm = TRUE)
tau_50 <- quantile(df_2019$D_rt, params$q_mid,  na.rm = TRUE)
tau_75 <- quantile(df_2019$D_rt, params$q_high, na.rm = TRUE)

message(sprintf("   Thresholds (2019): tau_20=%.4f | tau_50=%.4f | tau_75=%.4f",
                tau_20, tau_50, tau_75))

# 4.5 Classify dose bins and compute treatment timing g
# v5 NEW: Mid bin (tau_50, tau_75] added alongside Low and High
df_final <- df_panel %>%
  left_join(
    df_2019 %>% select(id_microrregiao, D_rt_2019 = D_rt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    dose_bin = case_when(
      D_rt_2019 >  tau_20 & D_rt_2019 <= tau_50 ~ "Low",
      D_rt_2019 >  tau_50 & D_rt_2019 <= tau_75 ~ "Mid",
      D_rt_2019 >  tau_75                         ~ "High",
      TRUE                                         ~ NA_character_
    ),
    g = {
      crossing <- ano[D_rt > tau_20]
      if (length(crossing) == 0) 0L else as.integer(min(crossing))
    }
  ) %>%
  ungroup() %>%
  mutate(g = ifelse(is.infinite(g), 0L, g)) %>%
  add_micro_names(df_geo)

message(sprintf("   Microregions: %d total | %d with g > 0 | %d with g = 0",
                n_distinct(df_final$id_microrregiao),
                n_distinct(df_final$id_microrregiao[df_final$g > 0]),
                n_distinct(df_final$id_microrregiao[df_final$g == 0])))

# ==============================================================================
# SECTION 5: SANITY CHECKS (fail-fast)
# ==============================================================================
message(" [4/12] Running Sanity Checks...")

n_dup <- df_final %>% count(id_microrregiao, ano) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) stop(sprintf("FAIL: %d duplicate (id_microrregiao, ano) pairs.", n_dup))
message("   [OK] No duplicate (id_microrregiao, ano) pairs.")

pre_lag_bad <- df_final %>%
  filter(ano < (min(params$years) + params$lag_years), D_rt != 0) %>%
  nrow()
if (pre_lag_bad > 0)
  stop(sprintf("FAIL: %d rows with non-zero D_rt before %d.",
               as.integer(pre_lag_bad),
               as.integer(min(params$years) + params$lag_years)))
message(sprintf("   [OK] D_rt = 0 for all years before %d.",
                as.integer(min(params$years) + params$lag_years)))

bin_sizes <- df_final %>%
  filter(!is.na(dose_bin)) %>%
  distinct(id_microrregiao, dose_bin) %>%
  count(dose_bin, name = "n_regions")
message("   Bin sizes:")
print(bin_sizes)
for (b in c("Low", "High")) {
  n_b <- bin_sizes$n_regions[bin_sizes$dose_bin == b]
  if (length(n_b) == 0 || n_b <= 50)
    stop(sprintf("FAIL: Bin '%s' has %s microregions (need > 50).",
                 b, ifelse(length(n_b) == 0, "0", as.character(n_b))))
}
message("   [OK] Low and High bins have > 50 microregions.")

non_monotone <- df_final %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(d_diff = D_rt - dplyr::lag(D_rt)) %>%
  filter(!is.na(d_diff), d_diff < -1e-10) %>%
  nrow()
if (non_monotone > 0) {
  warning(sprintf("NOTE: %d rows where D_rt decreased -- check lag logic.", non_monotone))
} else {
  message("   [OK] D_rt is weakly increasing over time.")
}

message(" -> All sanity checks passed.")

# ---- CONSOLE PLOT: Dose distribution histogram ----------------------------------
p_dose_hist <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin)) %>%
  distinct(id_microrregiao, D_rt_2019, dose_bin) %>%
  ggplot(aes(x = D_rt_2019, fill = dose_bin)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.85) +
  geom_vline(xintercept = c(tau_20, tau_50, tau_75),
             linetype = "dashed", color = "grey30", linewidth = 0.7) +
  annotate("text", x = c(tau_20, tau_50, tau_75), y = Inf,
           label = c(sprintf("tau_20\n%.1f", tau_20),
                     sprintf("tau_50\n%.1f", tau_50),
                     sprintf("tau_75\n%.1f", tau_75)),
           hjust = 1.1, vjust = 1.2, size = 2.8, color = "grey30") +
  scale_fill_manual(values = DOSE_COLORS, name = "Dose bin") +
  labs(
    title    = "Dose Distribution: D_rt (2019 cross-section)",
    subtitle = sprintf("Low=%d | Mid=%d | High=%d microregions | Not-yet-treated (g=0): %d",
                       bin_sizes$n_regions[bin_sizes$dose_bin == "Low"],
                       bin_sizes$n_regions[bin_sizes$dose_bin == "Mid"],
                       bin_sizes$n_regions[bin_sizes$dose_bin == "High"],
                       n_distinct(df_final$id_microrregiao[df_final$g == 0])),
    x = "D_rt  (scholarships per 1,000 youth, 4-year lag)",
    y = "Count"
  )
print(p_dose_hist)

# ---- CONSOLE PLOT: Treatment timing histogram -----------------------------------
p_timing_hist <- df_final %>%
  filter(g > 0, !is.na(dose_bin)) %>%
  distinct(id_microrregiao, g, dose_bin) %>%
  ggplot(aes(x = g, fill = dose_bin)) +
  geom_bar(position = "dodge", color = "white") +
  scale_fill_manual(values = DOSE_COLORS, name = "Dose bin") +
  scale_x_continuous(breaks = params$years) +
  labs(
    title    = "Treatment Timing: First year D_rt > tau_20",
    subtitle = "Distribution of first-treatment cohort (g) by dose bin",
    x = "First treatment year (g)", y = "N microregions"
  )
print(p_timing_hist)
message("   Console plots: dose distribution + treatment timing printed.")

# ==============================================================================
# SECTION 6: DESCRIPTIVE STATISTICS & DR COVARIATE SETUP
# ==============================================================================
message(" [5/12] Generating Descriptive Tables...")

# Table 1: Bin summary
tab1 <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin)) %>%
  group_by(dose_bin) %>%
  summarise(
    N_Regions           = n_distinct(id_microrregiao),
    Avg_Youth_Pop       = round(mean(pop, na.rm = TRUE), 0),
    Avg_Density_2019    = round(mean(D_rt_2019, na.rm = TRUE), 2),
    Median_Treatment_Yr = median(g[g > 0], na.rm = TRUE),
    .groups = "drop"
  )
write_csv(tab1, "Tables_Final/tab_bin_summary_v5.csv")
print(xtable(tab1, caption = "Summary Statistics by Dose Bin (Low / Mid / High)"),
      file = "Tables_Final/tab_bin_summary_v5.tex")
message("   Table 1 saved.")

# Aggregate Census to microregion for covariate balance and DR
df_covars <- df_census %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao) %>%
  summarise(
    pop_18_24              = sum(pop_18_24, na.rm = TRUE),
    pop_nonwhite           = sum(pop_nonwhite, na.rm = TRUE),
    share_no_formal_income = mean(share_no_formal_income, na.rm = TRUE),
    avg_hh_income_pc_mw    = mean(avg_hh_income_pc_mw, na.rm = TRUE),
    literacy_rate          = mean(literacy_rate, na.rm = TRUE),
    share_in_school        = mean(share_in_school, na.rm = TRUE),
    .groups = "drop"
  )

# Table 2: Covariate balance
tab2 <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin)) %>%
  distinct(id_microrregiao, dose_bin) %>%
  inner_join(df_covars, by = "id_microrregiao") %>%
  group_by(dose_bin) %>%
  summarise(
    across(
      c(pop_18_24, pop_nonwhite, share_no_formal_income,
        avg_hh_income_pc_mw, literacy_rate, share_in_school),
      list(mean = ~round(mean(.x, na.rm = TRUE), 3),
           sd   = ~round(sd(.x, na.rm = TRUE), 3))
    ),
    .groups = "drop"
  )
write_csv(tab2, "Tables_Final/tab_balance_v5.csv")
print(xtable(tab2, caption = "Pre-Treatment Covariate Balance by Dose Bin"),
      file = "Tables_Final/tab_balance_v5.tex")
message("   Table 2 saved.")

# v5 FIX 3: Enriched DR covariate dataframe
# Low-dose uses 4 covariates (enriched per Fix 3 spec).
# High-dose uses 2 covariates max (avoids singular matrix).
df_covars_dr <- df_covars %>%
  mutate(
    id_num         = as.integer(id_microrregiao),
    log_pop_18_24  = log(pmax(pop_18_24, 1)),
    share_nonwhite = ifelse(pop_18_24 > 0, pop_nonwhite / pop_18_24, 0)
  ) %>%
  select(id_num, pop_18_24, log_pop_18_24, share_nonwhite,
         share_no_formal_income, avg_hh_income_pc_mw,
         literacy_rate, share_in_school)

# Low-Dose DR formula: 4 covariates (Fix 3)
dr_xformla_low  <- ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate + share_in_school

# High-Dose DR formula: 2 covariates only (avoids singular matrix, Fix 3)
dr_xformla_high <- ~ log_pop_18_24 + avg_hh_income_pc_mw

message(sprintf("   v5 DR covariates: %d microregions | Low spec: 4 vars | High spec: 2 vars",
                nrow(df_covars_dr)))


# ==============================================================================
# SECTION 7: CASE STUDIES (4 auto-selected representative microregions)
# ==============================================================================
message(" [6/12] Selecting Case Study Microregions...")

candidates <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin), g > 0) %>%
  mutate(region_type = if_else(sigla_uf %in% params$central_states, "Central", "Marginal"))

case_micro <- bind_rows(
  candidates %>% filter(region_type == "Central",  dose_bin == "Low")  %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Central",  dose_bin == "High") %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Marginal", dose_bin == "Low")  %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Marginal", dose_bin == "High") %>% slice_sample(n = 1)
) %>%
  select(id_microrregiao, nome_microrregiao, sigla_uf, region_type, dose_bin, g, D_rt_2019)

stopifnot("Need exactly 4 case study units" = nrow(case_micro) == 4)
write_csv(case_micro, "Tables_Final/tab_casestudy_units_v5.csv")
message("   Selected case study microregions:")
print(case_micro)

# ==============================================================================
# SECTION 7B: GEOGRAPHIC DISTRIBUTION HEATMAP
# ==============================================================================
message(" Generating Geographic Heatmap...")

tryCatch({
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("geobr", quietly = TRUE)) {
    stop("sf or geobr package not available")
  }

  micro_sf <- geobr::read_micro_region(year = 2010, simplified = TRUE)

  map_data <- micro_sf %>%
    mutate(code_micro = as.character(code_micro)) %>%
    left_join(
      df_final %>%
        filter(ano == max(params$years)) %>%
        distinct(id_microrregiao, D_rt_2019, dose_bin) %>%
        mutate(id_microrregiao = as.character(id_microrregiao)),
      by = c("code_micro" = "id_microrregiao")
    )

  # Continuous heatmap: scholarship density
  p_heat <- ggplot(map_data) +
    geom_sf(aes(fill = D_rt_2019), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      name      = expression(D[r*","*2019]),
      option    = "inferno",
      direction = -1,
      na.value  = "grey90",
      trans     = "sqrt",
      breaks    = c(0, 10, 50, 100, 200),
      labels    = c("0", "10", "50", "100", "200")
    ) +
    labs(
      title    = "ProUni Scholarship Density by Microregion",
      subtitle = "Cumulative scholarships per 1,000 youth (lagged 4 years), 2019"
    ) +
    theme_void() +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      legend.position  = c(0.15, 0.25),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.4, "cm")
    )
  ggsave("Figures_Final/fig_map_dose_continuous_v5.pdf", p_heat,
         width = 8, height = 9, dpi = 300, bg = "white")

  # Categorical map: dose bins (Low / Mid / High / Baseline)
  map_data <- map_data %>%
    mutate(
      bin_label = case_when(
        dose_bin == "High" ~ "High Dose (> p75)",
        dose_bin == "Mid"  ~ "Mid Dose (p50-p75)",
        dose_bin == "Low"  ~ "Low Dose (p20-p50)",
        is.na(dose_bin) & !is.na(D_rt_2019) ~ "Baseline / Middle",
        TRUE ~ "No data"
      ),
      bin_label = factor(bin_label, levels = c(
        "High Dose (> p75)", "Mid Dose (p50-p75)", "Low Dose (p20-p50)",
        "Baseline / Middle", "No data"
      ))
    )

  p_bins <- ggplot(map_data) +
    geom_sf(aes(fill = bin_label), color = "white", linewidth = 0.05) +
    scale_fill_manual(
      name   = "Dose Bin",
      values = c(
        "High Dose (> p75)"  = "#d73027",
        "Mid Dose (p50-p75)" = "#fdae61",
        "Low Dose (p20-p50)" = "#fee08b",
        "Baseline / Middle"  = "#e0e0e0",
        "No data"            = "white"
      ),
      drop = FALSE
    ) +
    labs(
      title    = "Treatment Assignment by Dose Bin",
      subtitle = "Based on 2019 cross-sectional distribution of lagged scholarship density"
    ) +
    theme_void() +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = c(0.15, 0.25)
    )
  ggsave("Figures_Final/fig_map_dose_bins_v5.pdf", p_bins,
         width = 8, height = 9, dpi = 300, bg = "white")

  message("   Heatmap figures saved.")
}, error = function(e) message(sprintf("   WARNING: Geographic heatmap failed: %s", conditionMessage(e))))



# ==============================================================================
# SECTION 8: PREPARE REGRESSION DATASETS
# ==============================================================================
message(" [7/12] Building Regression Datasets...")

# 8.1 Build the FULL merged dataset (all races, both sexes)
df_reg_all <- df_rais %>%
  inner_join(
    df_geo %>% select(id_municipio, id_microrregiao),
    by = "id_municipio"
  ) %>%
  mutate(
    ano    = as.integer(ano),
    id_num = as.integer(id_microrregiao)
  ) %>%
  inner_join(
    df_final %>% select(id_microrregiao, ano, g, dose_bin, sigla_uf),
    by = c("id_microrregiao", "ano")
  )

# 8.2 Demographic subsets
# RAIS MTE race coding: 1=Indigena, 2=Branca, 4=Preta, 6=Amarela, 8=Parda
# Non-white = codes 1 (Indigena), 4 (Preta), 8 (Parda); White = code 2 (Branca)
df_nonwhite <- df_reg_all %>% filter(raca_cor %in% c("1", "4", "8"))
df_white    <- df_reg_all %>% filter(raca_cor == "2")

# Gender subsets (from non-white, matching main population)
# RAIS sexo: 1=Masculino, 2=Feminino
df_male   <- df_nonwhite %>% filter(sexo == "1")
df_female <- df_nonwhite %>% filter(sexo == "2")

# 8.3 Dose-bin subsets for MAIN RESULTS (non-white, both genders)
data_low  <- df_nonwhite %>% filter(dose_bin == "Low"  | g == 0)
data_mid  <- df_nonwhite %>% filter(dose_bin == "Mid"  | g == 0)   # NEW v5
data_high <- df_nonwhite %>% filter(dose_bin == "High" | g == 0)

# 8.4 Regional heterogeneity subsets (non-white, both genders)
data_cen <- df_nonwhite %>% filter(sigla_uf %in% params$central_states,  !is.na(dose_bin) | g == 0)
data_mar <- df_nonwhite %>% filter(!sigla_uf %in% params$central_states, !is.na(dose_bin) | g == 0)

# 8.4b FIVE OFFICIAL MACRO-REGIONS (IBGE) -- N, NE, CO, SE, S
# Supports the Belíndia-framework claim (Discussion 4.3) that positive Low-Dose
# effects concentrate in Centre-West / North / Northeast while null / negative
# High-Dose effects concentrate in South / Southeast.
uf_region_map <- tibble(
  sigla_uf = c("AC","AM","AP","PA","RO","RR","TO",                 # Norte
               "AL","BA","CE","MA","PB","PE","PI","RN","SE",       # Nordeste
               "DF","GO","MT","MS",                                # Centro-Oeste
               "ES","MG","RJ","SP",                                # Sudeste
               "PR","RS","SC"),                                    # Sul
  macro_region = c(rep("Norte", 7), rep("Nordeste", 9),
                   rep("Centro-Oeste", 4), rep("Sudeste", 4),
                   rep("Sul", 3))
)
df_nonwhite_reg <- df_nonwhite %>% left_join(uf_region_map, by = "sigla_uf")
data_region_N  <- df_nonwhite_reg %>% filter(macro_region == "Norte"        | g == 0)
data_region_NE <- df_nonwhite_reg %>% filter(macro_region == "Nordeste"     | g == 0)
data_region_CO <- df_nonwhite_reg %>% filter(macro_region == "Centro-Oeste" | g == 0)
data_region_SE <- df_nonwhite_reg %>% filter(macro_region == "Sudeste"      | g == 0)
data_region_S  <- df_nonwhite_reg %>% filter(macro_region == "Sul"          | g == 0)

# 8.5 Gender heterogeneity subsets
data_low_male    <- df_male   %>% filter(dose_bin == "Low"  | g == 0)
data_low_female  <- df_female %>% filter(dose_bin == "Low"  | g == 0)
data_high_male   <- df_male   %>% filter(dose_bin == "High" | g == 0)
data_high_female <- df_female %>% filter(dose_bin == "High" | g == 0)

# 8.6 Race heterogeneity subsets
data_low_white     <- df_white    %>% filter(dose_bin == "Low"  | g == 0)
data_low_nonwhite  <- df_nonwhite %>% filter(dose_bin == "Low"  | g == 0)
data_high_white    <- df_white    %>% filter(dose_bin == "High" | g == 0)
data_high_nonwhite <- df_nonwhite %>% filter(dose_bin == "High" | g == 0)

message(sprintf("   Main:    low=%d | mid=%d | high=%d rows",
                nrow(data_low), nrow(data_mid), nrow(data_high)))
message(sprintf("   Region:  cen=%d | mar=%d rows",  nrow(data_cen), nrow(data_mar)))
message(sprintf("   Gender:  low_m=%d | low_f=%d | high_m=%d | high_f=%d",
                nrow(data_low_male), nrow(data_low_female),
                nrow(data_high_male), nrow(data_high_female)))

# ---- DESCRIPTIVE CHECKS: Regression dataset quality --------------------------
message("\n=== REGRESSION DATASET DIAGNOSTICS ===")

# Combined dataset built once and reused for diagnostics + plot below
combined_bin_data <- bind_rows(
  data_low  %>% mutate(bin = "Low"),
  data_mid  %>% mutate(bin = "Mid"),
  data_high %>% mutate(bin = "High")
)

# Sample sizes and microregion counts
reg_diag <- combined_bin_data %>%
  group_by(bin) %>%
  summarise(n_obs = n(), n_regions = n_distinct(id_microrregiao), .groups = "drop")
print(reg_diag)

# Wage distribution by bin
wage_diag <- combined_bin_data %>%
  group_by(bin) %>%
  summarise(
    n_obs     = n(),
    mean_wage = round(mean(log_wage_sm, na.rm = TRUE), 4),
    sd_wage   = round(sd(log_wage_sm,   na.rm = TRUE), 4),
    pct_NA    = round(100 * mean(is.na(log_wage_sm)), 2),
    .groups   = "drop"
  )
message("  Log-wage summary by dose bin:")
print(wage_diag)

# Extreme outlier check
n_extreme <- combined_bin_data %>%
  filter(!is.na(log_wage_sm)) %>%
  summarise(n_neg = sum(log_wage_sm < -3), n_high = sum(log_wage_sm > 5))
if (n_extreme$n_neg > 0 || n_extreme$n_high > 0) {
  message(sprintf("  NOTE: %d obs with log_wage < -3 and %d obs with log_wage > 5 (check outliers)",
                  n_extreme$n_neg, n_extreme$n_high))
} else {
  message("  [OK] No extreme log-wage outliers (<-3 or >5).")
}
message("=== END REGRESSION DIAGNOSTICS ===\n")

# ---- CONSOLE PLOT: Wage distribution by dose bin --------------------------------
p_wage_dist <- combined_bin_data %>%
  filter(!is.na(log_wage_sm)) %>%
  ggplot(aes(x = factor(bin, levels = c("Low", "Mid", "High")),
             y = log_wage_sm, fill = bin)) +
  geom_violin(alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.10, outlier.size = 0.4, color = "grey20",
               outlier.alpha = 0.3) +
  scale_fill_manual(values = DOSE_COLORS) +
  labs(
    title    = "Log-Wage Distribution by Dose Bin",
    subtitle = "Formal sector, ages 22-35, non-white workers  |  RAIS 2005-2019",
    x = "Dose bin", y = "Log wage (min. wages)"
  ) +
  theme(legend.position = "none")
print(p_wage_dist)

# 8.7 Case study wage trajectory plot
cs_ids <- case_micro$id_microrregiao
df_case_plot <- df_nonwhite %>%
  filter(id_microrregiao %in% cs_ids) %>%
  inner_join(
    case_micro %>% select(id_microrregiao, nome_microrregiao, region_type),
    by = "id_microrregiao"
  ) %>%
  mutate(event_time = ano - g) %>%
  group_by(id_microrregiao, nome_microrregiao, dose_bin, region_type, event_time) %>%
  summarise(avg_log_wage = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
            .groups = "drop") %>%
  filter(event_time >= -5, event_time <= 5)

p_cases <- ggplot(df_case_plot, aes(
  x = event_time, y = avg_log_wage,
  color = interaction(dose_bin, region_type, sep = " / "),
  group = id_microrregiao
)) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0,  linetype = "dotted", color = "black") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  facet_wrap(~nome_microrregiao, scales = "free_y") +
  labs(
    title    = "Case Study: Individual Wage Trajectories",
    subtitle = "Event time relative to first treatment year (e = 0)",
    x = "Event time (e)", y = "Avg. log wage (min. wages)", color = "Bin / Region"
  ) +
  scale_x_continuous(breaks = -5:5)
ggsave("Figures_Final/fig_0_casestudy_v5.pdf", p_cases, width = 10, height = 7, dpi = 300)
message("   Case study figure saved.")

# ==============================================================================
# SECTION 8B: COHORT-BASED OUTCOME CONSTRUCTION (Duflo-style)
# ==============================================================================
message(" [NEW] Building cohort-based outcomes...")

df_cohort_reg <- df_rais_cohort %>%
  inner_join(
    df_geo %>% select(id_municipio, id_microrregiao),
    by = "id_municipio"
  ) %>%
  mutate(
    ano    = as.integer(ano),
    id_num = as.integer(id_microrregiao)
  ) %>%
  inner_join(
    df_final %>% select(id_microrregiao, ano, g, dose_bin, sigla_uf),
    by = c("id_microrregiao", "ano")
  )

df_cohort_nw <- df_cohort_reg %>% filter(raca_cor %in% c("1", "4", "8"))

df_cohort_micro <- df_cohort_nw %>%
  group_by(id_num, id_microrregiao, ano, g, dose_bin, age_group) %>%
  summarise(
    avg_log_wage = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
    total_emp    = sum(n_vinculos, na.rm = TRUE),
    .groups = "drop"
  )

df_cohort_wide <- df_cohort_micro %>%
  pivot_wider(
    id_cols     = c(id_num, id_microrregiao, ano, g, dose_bin),
    names_from  = age_group,
    values_from = c(avg_log_wage, total_emp),
    names_sep   = "_"
  ) %>%
  filter(!is.na(avg_log_wage_young), !is.na(avg_log_wage_old)) %>%
  mutate(
    wage_gap        = avg_log_wage_young - avg_log_wage_old,
    emp_ratio_young = total_emp_young / (total_emp_young + total_emp_old)
  )

df_young_only <- df_cohort_micro %>% filter(age_group == "young") %>%
  rename(y_wage = avg_log_wage, y_emp = total_emp)
df_old_only   <- df_cohort_micro %>% filter(age_group == "old") %>%
  rename(y_wage = avg_log_wage, y_emp = total_emp)

message(sprintf("   Cohort panel: %d microregion-years with both cohorts", nrow(df_cohort_wide)))

# ==============================================================================
# SECTION 8C: COHORT REGRESSION DATASETS
# ==============================================================================
gap_low  <- df_cohort_wide %>% filter(dose_bin == "Low"  | g == 0) %>%
  mutate(y = wage_gap) %>% select(id_num, ano, g, y)
gap_high <- df_cohort_wide %>% filter(dose_bin == "High" | g == 0) %>%
  mutate(y = wage_gap) %>% select(id_num, ano, g, y)

young_low  <- df_young_only %>% filter(dose_bin == "Low"  | g == 0) %>%
  mutate(y = y_wage) %>% select(id_num, ano, g, y)
young_high <- df_young_only %>% filter(dose_bin == "High" | g == 0) %>%
  mutate(y = y_wage) %>% select(id_num, ano, g, y)

old_low  <- df_old_only %>% filter(dose_bin == "Low"  | g == 0) %>%
  mutate(y = y_wage) %>% select(id_num, ano, g, y)
old_high <- df_old_only %>% filter(dose_bin == "High" | g == 0) %>%
  mutate(y = y_wage) %>% select(id_num, ano, g, y)

message(sprintf("   Gap:   low=%d | high=%d", nrow(gap_low), nrow(gap_high)))
message(sprintf("   Young: low=%d | high=%d", nrow(young_low), nrow(young_high)))
message(sprintf("   Old:   low=%d | high=%d", nrow(old_low), nrow(old_high)))



# ==============================================================================
# SECTION 8D: EDUCATION OUTCOMES (First Stage for LATE)
# ==============================================================================
message(" SECTION 8D: Building education outcomes...")

df_educ_reg <- df_rais_educ %>%
  inner_join(df_geo %>% select(id_municipio, id_microrregiao), by = "id_municipio") %>%
  mutate(ano = as.integer(ano), id_num = as.integer(id_microrregiao)) %>%
  inner_join(
    df_final %>% select(id_microrregiao, ano, g, dose_bin),
    by = c("id_microrregiao", "ano")
  ) %>%
  filter(raca_cor %in% c("1", "4", "8"))

df_educ_micro <- df_educ_reg %>%
  group_by(id_num, id_microrregiao, ano, g, dose_bin, age_group) %>%
  summarise(
    total_workers = sum(n_vinculos, na.rm = TRUE),
    total_college = sum(n_college, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(share_college = total_college / total_workers)

df_educ_wide <- df_educ_micro %>%
  pivot_wider(
    id_cols    = c(id_num, id_microrregiao, ano, g, dose_bin),
    names_from = age_group,
    values_from = c(share_college, total_workers),
    names_sep = "_"
  ) %>%
  filter(!is.na(share_college_young), !is.na(share_college_old)) %>%
  mutate(college_gap = share_college_young - share_college_old)

message(sprintf("   Education panel: %d microregion-years", nrow(df_educ_wide)))

# ==============================================================================
# SECTION 8E: EXTENDED REGRESSION DATASETS
# ==============================================================================
message(" SECTION 8E: Building extended regression datasets...")

make_subset <- function(df, outcome_col, bin) {
  df %>%
    filter(dose_bin == bin | g == 0) %>%
    mutate(y = .data[[outcome_col]]) %>%
    select(id_num, ano, g, y)
}

college_gap_low  <- make_subset(df_educ_wide, "college_gap", "Low")
college_gap_high <- make_subset(df_educ_wide, "college_gap", "High")

# Gender heterogeneity on wage gap
for (sex_val in c("1", "2")) {
  sex_label <- ifelse(sex_val == "1", "male", "female")
  df_sex <- df_cohort_reg %>%
    filter(sexo == sex_val, raca_cor %in% c("1", "4", "8")) %>%
    group_by(id_num, id_microrregiao, ano, g, dose_bin, age_group) %>%
    summarise(avg_log_wage = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(id_cols = c(id_num, id_microrregiao, ano, g, dose_bin),
                names_from = age_group, values_from = avg_log_wage, names_sep = "_") %>%
    filter(!is.na(young), !is.na(old)) %>%
    mutate(wage_gap = young - old)

  assign(paste0("gap_", sex_label, "_low"),
         df_sex %>% filter(dose_bin == "Low" | g == 0) %>%
           mutate(y = wage_gap) %>% select(id_num, ano, g, y))
  assign(paste0("gap_", sex_label, "_high"),
         df_sex %>% filter(dose_bin == "High" | g == 0) %>%
           mutate(y = wage_gap) %>% select(id_num, ano, g, y))
}

# Race heterogeneity on wage gap -- White subsample
df_cohort_white <- df_cohort_reg %>%
  filter(raca_cor %in% c("2")) %>%
  group_by(id_num, id_microrregiao, ano, g, dose_bin, age_group) %>%
  summarise(avg_log_wage = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(id_cols = c(id_num, id_microrregiao, ano, g, dose_bin),
              names_from = age_group, values_from = avg_log_wage, names_sep = "_") %>%
  filter(!is.na(young), !is.na(old)) %>%
  mutate(wage_gap = young - old)

gap_white_low  <- df_cohort_white %>% filter(dose_bin == "Low" | g == 0) %>%
  mutate(y = wage_gap) %>% select(id_num, ano, g, y)
gap_white_high <- df_cohort_white %>% filter(dose_bin == "High" | g == 0) %>%
  mutate(y = wage_gap) %>% select(id_num, ano, g, y)

gap_nonwhite_low  <- gap_low
gap_nonwhite_high <- gap_high

message("   Extended regression datasets built.")

# ==============================================================================
# SECTION 8F: BIRTH-COHORT OUTCOMES (Duflo-strict)
# ==============================================================================
message(" SECTION 8F: Building birth-cohort outcomes...")

bc_gap_low <- NULL; bc_gap_high <- NULL
bc_exposed_low <- NULL; bc_control_low <- NULL; df_bc_wide <- NULL

if (!is.null(df_rais_bc) && nrow(df_rais_bc) > 0) {
  df_bc_reg <- df_rais_bc %>%
    filter(!is.na(birth_cohort)) %>%
    inner_join(df_geo %>% select(id_municipio, id_microrregiao), by = "id_municipio") %>%
    mutate(ano = as.integer(ano), id_num = as.integer(id_microrregiao)) %>%
    inner_join(
      df_final %>% select(id_microrregiao, ano, g, dose_bin, D_rt, D_rt_2019),
      by = c("id_microrregiao", "ano")
    )

  df_bc_nw <- df_bc_reg %>% filter(raca_cor %in% c("1", "4", "8"))

  df_bc_micro <- df_bc_nw %>%
    group_by(id_num, id_microrregiao, ano, g, dose_bin, D_rt, D_rt_2019, birth_cohort) %>%
    summarise(
      avg_log_wage  = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
      total_emp     = sum(n_vinculos, na.rm = TRUE),
      total_college = sum(n_college, na.rm = TRUE),
      share_college = sum(n_college, na.rm = TRUE) / sum(n_vinculos, na.rm = TRUE),
      .groups = "drop"
    )

  df_bc_wide <- df_bc_micro %>%
    pivot_wider(
      id_cols     = c(id_num, id_microrregiao, ano, g, dose_bin, D_rt, D_rt_2019),
      names_from  = birth_cohort,
      values_from = c(avg_log_wage, total_emp, share_college),
      names_sep   = "_"
    ) %>%
    filter(!is.na(avg_log_wage_exposed), !is.na(avg_log_wage_control)) %>%
    mutate(
      wage_gap_bc    = avg_log_wage_exposed - avg_log_wage_control,
      college_gap_bc = share_college_exposed - share_college_control
    ) %>%
    filter(ano >= 2009)

  message(sprintf("   Birth-cohort panel (ano>=2009): %d microregion-years", nrow(df_bc_wide)))

  bc_gap_low  <- df_bc_wide %>% filter(dose_bin == "Low" | g == 0) %>%
    mutate(y = wage_gap_bc) %>% select(id_num, ano, g, y)
  bc_gap_high <- df_bc_wide %>% filter(dose_bin == "High" | g == 0) %>%
    mutate(y = wage_gap_bc) %>% select(id_num, ano, g, y)

  df_bc_exposed <- df_bc_micro %>% filter(birth_cohort == "exposed")
  df_bc_control <- df_bc_micro %>% filter(birth_cohort == "control")

  bc_exposed_low <- df_bc_exposed %>% filter(dose_bin == "Low" | g == 0) %>%
    mutate(y = avg_log_wage) %>% select(id_num, ano, g, y)
  bc_control_low <- df_bc_control %>% filter(dose_bin == "Low" | g == 0) %>%
    mutate(y = avg_log_wage) %>% select(id_num, ano, g, y)

  message("   Birth-cohort regression datasets built.")
} else {
  message("   SKIPPED: Birth-cohort data not available.")
}

# ==============================================================================
# v5 FIX 4: TRIPLE INTERACTION -- Young workers in Low-Dose + below-median income
# ==============================================================================
message(" [v5 FIX 4] Building triple-interaction dataset (young x low x poor)...")

res_triple_young_poor <- NULL
tryCatch({
  income_median <- median(df_covars_dr$avg_hh_income_pc_mw, na.rm = TRUE)
  poor_ids <- df_covars_dr %>%
    filter(avg_hh_income_pc_mw <= income_median) %>%
    pull(id_num)

  triple_data <- young_low %>%
    filter(id_num %in% poor_ids)

  message(sprintf("   Triple interaction: %d obs in %d poor Low-Dose microregions",
                  nrow(triple_data), n_distinct(triple_data$id_num)))

  res_triple_young_poor <- run_cs_preagg(triple_data, "Young/Low/Poor")
}, error = function(e) {
  message(sprintf("   WARNING: Triple interaction failed: %s", safe_err_msg(e)))
})



# Defensive NULL initialization for all result objects -- ensures Section 9L
# can build the all_results_v5 list even if earlier sections error out.
# (Variables are also assigned in their respective sections; these are guards only.)
for (.v in c(
  "res_low", "res_mid", "res_high", "res_cen", "res_mar",
  "res_low_emp", "res_high_emp", "res_cen_emp", "res_mar_emp",
  "res_male_low", "res_female_low", "res_male_high", "res_female_high",
  "res_white_low", "res_nonwhite_low", "res_white_high", "res_nonwhite_high",
  "res_low_lag3", "res_high_lag3",
  "res_low_dr", "res_high_dr", "res_low_dr_censo", "res_high_dr_censo",
  "res_low_buffer", "res_high_buffer",
  "res_placebo_low", "res_placebo_high",
  "res_low_25_75", "res_high_25_75", "res_low_33_67", "res_high_33_67",
  "res_gap_low", "res_gap_high", "res_young_low", "res_young_high",
  "res_old_low", "res_old_high",
  "res_college_low", "res_college_high",
  "res_gap_male_low", "res_gap_male_high",
  "res_gap_female_low", "res_gap_female_high",
  "res_gap_white_low", "res_gap_white_high",
  "res_gap_nonwhite_low", "res_gap_nonwhite_high",
  "res_bc_gap_low", "res_bc_gap_high",
  "res_bc_exposed_low", "res_bc_control_low",
  "res_triple_young_poor",
  "res_region_N", "res_region_NE", "res_region_CO",
  "res_region_SE", "res_region_S",
  "res_low_hq", "res_low_noq",
  "res_low_dr_policy", "res_high_dr_policy",
  "twfe_fies", "bacon_decomp", "tab_spt_balance"
)) {
  if (!exists(.v)) assign(.v, NULL)
}
rm(.v)

# ==============================================================================
# SECTION 9A: MAIN ESTIMATION -- WAGES (non-white, 4 models)
# ==============================================================================
message(" [8/12] Running Main Wage Models (Callaway & Sant'Anna)...")

message("   -> Low Dose / Wages...")
res_low  <- run_cs(data_low,  "log_wage_sm", label = "Low/Wages")

message("   -> Mid Dose / Wages (NEW v5)...")
res_mid  <- tryCatch(run_cs(data_mid, "log_wage_sm", label = "Mid/Wages"),
                     error = function(e) { message(sprintf("   WARNING: Mid/Wages failed: %s", safe_err_msg(e))); NULL })

message("   -> High Dose / Wages...")
res_high <- run_cs(data_high, "log_wage_sm", label = "High/Wages")

message("   -> Central / Wages...")
res_cen  <- run_cs(data_cen,  "log_wage_sm", label = "Central/Wages")

message("   -> Marginal / Wages...")
res_mar  <- run_cs(data_mar,  "log_wage_sm", label = "Marginal/Wages")

# 9A.b  Five IBGE macro-regions (Belíndia framework test)
message("   -> IBGE 5 macro-regions (Norte / Nordeste / Centro-Oeste / Sudeste / Sul)...")
res_region_N  <- tryCatch(run_cs(data_region_N,  "log_wage_sm", label = "Region/Norte"),
                          error = function(e) { message(sprintf("     Norte failed: %s", safe_err_msg(e))); NULL })
res_region_NE <- tryCatch(run_cs(data_region_NE, "log_wage_sm", label = "Region/Nordeste"),
                          error = function(e) { message(sprintf("     Nordeste failed: %s", safe_err_msg(e))); NULL })
res_region_CO <- tryCatch(run_cs(data_region_CO, "log_wage_sm", label = "Region/CentroOeste"),
                          error = function(e) { message(sprintf("     CentroOeste failed: %s", safe_err_msg(e))); NULL })
res_region_SE <- tryCatch(run_cs(data_region_SE, "log_wage_sm", label = "Region/Sudeste"),
                          error = function(e) { message(sprintf("     Sudeste failed: %s", safe_err_msg(e))); NULL })
res_region_S  <- tryCatch(run_cs(data_region_S,  "log_wage_sm", label = "Region/Sul"),
                          error = function(e) { message(sprintf("     Sul failed: %s", safe_err_msg(e))); NULL })

# ---- CONSOLE: Main CS estimates quick summary -----------------------------------
cat("\n=== MAIN CS ESTIMATES SUMMARY ===\n")
for (r in Filter(Negate(is.null),
                 list(res_low, res_mid, res_high, res_cen, res_mar,
                      res_region_N, res_region_NE, res_region_CO,
                      res_region_SE, res_region_S))) {
  row <- build_att_row(r)
  cat(sprintf("  %-22s  ATT=%+.4f  SE=%.4f  p=%.3f  pre_p=%.3f\n",
              row$Model, row$ATT, row$SE, row$p_val, row$Pretrend_p))
}
cat("=================================\n\n")

# ==============================================================================
# SECTION 9B: EMPLOYMENT MARGIN -- n_vinculos (non-white, 4 models)
# ==============================================================================
message("   -> Low Dose / Employment...")
res_low_emp  <- run_cs(data_low,  "n_vinculos", label = "Low/Emp",  agg = "logsum")

message("   -> High Dose / Employment...")
res_high_emp <- run_cs(data_high, "n_vinculos", label = "High/Emp", agg = "logsum")

message("   -> Central / Employment...")
res_cen_emp  <- run_cs(data_cen,  "n_vinculos", label = "Central/Emp", agg = "logsum")

message("   -> Marginal / Employment...")
res_mar_emp  <- run_cs(data_mar,  "n_vinculos", label = "Marginal/Emp", agg = "logsum")

# ==============================================================================
# SECTION 9C: GENDER HETEROGENEITY -- log wages (4 models)
# ==============================================================================
message(" [9/12] Running Gender Heterogeneity Models...")

message("   -> Male / Low Dose...")
res_male_low    <- run_cs(data_low_male,    "log_wage_sm", label = "Male/Low")

message("   -> Female / Low Dose...")
res_female_low  <- run_cs(data_low_female,  "log_wage_sm", label = "Female/Low")

message("   -> Male / High Dose...")
res_male_high   <- run_cs(data_high_male,   "log_wage_sm", label = "Male/High")

message("   -> Female / High Dose...")
res_female_high <- run_cs(data_high_female, "log_wage_sm", label = "Female/High")

# ==============================================================================
# SECTION 9D: RACE HETEROGENEITY -- log wages (4 models)
# ==============================================================================
message("   Running Race Heterogeneity Models...")

message("   -> White / Low Dose...")
res_white_low     <- run_cs(data_low_white,     "log_wage_sm", label = "White/Low")

message("   -> Non-white / Low Dose...")
res_nonwhite_low  <- run_cs(data_low_nonwhite,  "log_wage_sm", label = "Nonwhite/Low")

message("   -> White / High Dose...")
res_white_high    <- run_cs(data_high_white,    "log_wage_sm", label = "White/High")

message("   -> Non-white / High Dose...")
res_nonwhite_high <- run_cs(data_high_nonwhite, "log_wage_sm", label = "Nonwhite/High")

# ==============================================================================
# SECTION 9E: ROBUSTNESS -- ALTERNATIVE LAG (lag = 3 years, 2 models)
# ==============================================================================
message(" [10/12] Running Robustness: Alternative Lag (3 years)...")

lag_alt <- 3L
df_panel_lag3 <- df_panel %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(D_rt_alt = (dplyr::lag(cum_schol, lag_alt, default = 0) / pop) * 1000) %>%
  ungroup()

df_2019_lag3 <- df_panel_lag3 %>% filter(ano == max(params$years))
tau_20_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_low,  na.rm = TRUE)
tau_50_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_mid,  na.rm = TRUE)
tau_75_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_high, na.rm = TRUE)

df_final_lag3 <- df_panel_lag3 %>%
  left_join(
    df_2019_lag3 %>% select(id_microrregiao, D_rt_2019_alt = D_rt_alt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    dose_bin_alt = case_when(
      D_rt_2019_alt >  tau_20_alt & D_rt_2019_alt <= tau_50_alt ~ "Low",
      D_rt_2019_alt >  tau_75_alt                                ~ "High",
      TRUE                                                        ~ NA_character_
    ),
    g_alt = {
      crossing <- ano[D_rt_alt > tau_20_alt]
      if (length(crossing) == 0) 0L else as.integer(min(crossing))
    }
  ) %>%
  ungroup() %>%
  mutate(g_alt = ifelse(is.infinite(g_alt), 0L, g_alt))

df_rob <- df_nonwhite %>%
  select(-g, -dose_bin) %>%
  inner_join(
    df_final_lag3 %>% select(id_microrregiao, ano, g_alt, dose_bin_alt),
    by = c("id_microrregiao", "ano")
  ) %>%
  rename(g = g_alt, dose_bin = dose_bin_alt)

data_low_lag3  <- df_rob %>% filter(dose_bin == "Low"  | g == 0)
data_high_lag3 <- df_rob %>% filter(dose_bin == "High" | g == 0)

message("   -> Low Dose / Lag=3...")
res_low_lag3  <- run_cs(data_low_lag3,  "log_wage_sm", label = "Low/Lag3")

message("   -> High Dose / Lag=3...")
res_high_lag3 <- run_cs(data_high_lag3, "log_wage_sm", label = "High/Lag3")

message(" -> All 18 baseline + 2 lag-robustness models completed.")



# ==============================================================================
# SECTION 9E2: COHORT-BASED ESTIMATION (Duflo-style)
# ==============================================================================
message(" [NEW] Running Cohort-Based Models...")

res_gap_low   <- NULL; res_gap_high   <- NULL
res_young_low <- NULL; res_young_high <- NULL
res_old_low   <- NULL; res_old_high   <- NULL

message("   -> Wage Gap (young - old) / Low Dose...")
res_gap_low  <- tryCatch(run_cs_preagg(gap_low,  "Gap/Low"),
                         error = function(e) { message(sprintf("   WARNING: Gap/Low failed: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap (young - old) / High Dose...")
res_gap_high <- tryCatch(run_cs_preagg(gap_high, "Gap/High"),
                         error = function(e) { message(sprintf("   WARNING: Gap/High failed: %s", safe_err_msg(e))); NULL })

message("   -> Young (22-28) / Low Dose...")
res_young_low  <- tryCatch(run_cs_preagg(young_low,  "Young/Low"),
                           error = function(e) { message(sprintf("   WARNING: Young/Low failed: %s", safe_err_msg(e))); NULL })

message("   -> Young (22-28) / High Dose...")
res_young_high <- tryCatch(run_cs_preagg(young_high, "Young/High"),
                           error = function(e) { message(sprintf("   WARNING: Young/High failed: %s", safe_err_msg(e))); NULL })

message("   -> Old (29-35) / Low Dose [PLACEBO]...")
res_old_low  <- tryCatch(run_cs_preagg(old_low,  "Old/Low"),
                         error = function(e) { message(sprintf("   WARNING: Old/Low failed: %s", safe_err_msg(e))); NULL })

message("   -> Old (29-35) / High Dose [PLACEBO]...")
res_old_high <- tryCatch(run_cs_preagg(old_high, "Old/High"),
                         error = function(e) { message(sprintf("   WARNING: Old/High failed: %s", safe_err_msg(e))); NULL })

message(" -> Cohort-based estimation completed.")

# ==============================================================================
# SECTION 9E3: EXTENDED COHORT MODELS (gender, race, college gap)
# ==============================================================================
message(" [NEW] Running Extended Cohort Models...")

res_college_low <- NULL; res_college_high <- NULL
res_gap_male_low <- NULL; res_gap_male_high <- NULL
res_gap_female_low <- NULL; res_gap_female_high <- NULL
res_gap_white_low <- NULL; res_gap_white_high <- NULL
res_gap_nonwhite_low <- NULL; res_gap_nonwhite_high <- NULL

message("   -> College Gap / Low Dose...")
res_college_low  <- tryCatch(run_cs_preagg(college_gap_low,  "CollegeGap/Low"),
                             error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> College Gap / High Dose...")
res_college_high <- tryCatch(run_cs_preagg(college_gap_high, "CollegeGap/High"),
                             error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Male / Low Dose...")
res_gap_male_low    <- tryCatch(run_cs_preagg(gap_male_low,    "WageGap/Male/Low"),
                                error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Male / High Dose...")
res_gap_male_high   <- tryCatch(run_cs_preagg(gap_male_high,   "WageGap/Male/High"),
                                error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Female / Low Dose...")
res_gap_female_low  <- tryCatch(run_cs_preagg(gap_female_low,  "WageGap/Female/Low"),
                                error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Female / High Dose...")
res_gap_female_high <- tryCatch(run_cs_preagg(gap_female_high, "WageGap/Female/High"),
                                error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / White / Low Dose...")
res_gap_white_low    <- tryCatch(run_cs_preagg(gap_white_low,    "WageGap/White/Low"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / White / High Dose...")
res_gap_white_high   <- tryCatch(run_cs_preagg(gap_white_high,   "WageGap/White/High"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Nonwhite / Low Dose...")
res_gap_nonwhite_low <- tryCatch(run_cs_preagg(gap_nonwhite_low, "WageGap/Nonwhite/Low"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message("   -> Wage Gap / Nonwhite / High Dose...")
res_gap_nonwhite_high<- tryCatch(run_cs_preagg(gap_nonwhite_high,"WageGap/Nonwhite/High"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

message(" -> Extended cohort models completed.")

# ==============================================================================
# SECTION 9F: CONDITIONAL (DOUBLY-ROBUST) MODELS -- v5 Fix 3
# Enriched Low spec (4 covars); reduced High spec (2 covars, avoids singular matrix)
# ==============================================================================
message(" [NEW] Running Conditional (Doubly-Robust) Wage Models...")

res_low_dr  <- tryCatch({
  message("   -> Low Dose / DR Wages (4 covars -- Sant'Anna-Zhao DR)...")
  run_cs(data_low, "log_wage_sm", label = "Low/DR",
         xformla = dr_xformla_low, covars_df = df_covars_dr,
         est_method = "dr")
}, error = function(e) { message(sprintf("   WARNING: Low/DR failed: %s", safe_err_msg(e))); NULL })

res_high_dr <- tryCatch({
  message("   -> High Dose / DR Wages (2 covars -- Sant'Anna-Zhao DR)...")
  run_cs(data_high, "log_wage_sm", label = "High/DR",
         xformla = dr_xformla_high, covars_df = df_covars_dr,
         est_method = "dr")
}, error = function(e) { message(sprintf("   WARNING: High/DR failed: %s", safe_err_msg(e))); NULL })

# Flag divergence between unconditional and conditional estimates
dr_divergence_flag <- FALSE
if (!is.null(res_low_dr) && !is.null(res_high_dr)) {
  for (pair in list(
    list(unc = res_low,  cond = res_low_dr,  lbl = "Low"),
    list(unc = res_high, cond = res_high_dr, lbl = "High")
  )) {
    diff_att <- abs(pair$unc$simple$overall.att - pair$cond$simple$overall.att)
    se_pool  <- sqrt(pair$unc$simple$overall.se^2 + pair$cond$simple$overall.se^2)
    if (se_pool > 0 && diff_att / se_pool > 2) {
      message(sprintf("   FLAG: %s unconditional vs DR ATTs diverge by %.1f pooled SEs!",
                      pair$lbl, diff_att / se_pool))
      dr_divergence_flag <- TRUE
    }
  }
  if (!dr_divergence_flag) message("   OK: Unconditional and DR estimates are broadly consistent.")
} else {
  message("   WARNING: DR comparison skipped (one or both models failed).")
}



# ==============================================================================
# SECTION 9G: HonestDiD SENSITIVITY ANALYSIS -- v5 Fix 2 (run_honest_v5)
# Runs only on Low/Wages and Low/DR (avoid singular matrix in High)
# ==============================================================================
message(" [NEW] Running HonestDiD Sensitivity Analysis (v5 Fix 2)...")

honest_results <- list()
if (requireNamespace("HonestDiD", quietly = TRUE)) {
  for (spec in list(
    list(es = if (!is.null(res_low))    res_low$es    else NULL, lbl = "Low/Wages"),
    list(es = if (!is.null(res_low_dr)) res_low_dr$es else NULL, lbl = "Low/DR")
  )) {
    if (is.null(spec$es)) {
      message(sprintf("   SKIP HonestDiD for %s: model result is NULL", spec$lbl))
      next
    }
    honest_results[[spec$lbl]] <- run_honest_v5(spec$es, spec$lbl)
  }

  # Save sensitivity plots
  for (nm in names(honest_results)) {
    hr <- honest_results[[nm]]
    if (!is.null(hr$results) && nrow(hr$results) > 0) {
      tryCatch({
        p_honest <- ggplot(hr$results, aes(x = Mbar, y = estimate)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.20, fill = "#0072B2") +
          geom_line(color = "#0072B2", linewidth = 1) +
          geom_point(color = "#0072B2", size = 2) +
          labs(
            title    = sprintf("HonestDiD Sensitivity: %s", hr$label),
            subtitle = sprintf("Breakdown Mbar = %.2f | VCov: %s",
                               hr$breakdown,
                               if (isTRUE(hr$used_full_vcov)) "full matrix" else "diagonal fallback"),
            x = expression(bar(M) ~ "(relative magnitudes)"),
            y = "ATT (robust CI)"
          )
        fname <- sprintf("Figures_Final/fig_honestdid_%s_v5.pdf",
                         gsub("[/ ]", "_", tolower(hr$label)))
        ggsave(fname, p_honest, width = 7, height = 5, dpi = 300)
        message(sprintf("   Saved: %s", fname))
      }, error = function(e) {
        message(sprintf("   WARNING: HonestDiD plot failed for %s: %s", nm, safe_err_msg(e)))
      })
    }
  }
} else {
  message("   SKIPPED: HonestDiD package not available.")
}

# ==============================================================================
# SECTION 9H: BUFFER-ZONE EXCLUSION TEST FOR SUTVA
# ==============================================================================
message(" [NEW] Running Buffer-Zone Exclusion Test (SUTVA)...")

res_low_buffer  <- NULL
res_high_buffer <- NULL
border_ids      <- integer(0)

if (requireNamespace("sf", quietly = TRUE) && requireNamespace("geobr", quietly = TRUE)) {
  tryCatch({
    micro_sf_buf <- geobr::read_micro_region(year = 2010, simplified = TRUE)
    neighbors    <- sf::st_touches(micro_sf_buf)
    high_dose_ids <- unique(df_final$id_microrregiao[df_final$dose_bin == "High"])

    for (i in seq_along(micro_sf_buf$code_micro)) {
      code_i <- as.character(micro_sf_buf$code_micro[i])
      if (code_i %in% high_dose_ids) next
      neighbor_codes <- as.character(micro_sf_buf$code_micro[neighbors[[i]]])
      if (any(neighbor_codes %in% high_dose_ids)) border_ids <- c(border_ids, code_i)
    }
    border_ids <- unique(border_ids)
    message(sprintf("   %d comparison microregions border High-Dose units (excluded).",
                    length(border_ids)))

    data_low_nobuffer  <- data_low  %>% filter(!(id_microrregiao %in% border_ids))
    data_high_nobuffer <- data_high %>% filter(!(id_microrregiao %in% border_ids))

    message("   -> Low Dose / Buffer exclusion...")
    res_low_buffer  <- run_cs(data_low_nobuffer,  "log_wage_sm", label = "Low/Buffer")

    message("   -> High Dose / Buffer exclusion...")
    res_high_buffer <- run_cs(data_high_nobuffer, "log_wage_sm", label = "High/Buffer")

  }, error = function(e) message(sprintf("   WARNING: Buffer-zone test failed: %s", safe_err_msg(e))))
} else {
  message("   SKIPPED: sf/geobr packages not available.")
}

# ==============================================================================
# SECTION 9I: PLACEBO COHORT TEST (ages 36-50) -- v5 Fix 6: ACTIVATED
# ==============================================================================
message(" [v5 Fix 6] Placebo Cohort Test (ages 36-50) -- ACTIVATED")

res_placebo_low  <- NULL
res_placebo_high <- NULL

tryCatch({
  placebo_list <- list()
  for (y in params$years) {
    q_placebo <- paste0("
      SELECT
        ano,
        id_municipio,
        raca_cor,
        sexo,
        COUNT(*)                                             AS n_vinculos,
        AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1)))  AS log_wage_sm
      FROM `basedosdados.br_me_rais.microdados_vinculos`
      WHERE
        ano = ", y, "
        AND idade BETWEEN 36 AND 50
        AND valor_remuneracao_media_sm > 0
        AND tipo_vinculo = '10'
        AND vinculo_ativo_3112 = '1'
      GROUP BY ano, id_municipio, raca_cor, sexo
    ")
    tryCatch({
      placebo_list[[as.character(y)]] <- basedosdados::read_sql(q_placebo)
      message(sprintf("     Placebo year %d OK", y))
    }, error = function(e) message(sprintf("     Error placebo year %d: %s", y, safe_err_msg(e))))
    gc()
  }

  df_rais_placebo <- bind_rows(placebo_list)
  rm(placebo_list); gc()

  df_placebo_reg <- df_rais_placebo %>%
    filter(raca_cor %in% c("1", "4", "8")) %>%
    inner_join(df_geo %>% select(id_municipio, id_microrregiao), by = "id_municipio") %>%
    mutate(ano = as.integer(ano), id_num = as.integer(id_microrregiao)) %>%
    inner_join(df_final %>% select(id_microrregiao, ano, g, dose_bin),
               by = c("id_microrregiao", "ano"))

  data_placebo_low  <- df_placebo_reg %>% filter(dose_bin == "Low"  | g == 0)
  data_placebo_high <- df_placebo_reg %>% filter(dose_bin == "High" | g == 0)

  message("   -> Placebo Low Dose / Wages...")
  res_placebo_low  <- run_cs(data_placebo_low,  "log_wage_sm", label = "Placebo/Low")

  message("   -> Placebo High Dose / Wages...")
  res_placebo_high <- run_cs(data_placebo_high, "log_wage_sm", label = "Placebo/High")

}, error = function(e) {
  message(sprintf("   WARNING: Placebo cohort test failed: %s", safe_err_msg(e)))
})

# ==============================================================================
# SECTION 9J: TWFE BENCHMARK
# ==============================================================================
message(" [NEW] Running TWFE Benchmark...")

twfe_fit <- NULL
twfe_emp <- NULL
if (requireNamespace("fixest", quietly = TRUE)) {
  tryCatch({
    twfe_data <- df_nonwhite %>%
      group_by(id_num, ano, g) %>%
      summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(treated = as.integer(g > 0 & ano >= g))

    twfe_fit <- fixest::feols(y ~ treated | id_num + ano, data = twfe_data, cluster = ~id_num)
    message(sprintf("   TWFE Wages:  coef=%.4f (SE=%.4f)",
                    coef(twfe_fit)["treated"], fixest::se(twfe_fit)["treated"]))

    twfe_emp_data <- df_nonwhite %>%
      group_by(id_num, ano, g) %>%
      summarise(y = log(sum(n_vinculos, na.rm = TRUE)), .groups = "drop") %>%
      mutate(treated = as.integer(g > 0 & ano >= g))

    twfe_emp <- fixest::feols(y ~ treated | id_num + ano, data = twfe_emp_data, cluster = ~id_num)
    message(sprintf("   TWFE Employ: coef=%.4f (SE=%.4f)",
                    coef(twfe_emp)["treated"], fixest::se(twfe_emp)["treated"]))

    twfe_summary <- tibble(
      Model = c("TWFE/Wages", "TWFE/Employment"),
      Coef  = c(coef(twfe_fit)["treated"], coef(twfe_emp)["treated"]),
      SE    = c(fixest::se(twfe_fit)["treated"], fixest::se(twfe_emp)["treated"]),
      N     = c(stats::nobs(twfe_fit), stats::nobs(twfe_emp))
    )
    write_csv(twfe_summary, "Tables_Final/tab_twfe_benchmark_v5.csv")
    message("   TWFE benchmark saved.")
  }, error = function(e) message(sprintf("   WARNING: TWFE failed: %s", safe_err_msg(e))))
} else {
  message("   SKIPPED: fixest package not available.")
}

# ==============================================================================
# SECTION 9K: ALTERNATIVE DISCRETIZATION CUTOFFS
# ==============================================================================
message(" [NEW] Running Alternative Cutoff Robustness...")

run_alt_cutoffs <- function(q_lo, q_hi, suffix) {
  tau_lo  <- quantile(df_2019$D_rt, q_lo, na.rm = TRUE)
  tau_mid <- quantile(df_2019$D_rt, mean(c(q_lo, q_hi)), na.rm = TRUE)
  tau_hi  <- quantile(df_2019$D_rt, q_hi, na.rm = TRUE)
  message(sprintf("   Cutoffs %s: tau_lo=%.4f | tau_mid=%.4f | tau_hi=%.4f",
                  suffix, tau_lo, tau_mid, tau_hi))

  df_alt <- df_final %>%
    group_by(id_microrregiao) %>%
    mutate(
      dose_bin_alt = case_when(
        D_rt_2019 >  tau_lo & D_rt_2019 <= tau_mid ~ "Low",
        D_rt_2019 >  tau_hi                          ~ "High",
        TRUE                                          ~ NA_character_
      ),
      g_alt = {
        crossing <- ano[D_rt > tau_lo]
        if (length(crossing) == 0) 0L else as.integer(min(crossing))
      }
    ) %>%
    ungroup() %>%
    mutate(g_alt = ifelse(is.infinite(g_alt), 0L, g_alt))

  data_low_alt <- df_nonwhite %>%
    select(-dose_bin, -g) %>%
    inner_join(df_alt %>% select(id_microrregiao, ano, dose_bin_alt, g_alt),
               by = c("id_microrregiao", "ano")) %>%
    rename(dose_bin = dose_bin_alt, g = g_alt) %>%
    filter(dose_bin == "Low" | g == 0)

  data_high_alt <- df_nonwhite %>%
    select(-dose_bin, -g) %>%
    inner_join(df_alt %>% select(id_microrregiao, ano, dose_bin_alt, g_alt),
               by = c("id_microrregiao", "ano")) %>%
    rename(dose_bin = dose_bin_alt, g = g_alt) %>%
    filter(dose_bin == "High" | g == 0)

  res_lo <- tryCatch(run_cs(data_low_alt,  "log_wage_sm", label = paste0("Low/",  suffix)),
                     error = function(e) { message(sprintf("   WARNING: %s Low failed: %s", suffix, safe_err_msg(e))); NULL })
  res_hi <- tryCatch(run_cs(data_high_alt, "log_wage_sm", label = paste0("High/", suffix)),
                     error = function(e) { message(sprintf("   WARNING: %s High failed: %s", suffix, safe_err_msg(e))); NULL })
  list(low = res_lo, high = res_hi)
}

message("   -> tau_25 / tau_75 cutoffs...")
alt_25_75 <- run_alt_cutoffs(0.25, 0.75, "q25_75")
res_low_25_75  <- alt_25_75$low
res_high_25_75 <- alt_25_75$high

message("   -> tau_33 / tau_67 cutoffs...")
alt_33_67 <- run_alt_cutoffs(1/3, 2/3, "q33_67")
res_low_33_67  <- alt_33_67$low
res_high_33_67 <- alt_33_67$high

message(" -> All robustness and sensitivity sections completed.")



# SECTION 9M removed in final draft: weak first stage precludes
#  Wald/LATE inference; reported as reduced-form throughout.
# ==============================================================================
# SECTION 9M: FIRST-STAGE DIAGNOSTICS  -- DISABLED IN FINAL DRAFT
# ------------------------------------------------------------------------------
# message(" SECTION 9M: First-stage diagnostics (LATE not reported; F documented)...")
#
# fs_rows <- list()
# for (r in list(list(fs = res_college_low,  rf = res_low,  lbl = "Low"),
#                list(fs = res_college_high, rf = res_high, lbl = "High"))) {
#   if (is.null(r$fs) || is.null(r$rf)) next
#   beta_fs <- r$fs$simple$overall.att
#   se_fs   <- r$fs$simple$overall.se
#   beta_rf <- r$rf$simple$overall.att
#   se_rf   <- r$rf$simple$overall.se
#   F_fs    <- if (!is.na(se_fs) && se_fs > 0) (beta_fs / se_fs)^2 else NA_real_
#   ar_status <- if (!is.na(F_fs) && F_fs < 1)
#                  "unbounded" else
#                  "computed"
#   fs_rows[[r$lbl]] <- tibble(
#     Stratum = r$lbl,
#     FS_ATT  = round(beta_fs, 5),
#     FS_SE   = round(se_fs,   5),
#     FS_F    = round(F_fs,    3),
#     RF_ATT  = round(beta_rf, 4),
#     RF_SE   = round(se_rf,   4),
#     AR_status = ar_status
#   )
# }
# tab_first_stage <- bind_rows(fs_rows)
# if (nrow(tab_first_stage) > 0) {
#   write_csv(tab_first_stage, "Tables_Final/tab_first_stage_v5.csv")
#   tryCatch(print(xtable(tab_first_stage,
#                         caption = "First-stage diagnostics by dose stratum.",
#                         label   = "tab:first_stage",
#                         digits  = c(0, 0, 5, 5, 3, 4, 4, 0)),
#                  file = "Tables_Final/tab_first_stage_v5.tex",
#                  include.rownames = FALSE),
#            error = function(e) NULL)
#   print(tab_first_stage)
# }

# ==============================================================================
# SECTION 9N: TWFE BENCHMARK ON COHORT OUTCOMES
# ==============================================================================
message(" SECTION 9N: TWFE Benchmark on cohort outcomes...")

twfe_gap   <- NULL
twfe_young <- NULL
tryCatch({
  twfe_gap_data <- df_cohort_wide %>%
    mutate(treated = as.integer(g > 0 & ano >= g))
  twfe_gap <- fixest::feols(wage_gap ~ treated | id_num + ano,
                            data = twfe_gap_data, cluster = ~id_num)

  twfe_young_data <- df_young_only %>%
    mutate(treated = as.integer(g > 0 & ano >= g))
  twfe_young <- fixest::feols(y_wage ~ treated | id_num + ano,
                              data = twfe_young_data, cluster = ~id_num)

  twfe_cohort_summary <- tibble(
    Model = c("TWFE/WageGap", "TWFE/Young"),
    Coef  = c(coef(twfe_gap)["treated"], coef(twfe_young)["treated"]),
    SE    = c(coeftable(twfe_gap)["treated", 2], coeftable(twfe_young)["treated", 2]),
    N     = c(nobs(twfe_gap), nobs(twfe_young))
  )
  write_csv(twfe_cohort_summary, "Tables_Final/tab_twfe_cohort_v5.csv")
  message("   TWFE cohort benchmark saved.")
}, error = function(e) message(sprintf("   TWFE cohort SKIPPED: %s", e$message)))

# ==============================================================================
# SECTION 9O: ALTERNATIVE CUTOFFS ON WAGE GAP
# ==============================================================================
message(" SECTION 9O: Alternative cutoffs on wage gap...")

alt_cohort_results <- list()
for (alt_name in c("q25_q75", "q33_q67")) {
  qs <- if (alt_name == "q25_q75") c(0.25, 0.75) else c(0.33, 0.67)
  tau_alt_low  <- quantile(df_2019$D_rt, qs[1], na.rm = TRUE)
  tau_alt_high <- quantile(df_2019$D_rt, qs[2], na.rm = TRUE)

  df_alt <- df_cohort_wide %>%
    left_join(
      df_final %>% select(id_microrregiao, ano, D_rt_2019) %>% distinct(),
      by = c("id_microrregiao", "ano")
    ) %>%
    mutate(
      alt_bin = case_when(
        D_rt_2019 > tau_alt_low & D_rt_2019 <= quantile(df_2019$D_rt, 0.50, na.rm = TRUE) ~ "Low",
        D_rt_2019 > tau_alt_high ~ "High",
        TRUE ~ NA_character_
      )
    )

  d_alt_low <- df_alt %>% filter(alt_bin == "Low" | g == 0) %>%
    mutate(y = wage_gap) %>% select(id_num, ano, g, y)
  res <- tryCatch(run_cs_preagg(d_alt_low, paste0("AltCut/", alt_name, "/Low")),
                  error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })
  if (!is.null(res)) alt_cohort_results[[paste0(alt_name, "_low")]] <- res

  d_alt_high <- df_alt %>% filter(alt_bin == "High" | g == 0) %>%
    mutate(y = wage_gap) %>% select(id_num, ano, g, y)
  res <- tryCatch(run_cs_preagg(d_alt_high, paste0("AltCut/", alt_name, "/High")),
                  error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })
  if (!is.null(res)) alt_cohort_results[[paste0(alt_name, "_high")]] <- res
}
message("   Alternative cutoff models on wage gap done.")

# ==============================================================================
# SECTION 9P: CONTINUOUS DOSE-RESPONSE -- contdid (v5 Fix 1: scalar target_parameter)
# ==============================================================================
message(" SECTION 9P: contdid (v5 Fix 1)...")

res_contdid_level <- NULL
res_contdid_slope <- NULL

tryCatch({
  if (!requireNamespace("contdid", quietly = TRUE)) {
    devtools::install_github("bcallaway11/contdid")
  }
  library(contdid)

  contdid_panel <- df_nonwhite %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    inner_join(
      df_final %>% mutate(id_num = as.integer(id_microrregiao)) %>%
        select(id_num, ano, D_rt),
      by = c("id_num", "ano")
    ) %>%
    filter(!is.na(y), !is.na(D_rt)) %>%
    rename(id = id_num, time = ano, dose = D_rt, G = g)

  message(sprintf("   contdid panel: %d obs, %d units, dose range [%.1f, %.1f]",
                  nrow(contdid_panel), n_distinct(contdid_panel$id),
                  min(contdid_panel$dose), max(contdid_panel$dose)))

  # v5 Fix 1: pass target_parameter as scalar string (NOT a vector)
  message("   -> contdid level...")
  res_contdid_level <- cont_did(
    yname = "y", tname = "time", idname = "id",
    dname = "dose", gname = "G",
    data = as.data.frame(contdid_panel),
    control_group = "notyettreated",
    target_parameter = "level"
  )
  tryCatch({
    p_cd_level <- plot(res_contdid_level)
    ggsave("Figures_Final/fig_contdid_level_v5.pdf", p_cd_level,
           width = 9, height = 5.5, dpi = 300)
    message("   contdid level figure saved.")
  }, error = function(e) message(sprintf("   contdid level plot failed: %s", e$message)))

  message("   -> contdid slope...")
  res_contdid_slope <- cont_did(
    yname = "y", tname = "time", idname = "id",
    dname = "dose", gname = "G",
    data = as.data.frame(contdid_panel),
    control_group = "notyettreated",
    target_parameter = "slope"
  )
  tryCatch({
    p_cd_slope <- plot(res_contdid_slope)
    ggsave("Figures_Final/fig_contdid_slope_v5.pdf", p_cd_slope,
           width = 9, height = 5.5, dpi = 300)
    message("   contdid slope figure saved.")
  }, error = function(e) message(sprintf("   contdid slope plot failed: %s", e$message)))

  message("   contdid succeeded.")
}, error = function(e) {
  message(sprintf("   contdid SKIPPED: %s -- falling back to fixest results.", safe_err_msg(e)))
})



# ==============================================================================
# SECTION 9Q: CONTINUOUS DOSE MODELS (fixest) -- MAXIMUM POWER
# ==============================================================================
message(" SECTION 9Q: Continuous dose models (fixest)...")

fit_lin      <- NULL; fit_log    <- NULL; fit_quad   <- NULL
fit_lin_dr   <- NULL; fit_young_log <- NULL
fit_bc_log   <- NULL; fit_bc_quad   <- NULL
cont_results <- NULL

tryCatch({
  library(fixest)

  # --- Pooled wages (ages 22-35, non-white) ---
  cont_data <- df_nonwhite %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    inner_join(
      df_final %>% mutate(id_num = as.integer(id_microrregiao)) %>%
        select(id_num, ano, D_rt, D_rt_2019),
      by = c("id_num", "ano")
    ) %>%
    mutate(post = as.integer(g > 0 & ano >= g))

  fit_lin  <- feols(y ~ D_rt:post | id_num + ano,
                    data = cont_data, cluster = ~id_num)
  fit_log  <- feols(y ~ I(log(1 + D_rt)):post | id_num + ano,
                    data = cont_data, cluster = ~id_num)
  fit_quad <- feols(y ~ D_rt:post + I(D_rt^2):post | id_num + ano,
                    data = cont_data, cluster = ~id_num)

  cont_data_dr <- cont_data %>% left_join(df_covars_dr, by = "id_num")
  fit_lin_dr <- feols(
    y ~ D_rt:post + post:log_pop_18_24 + post:share_nonwhite +
      post:share_no_formal_income + post:avg_hh_income_pc_mw + post:literacy_rate |
      id_num + ano,
    data = cont_data_dr, cluster = ~id_num
  )

  message("\n   === CONTINUOUS DOSE RESULTS (POOLED WAGES) ===")
  ct_lin  <- coeftable(fit_lin);  ct_log  <- coeftable(fit_log)
  ct_quad <- coeftable(fit_quad); ct_dr   <- coeftable(fit_lin_dr)

  cont_results <- bind_rows(
    tibble(Model="Linear D*post",    Coef=ct_lin[1,1],  SE=ct_lin[1,2],
           t=ct_lin[1,3],  p=ct_lin[1,4],  N=nobs(fit_lin)),
    tibble(Model="Log(1+D)*post",    Coef=ct_log[1,1],  SE=ct_log[1,2],
           t=ct_log[1,3],  p=ct_log[1,4],  N=nobs(fit_log)),
    tibble(Model="Quad D*post",      Coef=ct_quad[1,1], SE=ct_quad[1,2],
           t=ct_quad[1,3], p=ct_quad[1,4], N=nobs(fit_quad)),
    tibble(Model="Quad D^2*post",    Coef=ct_quad[2,1], SE=ct_quad[2,2],
           t=ct_quad[2,3], p=ct_quad[2,4], N=nobs(fit_quad)),
    tibble(Model="Linear+DR D*post", Coef=ct_dr[1,1],   SE=ct_dr[1,2],
           t=ct_dr[1,3],   p=ct_dr[1,4],   N=nobs(fit_lin_dr))
  )
  write_csv(cont_results, "Tables_Final/tab_continuous_dose_v5.csv")
  message("   Continuous dose results saved.")
  print(cont_results)

  # Delta-method: peak dose D* = -beta1 / (2 * beta2) for quadratic model
  if (!is.null(fit_quad)) {
    b1 <- ct_quad[1, 1]
    b2 <- ct_quad[2, 1]
    if (abs(b2) > 1e-10 && b2 < 0) {
      D_star <- -b1 / (2 * b2)
      # Delta-method SE: gradient = [-1/(2b2), b1/(2b2^2)]
      vcv_quad <- vcov(fit_quad)[1:2, 1:2]
      grad <- c(-1 / (2 * b2), b1 / (2 * b2^2))
      D_star_se <- sqrt(as.numeric(t(grad) %*% vcv_quad %*% grad))
      message(sprintf("   Peak dose D* = %.2f (SE = %.2f) [Delta method]", D_star, D_star_se))
    }
  }

  # --- Young-only wages (ages 22-28) ---
  if (exists("df_young_only") && nrow(df_young_only) > 0) {
    cont_young_data <- df_young_only %>%
      inner_join(
        df_final %>% mutate(id_num = as.integer(id_microrregiao)) %>%
          select(id_num, ano, D_rt),
        by = c("id_num", "ano")
      ) %>%
      mutate(post = as.integer(g > 0 & ano >= g), y = y_wage)

    fit_young_log <- feols(y ~ I(log(1 + D_rt)):post | id_num + ano,
                           data = cont_young_data, cluster = ~id_num)
    ct_yl <- coeftable(fit_young_log)
    message(sprintf("   Young/Log: coef=%.6f SE=%.6f t=%.3f p=%.4f",
                    ct_yl[1,1], ct_yl[1,2], ct_yl[1,3], ct_yl[1,4]))
  }

  # --- Birth-cohort wage gap ---
  if (!is.null(df_bc_wide) && nrow(df_bc_wide) > 0) {
    cont_bc_data <- df_bc_wide %>%
      mutate(post = as.integer(g > 0 & ano >= g))

    fit_bc_log <- feols(wage_gap_bc ~ I(log(1 + D_rt)):post | id_num + ano,
                        data = cont_bc_data, cluster = ~id_num)
    ct_bcl <- coeftable(fit_bc_log)
    message(sprintf("   BirthCohort/Log: coef=%.6f SE=%.6f t=%.3f p=%.4f",
                    ct_bcl[1,1], ct_bcl[1,2], ct_bcl[1,3], ct_bcl[1,4]))

    fit_bc_quad <- feols(wage_gap_bc ~ D_rt:post + I(D_rt^2):post | id_num + ano,
                         data = cont_bc_data, cluster = ~id_num)
    ct_bcq <- coeftable(fit_bc_quad)
    message(sprintf("   BirthCohort/Quad: D coef=%.6f (t=%.3f), D^2 coef=%.8f (t=%.3f)",
                    ct_bcq[1,1], ct_bcq[1,3], ct_bcq[2,1], ct_bcq[2,3]))
  }

}, error = function(e) message(sprintf("   SECTION 9Q FAILED: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9R: contdid ON COHORT WAGE GAP (secondary)
# ==============================================================================
message(" SECTION 9R: contdid on cohort wage gap...")

res_contdid_gap <- NULL
tryCatch({
  if (!requireNamespace("contdid", quietly = TRUE)) stop("contdid not installed")
  library(contdid)

  contdid_gap_panel <- df_cohort_wide %>%
    select(id_num, ano, g, wage_gap) %>%
    inner_join(
      df_final %>% mutate(id_num = as.integer(id_microrregiao)) %>%
        select(id_num, ano, D_rt),
      by = c("id_num", "ano")
    ) %>%
    filter(!is.na(wage_gap), !is.na(D_rt)) %>%
    rename(id = id_num, time = ano, dose = D_rt, G = g, y = wage_gap)

  res_contdid_gap <- cont_did(
    yname = "y", tname = "time", idname = "id",
    dname = "dose", gname = "G",
    data = as.data.frame(contdid_gap_panel),
    control_group = "notyettreated",
    target_parameter = "level"
  )
  tryCatch({
    p_cdg <- plot(res_contdid_gap)
    ggsave("Figures_Final/fig_contdid_gap_v5.pdf", p_cdg,
           width = 9, height = 5.5, dpi = 300)
    message("   contdid wage gap figure saved.")
  }, error = function(e) message(sprintf("   contdid gap plot failed: %s", e$message)))
}, error = function(e) {
  message(sprintf("   contdid gap SKIPPED: %s", safe_err_msg(e)))
})

# ==============================================================================
# SECTION 9S: BIRTH-COHORT CS MODELS (Duflo-strict)
# ==============================================================================
message(" SECTION 9S: Birth-cohort CS models...")

res_bc_gap_low <- NULL; res_bc_gap_high <- NULL
res_bc_exposed_low <- NULL; res_bc_control_low <- NULL

if (!is.null(bc_gap_low) && nrow(bc_gap_low) > 0) {
  message("   -> Birth-cohort gap / Low Dose...")
  res_bc_gap_low <- tryCatch(run_cs_preagg(bc_gap_low, "BC_Gap/Low"),
                             error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

  message("   -> Birth-cohort gap / High Dose...")
  res_bc_gap_high <- tryCatch(run_cs_preagg(bc_gap_high, "BC_Gap/High"),
                              error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

  message("   -> Birth-cohort exposed / Low Dose...")
  res_bc_exposed_low <- tryCatch(run_cs_preagg(bc_exposed_low, "BC_Exposed/Low"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

  message("   -> Birth-cohort control / Low Dose [PLACEBO]...")
  res_bc_control_low <- tryCatch(run_cs_preagg(bc_control_low, "BC_Control/Low"),
                                 error = function(e) { message(sprintf("   WARNING: %s", safe_err_msg(e))); NULL })

  if (!is.null(res_bc_gap_low) && !is.null(res_bc_gap_high)) {
    p <- plot_event_study(res_bc_gap_low$tidy, res_bc_gap_high$tidy,
                          "Low Dose", "High Dose", "#E69F00", "#0072B2",
                          "Birth-Cohort Wage Gap: Exposed (1987-93) vs Control (1975-82)",
                          ylab = "ATT on wage gap (exposed - control)")
    ggsave("Figures_Final/fig_birthcohort_gap_v5.pdf", p, width = 9, height = 5.5, dpi = 300)
    message("   Birth-cohort gap figure saved.")
  }

  if (!is.null(res_bc_exposed_low) && !is.null(res_bc_control_low)) {
    p <- plot_event_study(res_bc_exposed_low$tidy, res_bc_control_low$tidy,
                          "Exposed (born 1987-93)", "Control (born 1975-82) [Placebo]",
                          "#2166ac", "#999999",
                          "Birth-Cohort Decomposition: Low Dose",
                          ylab = "ATT on log(wages)")
    ggsave("Figures_Final/fig_birthcohort_decomp_v5.pdf", p, width = 9, height = 5.5, dpi = 300)
    message("   Birth-cohort decomposition figure saved.")
  }
} else {
  message("   SKIPPED: Birth-cohort data not available.")
}



# ==============================================================================
# SECTION 3H: CENSO DA EDUCACAO SUPERIOR (Higher Education Census -- INEP)
# ==============================================================================
# Source: basedosdados.br_inep_censo_educacao_superior.curso
# Available years: 2009-2024. We use 2009 as the pre-treatment covariate
# snapshot (earliest available year overlapping the study period).
#
# NOTE: Treatment timing g ranges 2006-2012 for most microregions. The 2009
# snapshot is post-treatment for some early-treated units; treat as a
# cross-sectional HE supply measure, not a clean pre-treatment variable.
# Use with caution and note this limitation in the text.
# ==============================================================================
message(" SECTION 3H: Downloading Higher Education Census (Censo Educacao Superior)...")

df_censo_educ <- load_or_fetch(
  "Documentation_Final/cache_censo_educ.rds",
  function() {
    q_censo <- "
      SELECT
        id_municipio,
        SUM(COALESCE(quantidade_ingressantes, 0))                      AS total_ingressantes,
        SUM(CASE WHEN rede = 'Privada'
                 THEN COALESCE(quantidade_matriculas, 0) ELSE 0 END)   AS private_matriculas,
        SUM(COALESCE(quantidade_matriculas, 0))                        AS total_matriculas,
        SUM(COALESCE(quantidade_ingressantes_financiamento_nao_reembolsavel_prouni_integral, 0) +
            COALESCE(quantidade_ingressantes_financiamento_nao_reembolsavel_prouni_parcial, 0))
                                                                       AS prouni_ingressantes,
        SUM(COALESCE(quantidade_ingressantes_preta,   0) +
            COALESCE(quantidade_ingressantes_parda,   0) +
            COALESCE(quantidade_ingressantes_indigena, 0))             AS nonwhite_ingressantes
      FROM `basedosdados.br_inep_censo_educacao_superior.curso`
      WHERE ano = 2009
        AND id_municipio IS NOT NULL
      GROUP BY id_municipio
    "
    basedosdados::read_sql(q_censo)
  },
  "censo_educ_2009"
)

message(sprintf("   Censo Educacao Superior 2009: %d municipalities", nrow(df_censo_educ)))
gc()

# ==============================================================================
# SECTION 6B: EDUCATION CENSUS COVARIATES (microregion-level)
# ==============================================================================
message(" SECTION 6B: Building Education Census covariates at microregion level...")

df_censo_micro <- df_censo_educ %>%
  inner_join(
    df_geo %>% distinct(id_municipio, id_microrregiao),
    by = "id_municipio"
  ) %>%
  group_by(id_microrregiao) %>%
  summarise(
    total_ingressantes    = sum(total_ingressantes,    na.rm = TRUE),
    private_matriculas    = sum(private_matriculas,    na.rm = TRUE),
    total_matriculas      = sum(total_matriculas,      na.rm = TRUE),
    prouni_ingressantes   = sum(prouni_ingressantes,   na.rm = TRUE),
    nonwhite_ingressantes = sum(nonwhite_ingressantes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Join youth population for per-1000 denominators
  left_join(
    df_pop_micro %>% rename(id_microrregiao = id_microrregiao),
    by = "id_microrregiao"
  ) %>%
  mutate(
    id_num = as.integer(id_microrregiao),
    # HE absorption capacity: new enrolments per 1,000 youth (2009)
    college_density = ifelse(pop > 0, total_ingressantes / pop * 1000, NA_real_),
    # Share of total enrolments in private institutions (supply-side control)
    private_enroll_share = ifelse(total_matriculas > 0,
                                  private_matriculas / total_matriculas, NA_real_),
    # ProUni take-up rate per 1,000 youth -- informative but NOTE: correlated with dose
    prouni_per_1000_censo = ifelse(pop > 0, prouni_ingressantes / pop * 1000, NA_real_),
    # Non-white share of new entrants (pre-treatment diversity of HE access)
    nonwhite_enroll_share = ifelse(total_ingressantes > 0,
                                   nonwhite_ingressantes / total_ingressantes, NA_real_)
  ) %>%
  select(id_num, college_density, private_enroll_share,
         prouni_per_1000_censo, nonwhite_enroll_share)

message(sprintf("   Education Census covariates: %d microregions with HE data",
                sum(!is.na(df_censo_micro$college_density))))

# Merge into the enriched DR covariate frame
df_covars_dr_enrich <- df_covars_dr %>%
  left_join(df_censo_micro, by = "id_num") %>%
  mutate(
    # Impute zeros for microregions with no HE institutions in 2009
    college_density      = replace_na(college_density,      0),
    private_enroll_share = replace_na(private_enroll_share, 0),
    nonwhite_enroll_share = replace_na(nonwhite_enroll_share, 0)
    # Note: prouni_per_1000_censo left as NA where missing to avoid
    # confounding with the treatment variable itself
  )

message(sprintf("   df_covars_dr_enrich: %d microregions | %d columns",
                nrow(df_covars_dr_enrich), ncol(df_covars_dr_enrich)))

# Enriched DR formulas (Education Census adds college_density + private_enroll_share)
# Low-Dose: 6 covariates -- feasible with ~167 treated + control units
dr_xformla_low_enrich  <- ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate +
                             share_in_school + college_density + private_enroll_share

# High-Dose: 3 covariates -- conservative to avoid singular matrix (~140 units)
dr_xformla_high_enrich <- ~ log_pop_18_24 + avg_hh_income_pc_mw + college_density

message("   Enriched DR formulas defined (Low: 6 covars | High: 3 covars).")

# ==============================================================================
# SECTION 9F_v2: ENRICHED DR MODELS WITH EDUCATION CENSUS COVARIATES
# ==============================================================================
message(" SECTION 9F_v2: Re-running DR models enriched with Education Census...")

res_low_dr_censo  <- tryCatch({
  message("   -> Low Dose / DR + Censo (6 covars, Sant'Anna-Zhao)...")
  run_cs(data_low, "log_wage_sm", label = "Low/DR/Censo",
         xformla  = dr_xformla_low_enrich,
         covars_df = df_covars_dr_enrich,
         est_method = "dr")
}, error = function(e) {
  message(sprintf("   WARNING: Low/DR/Censo failed: %s", safe_err_msg(e)))
  NULL
})

res_high_dr_censo <- tryCatch({
  message("   -> High Dose / DR + Censo (3 covars, Sant'Anna-Zhao)...")
  run_cs(data_high, "log_wage_sm", label = "High/DR/Censo",
         xformla  = dr_xformla_high_enrich,
         covars_df = df_covars_dr_enrich,
         est_method = "dr")
}, error = function(e) {
  message(sprintf("   WARNING: High/DR/Censo failed: %s", safe_err_msg(e)))
  NULL
})

# Flag divergence from the base DR estimates
if (!is.null(res_low_dr_censo) && !is.null(res_low_dr)) {
  diff_att <- abs(res_low_dr$simple$overall.att - res_low_dr_censo$simple$overall.att)
  se_pool  <- sqrt(res_low_dr$simple$overall.se^2 + res_low_dr_censo$simple$overall.se^2)
  message(sprintf("   Low/DR vs Low/DR/Censo ATT diff: %.4f (%.1f pooled SEs)",
                  diff_att, if (se_pool > 0) diff_att / se_pool else NA_real_))
}

# HonestDiD on enriched Low/DR/Censo (using same run_honest_v5 helper)
if (requireNamespace("HonestDiD", quietly = TRUE) && !is.null(res_low_dr_censo)) {
  hr_censo <- run_honest_v5(res_low_dr_censo$es, "Low/DR/Censo")
  honest_results[["Low/DR/Censo"]] <- hr_censo
  if (!is.null(hr_censo$results) && nrow(hr_censo$results) > 0) {
    tryCatch({
      p_hc <- ggplot(hr_censo$results, aes(x = Mbar, y = estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.20, fill = "#009E73") +
        geom_line(color = "#009E73", linewidth = 1) +
        geom_point(color = "#009E73", size = 2) +
        labs(
          title    = "HonestDiD Sensitivity: Low/DR/Censo",
          subtitle = sprintf("Breakdown Mbar = %.2f | Enriched with Education Census",
                             hr_censo$breakdown),
          x = expression(bar(M) ~ "(relative magnitudes)"),
          y = "ATT (robust CI)"
        )
      ggsave("Figures_Final/fig_honestdid_low_dr_censo_v5.pdf", p_hc,
             width = 7, height = 5, dpi = 300)
      message("   HonestDiD Low/DR/Censo figure saved.")
    }, error = function(e) {
      message(sprintf("   WARNING: HonestDiD censo plot failed: %s", safe_err_msg(e)))
    })
  }
}

message(" -> Education Census covariate integration complete.")


# ==============================================================================
# SECTION 9L: ATT SUMMARY TABLE (all models)
# ==============================================================================
message(" Building ATT summary table (all models)...")

all_results_v5 <- Filter(Negate(is.null), list(
  # Main wage models
  res_low, res_mid, res_high, res_cen, res_mar,
  # Belíndia 5-region split
  res_region_N, res_region_NE, res_region_CO, res_region_SE, res_region_S,
  # Employment models
  res_low_emp, res_high_emp, res_cen_emp, res_mar_emp,
  # Gender heterogeneity
  res_male_low, res_female_low, res_male_high, res_female_high,
  # Race heterogeneity
  res_white_low, res_nonwhite_low, res_white_high, res_nonwhite_high,
  # Lag robustness
  res_low_lag3, res_high_lag3,
  # Doubly-robust (base + enriched with Education Census)
  res_low_dr, res_high_dr,
  res_low_dr_censo, res_high_dr_censo,
  # SUTVA buffer
  res_low_buffer, res_high_buffer,
  # Placebo cohort (Fix 6)
  res_placebo_low, res_placebo_high,
  # Alternative cutoffs
  res_low_25_75, res_high_25_75, res_low_33_67, res_high_33_67,
  # Cohort Duflo-style: gap, young, old
  res_gap_low, res_gap_high,
  res_young_low, res_young_high,
  res_old_low, res_old_high,
  # Extended cohort: college gap, gender, race
  res_college_low, res_college_high,
  res_gap_male_low, res_gap_male_high,
  res_gap_female_low, res_gap_female_high,
  res_gap_white_low, res_gap_white_high,
  res_gap_nonwhite_low, res_gap_nonwhite_high,
  # Birth-cohort (Duflo-strict)
  res_bc_gap_low, res_bc_gap_high,
  res_bc_exposed_low, res_bc_control_low,
  # Triple interaction (Fix 4)
  res_triple_young_poor,
  # Alt cohort cutoffs
  alt_cohort_results[["q25_q75_low"]], alt_cohort_results[["q25_q75_high"]],
  alt_cohort_results[["q33_q67_low"]], alt_cohort_results[["q33_q67_high"]]
))

tab_att <- bind_rows(lapply(all_results_v5, build_att_row))

# Add TWFE rows
for (tw in list(list(fit = twfe_fit, label = "TWFE/Wages"),
                list(fit = twfe_emp,  label = "TWFE/Employment"),
                list(fit = twfe_gap,  label = "TWFE/WageGap"),
                list(fit = twfe_young,label = "TWFE/Young"))) {
  if (!is.null(tw$fit)) {
    tc <- tryCatch(fixest::coeftable(tw$fit), error = function(e) NULL)
    if (!is.null(tc)) {
      tab_att <- bind_rows(tab_att, tibble(
        Model      = tw$label,
        ATT        = round(tc[1, "Estimate"],   4),
        SE         = round(tc[1, "Std. Error"],  4),
        CI_low     = round(tc[1, "Estimate"] - 1.96 * tc[1, "Std. Error"], 4),
        CI_high    = round(tc[1, "Estimate"] + 1.96 * tc[1, "Std. Error"], 4),
        t_val      = round(tc[1, "t value"],    3),
        p_val      = round(tc[1, "Pr(>|t|)"],   4),
        Pretrend_p = NA_real_,
        N_obs      = as.integer(stats::nobs(tw$fit)),
        N_treated  = NA_integer_,
        N_notyet   = NA_integer_
      ))
    }
  }
}

write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
print(xtable(tab_att,
             caption = "Summary of ATT Estimates Across All Specifications (v5)",
             digits  = c(0, 0, 4, 4, 4, 4, 3, 4, 4, 0, 0, 0)),
      file = "Tables_Final/tab_att_v5.tex", include.rownames = FALSE)
message(sprintf("   ATT summary table saved: %d models.", nrow(tab_att)))
print(tab_att)

# ---- CONSOLE PLOT: ATT forest plot (key models) ---------------------------------
# tab_att is already built above; filter to the key model subset directly
key_models_in_forest <- c(
  "Low/Wages", "Mid/Wages", "High/Wages",
  "Male/Low", "Female/Low",
  "Non-white/Low", "White/Low",
  "Central/Wages", "Marginal/Wages",
  "Region/Norte", "Region/Nordeste", "Region/CentroOeste",
  "Region/Sudeste", "Region/Sul",
  "Low/DR", "Low/DR/Censo", "Low/DR/Policy",
  "Low/HighQualityHEI"
)
df_forest <- tab_att %>%
  filter(Model %in% key_models_in_forest) %>%
  select(label = Model, att = ATT, ci_lo = CI_low, ci_hi = CI_high)

if (nrow(df_forest) > 0) {
  p_forest <- ggplot(df_forest,
                     aes(x = att, y = reorder(label, att))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.3, color = "#0072B2", linewidth = 0.8) +
    geom_point(color = "#0072B2", size = 3) +
    labs(
      title    = "ATT Estimates: Key Models (95% CI)",
      subtitle = "Callaway-Sant'Anna staggered DiD | Not-yet-treated controls",
      x = "ATT (log wage effect)", y = NULL
    )
  print(p_forest)
  message("   Console plot: ATT forest plot printed.")
}

# ==============================================================================
# SECTION 9L2: ATT SUMMARY TABLE
# ------------------------------------------------------------------------------
# Write the current tab_att (raw ATT/SE/p-values) to disk.
# p_val uses the normal approximation: 2*pnorm(-|ATT/SE|), which is standard
# for large-sample clustered inference with many microregion clusters.
# ==============================================================================
message(" Writing ATT summary table...")

write_csv(tab_att, "Tables_Final/tab_att_v5.csv")

# Key-models view for the main text
tab_att_key <- tab_att %>%
  filter(Model %in% c("Low/Wages", "Mid/Wages", "High/Wages",
                      "Low/DR", "High/DR", "Low/DR/Censo", "High/DR/Censo",
                      "Male/Low", "Female/Low", "Non-white/Low", "White/Low",
                      "Central/Wages", "Marginal/Wages",
                      "Low/Buffer", "High/Buffer",
                      "Placebo/Low", "Placebo/High",
                      "TWFE/Wages"))
write_csv(tab_att_key, "Tables_Final/tab_att_key_v5.csv")
tryCatch(print(xtable(tab_att_key,
                      caption = "Key ATT Estimates (raw normal-approximation $p$-values)",
                      label   = "tab:att_key",
                      digits  = c(0, 0, 4, 4, 3, 4, 4, 4, 4, 0, 0, 0)),
               file = "Tables_Final/tab_att_key_v5.tex",
               include.rownames = FALSE),
         error = function(e) message(sprintf("   WARNING: tab_att_key tex export failed: %s",
                                             safe_err_msg(e))))
message(sprintf("   ATT summary written: %d specifications.", nrow(tab_att)))

# ==============================================================================
# SECTION 9L3: HonestDiD BREAKDOWN SUMMARY -- CSV export
# ------------------------------------------------------------------------------
# For each HonestDiD run we already computed `breakdown` = largest Mbar for
# which the lower CI bound remains > 0 (so the sign of the point estimate
# survives allowing for non-linear pre-trend deviations of magnitude Mbar
# times the largest observed pre-trend coefficient).  Export to a clean CSV.
# ==============================================================================
tab_honestdid_breakdown <- bind_rows(lapply(names(honest_results), function(nm) {
  hr <- honest_results[[nm]]
  tibble(
    Specification  = nm,
    Breakdown_Mbar = ifelse(is.null(hr$breakdown), NA_real_, hr$breakdown),
    Full_VCov      = isTRUE(hr$used_full_vcov),
    N_Mbar_Tested  = ifelse(is.null(hr$results), NA_integer_, nrow(hr$results))
  )
}))
if (nrow(tab_honestdid_breakdown) > 0) {
  write_csv(tab_honestdid_breakdown, "Tables_Final/tab_honestdid_breakdown_v5.csv")
  tryCatch(print(xtable(tab_honestdid_breakdown,
                        caption = "HonestDiD Breakdown $\\bar{M}$: largest relative-magnitude perturbation for which the robust 95\\% CI still excludes zero",
                        label   = "tab:honest_breakdown",
                        digits  = c(0, 0, 2, 0, 0)),
                 file = "Tables_Final/tab_honestdid_breakdown_v5.tex",
                 include.rownames = FALSE),
           error = function(e) NULL)
  message(sprintf("   HonestDiD breakdown summary: %d specs exported.",
                  nrow(tab_honestdid_breakdown)))
  print(tab_honestdid_breakdown)
}

# ==============================================================================
# SECTION 9L4 removed: SPT balance table relied on Callaway,
#  Goodman-Bacon and Sant'Anna (2024) Assumption 3, which is not
#  defended in the final draft.
# ------------------------------------------------------------------------------
# SECTION 9L4: STRONG PARALLEL TRENDS -- DISABLED IN FINAL DRAFT
# ------------------------------------------------------------------------------
# message(" SECTION 9L4: Strong-parallel-trends balance table...")
#
# tab_spt_balance <- tryCatch({
#   if (exists("df_covars_dr_enrich") && !is.null(df_covars_dr_enrich)) {
#     bin_lookup <- df_final %>% distinct(id_microrregiao, dose_bin) %>%
#       mutate(id_num = match(as.character(id_microrregiao),
#                              as.character(sort(unique(df_final$id_microrregiao)))))
#     df_bal <- df_covars_dr_enrich %>% left_join(bin_lookup %>% select(id_num, dose_bin),
#                                                  by = "id_num")
#     numeric_cols <- names(df_bal)[sapply(df_bal, is.numeric)]
#     numeric_cols <- setdiff(numeric_cols, c("id_num"))
#     bal_rows <- lapply(numeric_cols, function(cc) {
#       x <- df_bal[[cc]]
#       pooled_sd <- function(a, b) sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
#       std_diff <- function(a, b) {
#         s <- pooled_sd(a, b); if (is.na(s) || s == 0) NA_real_ else
#           (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / s
#       }
#       xa <- x[df_bal$dose_bin == "Low"]
#       xb <- x[df_bal$dose_bin == "High"]
#       xc <- x[df_bal$dose_bin == "Mid"]
#       xn <- x[is.na(df_bal$dose_bin)]
#       tibble(
#         Covariate           = cc,
#         Mean_Low            = round(mean(xa, na.rm = TRUE), 4),
#         Mean_Mid            = round(mean(xc, na.rm = TRUE), 4),
#         Mean_High           = round(mean(xb, na.rm = TRUE), 4),
#         Mean_NotYet         = round(mean(xn, na.rm = TRUE), 4),
#         StdDiff_Low_High    = round(std_diff(xa, xb), 3),
#         StdDiff_Low_NotYet  = round(std_diff(xa, xn), 3),
#         StdDiff_High_NotYet = round(std_diff(xb, xn), 3)
#       )
#     })
#     do.call(bind_rows, bal_rows)
#   } else NULL
# }, error = function(e) { NULL })
#
# if (!is.null(tab_spt_balance)) {
#   write_csv(tab_spt_balance, "Tables_Final/tab_spt_balance_v5.csv")
#   ...
# }

# ==============================================================================
# SECTION 9L5: GOODMAN-BACON DECOMPOSITION OF THE TWFE ESTIMATOR
# ------------------------------------------------------------------------------
# The "forbidden-comparison" bias of TWFE in staggered designs (Goodman-Bacon
# 2021) decomposes the TWFE coefficient into a weighted average of 2x2 DD
# estimates across every pair of treated/untreated subgroups.  The
# bacondecomp package provides this decomposition directly.  We report
# weights and sub-effects, and flag the share of the TWFE coefficient driven
# by forbidden "already-treated" comparisons.
# ==============================================================================
message(" SECTION 9L5: Goodman-Bacon decomposition of TWFE...")

bacon_decomp <- NULL
tryCatch({
  if (!requireNamespace("bacondecomp", quietly = TRUE)) {
    message("   bacondecomp not installed; attempting install...")
    tryCatch(install.packages("bacondecomp", repos = "https://cloud.r-project.org"),
             error = function(e) message("   bacondecomp install failed; skipping."))
  }
  if (requireNamespace("bacondecomp", quietly = TRUE)) {
    twfe_bacon_data <- df_final %>%
      mutate(treated = as.integer(g > 0 & ano >= g)) %>%
      filter(!is.na(log_wage_sm)) %>%
      group_by(id_microrregiao, ano) %>%
      summarise(log_wage_sm = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
                treated = max(treated, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(id_num = as.integer(factor(id_microrregiao)))
    # Balanced panel is required by bacondecomp
    n_years <- n_distinct(twfe_bacon_data$ano)
    complete_ids <- twfe_bacon_data %>%
      group_by(id_num) %>% summarise(nn = n(), .groups = "drop") %>%
      filter(nn == n_years) %>% pull(id_num)
    twfe_bacon_data <- twfe_bacon_data %>% filter(id_num %in% complete_ids)

    bacon_decomp <- bacondecomp::bacon(
      log_wage_sm ~ treated,
      data    = as.data.frame(twfe_bacon_data),
      id_var  = "id_num",
      time_var= "ano"
    )
    bacon_tab <- as_tibble(bacon_decomp)
    write_csv(bacon_tab, "Tables_Final/tab_bacon_v5.csv")
    # Summarise by type
    bacon_summary <- bacon_tab %>%
      group_by(type) %>%
      summarise(weight_sum = sum(weight, na.rm = TRUE),
                att_weighted = sum(estimate * weight, na.rm = TRUE) /
                               pmax(sum(weight, na.rm = TRUE), 1e-12),
                .groups = "drop")
    write_csv(bacon_summary, "Tables_Final/tab_bacon_summary_v5.csv")
    tryCatch(print(xtable(bacon_summary,
                          caption = "Goodman-Bacon (2021) decomposition of the naive TWFE wage coefficient.  The `Earlier vs Later Treated' row is the forbidden comparison that drives the sign bias.",
                          label   = "tab:bacon",
                          digits  = c(0, 0, 4, 4)),
                   file = "Tables_Final/tab_bacon_summary_v5.tex",
                   include.rownames = FALSE),
             error = function(e) NULL)
    print(bacon_summary)
    message("   Goodman-Bacon decomposition exported.")
  }
}, error = function(e) message(sprintf("   Bacon decomposition failed: %s",
                                        safe_err_msg(e))))

# ==============================================================================
# SECTION 10: PUBLICATION FIGURES (all with _v5 suffix)
# ==============================================================================
message(" [11/12] Generating Publication Figures...")

# Figure 1: Dose-Response wages -- High vs Low
p1 <- plot_event_study(
  res_low$tidy, res_high$tidy,
  "Low Dose (p20-p50)", "High Dose (>p75)",
  "#E69F00", "#0072B2",
  "Dose-Response: High vs. Low Exposure Intensity"
)
print(p1)
ggsave("Figures_Final/fig_1_dose_wages_v5.pdf", p1, width = 9, height = 5.5, dpi = 300)

# Figure 2: Regional divergence -- Central vs Marginal
if (!is.null(res_cen) && !is.null(res_mar)) {
  p2 <- plot_event_study(
    res_cen$tidy, res_mar$tidy,
    "Central (South/Southeast)", "Marginal (North/Northeast)",
    "#2980b9", "#e74c3c",
    "Regional Divergence: Centre vs. Periphery"
  )
  print(p2)
  ggsave("Figures_Final/fig_2_region_wages_v5.pdf", p2, width = 9, height = 5.5, dpi = 300)
}

# Figure 3: Employment -- High vs Low
p3 <- plot_event_study(
  res_low_emp$tidy, res_high_emp$tidy,
  "Low Dose", "High Dose",
  "#E69F00", "#0072B2",
  "Employment Margin: High vs. Low Exposure",
  ylab = "ATT on log(formal employment)"
)
print(p3)
ggsave("Figures_Final/fig_3_dose_emp_v5.pdf", p3, width = 9, height = 5.5, dpi = 300)

# Figure 4a-b: Gender heterogeneity
if (!is.null(res_male_low) && !is.null(res_female_low)) {
  p4a <- plot_event_study(res_male_low$tidy, res_female_low$tidy,
                          "Male", "Female", "#1b9e77", "#d95f02",
                          "Gender Heterogeneity: Low Dose")
  print(p4a)
  ggsave("Figures_Final/fig_4a_gender_low_v5.pdf", p4a, width = 9, height = 5.5, dpi = 300)
}
if (!is.null(res_male_high) && !is.null(res_female_high)) {
  p4b <- plot_event_study(res_male_high$tidy, res_female_high$tidy,
                          "Male", "Female", "#1b9e77", "#d95f02",
                          "Gender Heterogeneity: High Dose")
  print(p4b)
  ggsave("Figures_Final/fig_4b_gender_high_v5.pdf", p4b, width = 9, height = 5.5, dpi = 300)
}

# Figure 5a-b: Race heterogeneity
if (!is.null(res_white_low) && !is.null(res_nonwhite_low)) {
  p5a <- plot_event_study(res_white_low$tidy, res_nonwhite_low$tidy,
                          "White", "Non-white", "#7570b3", "#e7298a",
                          "Ethnicity Heterogeneity: Low Dose")
  print(p5a)
  ggsave("Figures_Final/fig_5a_race_low_v5.pdf", p5a, width = 9, height = 5.5, dpi = 300)
}
if (!is.null(res_white_high) && !is.null(res_nonwhite_high)) {
  p5b <- plot_event_study(res_white_high$tidy, res_nonwhite_high$tidy,
                          "White", "Non-white", "#7570b3", "#e7298a",
                          "Ethnicity Heterogeneity: High Dose")
  print(p5b)
  ggsave("Figures_Final/fig_5b_race_high_v5.pdf", p5b, width = 9, height = 5.5, dpi = 300)
}

# Figure 6: Lag robustness
if (!is.null(res_low_lag3)) {
  p6 <- plot_event_study(res_low$tidy, res_low_lag3$tidy,
                         "Baseline (lag=4)", "Robustness (lag=3)",
                         "#0072B2", "#D55E00",
                         "Robustness: Sensitivity to Lag Length (Low Dose)")
  print(p6)
  ggsave("Figures_Final/fig_6_lag_robustness_v5.pdf", p6, width = 9, height = 5.5, dpi = 300)
}

# Figure 7a-b: Unconditional vs DR (base + census-enriched)
if (!is.null(res_low_dr)) {
  p7a <- plot_event_study(res_low$tidy, res_low_dr$tidy,
                          "Unconditional", "Doubly-Robust (base)",
                          "#E69F00", "#009E73",
                          "Sensitivity: Unconditional vs. DR -- Low Dose")
  print(p7a)
  ggsave("Figures_Final/fig_7a_dr_low_v5.pdf", p7a, width = 9, height = 5.5, dpi = 300)
}
if (!is.null(res_low_dr_censo)) {
  p7b <- plot_event_study(res_low_dr$tidy, res_low_dr_censo$tidy,
                          "DR (base 4 covars)", "DR + Education Census",
                          "#009E73", "#CC79A7",
                          "Sensitivity: DR vs. DR + Education Census -- Low Dose")
  ggsave("Figures_Final/fig_7b_dr_censo_low_v5.pdf", p7b, width = 9, height = 5.5, dpi = 300)
}

# Figure 8: SUTVA buffer-zone exclusion
if (!is.null(res_low_buffer)) {
  p8 <- plot_event_study(res_low$tidy, res_low_buffer$tidy,
                         "Baseline", "Buffer Excluded",
                         "#E69F00", "#CC79A7",
                         "SUTVA Robustness: Buffer-Zone Exclusion (Low Dose)")
  print(p8)
  ggsave("Figures_Final/fig_8_buffer_low_v5.pdf", p8, width = 9, height = 5.5, dpi = 300)
}

# Figure 9: Alternative cutoffs
if (!is.null(res_low_25_75)) {
  p9 <- plot_event_study(res_low$tidy, res_low_25_75$tidy,
                         "Baseline (p20/p50)", "Alt (p25/p75)",
                         "#E69F00", "#56B4E9",
                         "Cutoff Robustness: Low Dose")
  print(p9)
  ggsave("Figures_Final/fig_9_altcut_low_v5.pdf", p9, width = 9, height = 5.5, dpi = 300)
}

# Figure 10: Placebo cohort (ages 36-50)
if (!is.null(res_placebo_low)) {
  p10 <- plot_event_study(res_low$tidy, res_placebo_low$tidy,
                          "Main (ages 22-35)", "Placebo (ages 36-50)",
                          "#0072B2", "#999999",
                          "Placebo Test: Main vs. Non-Exposed Cohort (Low Dose)")
  print(p10)
  ggsave("Figures_Final/fig_10_placebo_v5.pdf", p10, width = 9, height = 5.5, dpi = 300)
}

# Figure 11: Mid dose bin (NEW v5)
if (!is.null(res_mid)) {
  p11 <- plot_event_study(res_low$tidy, res_mid$tidy,
                          "Low Dose (p20-p50)", "Mid Dose (p50-p75)",
                          "#E69F00", "#F0E442",
                          "Dose Gradient: Low vs. Mid Exposure (NEW v5)")
  print(p11)
  ggsave("Figures_Final/fig_11_mid_dose_v5.pdf", p11, width = 9, height = 5.5, dpi = 300)
}

message("   Main publication figures saved.")

# ==============================================================================
# SECTION 10B: COHORT FIGURES (Duflo-style)
# ==============================================================================
message(" Generating Cohort Figures...")

# Headline: wage gap Low vs High
if (!is.null(res_gap_low) && !is.null(res_gap_high)) {
  p_gap <- plot_event_study(res_gap_low$tidy, res_gap_high$tidy,
                            "Low Dose", "High Dose", "#E69F00", "#0072B2",
                            "Cohort Wage Gap: Young (22-28) vs. Old (29-35)",
                            ylab = "ATT on wage gap (young - old)")
  print(p_gap)
  ggsave("Figures_Final/fig_cohort_gap_v5.pdf", p_gap, width = 9, height = 5.5, dpi = 300)
}

# Young vs Old within Low dose (placebo contrast)
if (!is.null(res_young_low) && !is.null(res_old_low)) {
  p_yo <- plot_event_study(res_young_low$tidy, res_old_low$tidy,
                           "Young (22-28)", "Old (29-35) [Placebo]",
                           "#2166ac", "#999999",
                           "Age-Specific Effects: Exposed vs. Non-Exposed (Low Dose)",
                           ylab = "ATT on log(wages)")
  print(p_yo)
  ggsave("Figures_Final/fig_cohort_young_vs_old_v5.pdf", p_yo, width = 9, height = 5.5, dpi = 300)
}

# Triple interaction: young + low-dose + poor
if (!is.null(res_triple_young_poor) && !is.null(res_young_low)) {
  p_tri <- plot_event_study(res_young_low$tidy, res_triple_young_poor$tidy,
                            "Young / Low Dose", "Young / Low Dose / Poor",
                            "#2166ac", "#d73027",
                            "Triple Interaction: Credit-Constrained Young Workers (Fix 4)",
                            ylab = "ATT on log(wages, young)")
  print(p_tri)
  ggsave("Figures_Final/fig_triple_young_poor_v5.pdf", p_tri, width = 9, height = 5.5, dpi = 300)
}

# College share gap (first stage)
if (!is.null(res_college_low) && !is.null(res_college_high)) {
  p_fs <- plot_event_study(res_college_low$tidy, res_college_high$tidy,
                           "Low Dose", "High Dose", "#E69F00", "#0072B2",
                           "First Stage: College Share Gap (Young - Old)",
                           ylab = "ATT on college share gap")
  print(p_fs)
  ggsave("Figures_Final/fig_first_stage_v5.pdf", p_fs, width = 9, height = 5.5, dpi = 300)
}

# Gender on wage gap
if (!is.null(res_gap_male_low) && !is.null(res_gap_female_low)) {
  p_gm <- plot_event_study(res_gap_male_low$tidy, res_gap_female_low$tidy,
                           "Male", "Female", "#009E73", "#D55E00",
                           "Gender Gap in Cohort Wage Premium: Low Dose",
                           ylab = "ATT on wage gap (young - old)")
  print(p_gm)
  ggsave("Figures_Final/fig_cohort_gender_low_v5.pdf", p_gm, width = 9, height = 5.5, dpi = 300)
}

# Race on wage gap
if (!is.null(res_gap_nonwhite_low) && !is.null(res_gap_white_low)) {
  p_race_low <- plot_event_study(res_gap_nonwhite_low$tidy, res_gap_white_low$tidy,
                                 "Non-white", "White", "#CC79A7", "#56B4E9",
                                 "Ethnicity Gap in Cohort Wage Premium: Low Dose",
                                 ylab = "ATT on wage gap (young - old)")
  print(p_race_low)
  ggsave("Figures_Final/fig_cohort_race_low_v5.pdf", p_race_low, width = 9, height = 5.5, dpi = 300)
}

# Birth-cohort figures
if (!is.null(res_bc_gap_low) && !is.null(res_bc_gap_high)) {
  p_bc <- plot_event_study(res_bc_gap_low$tidy, res_bc_gap_high$tidy,
                           "Low Dose", "High Dose", "#E69F00", "#0072B2",
                           "Birth-Cohort Wage Gap: Exposed (1987-93) vs Control (1975-82)",
                           ylab = "ATT on wage gap (exposed - control)")
  print(p_bc)
  ggsave("Figures_Final/fig_birthcohort_gap_v5.pdf", p_bc, width = 9, height = 5.5, dpi = 300)
}

message("   Cohort figures saved.")

# ==============================================================================
# SECTION 10C: EVENT-STUDY CSV EXPORTS (appendix data)
# ==============================================================================
message(" Exporting event-study CSVs...")

es_exports_v5 <- Filter(Negate(is.null), list(
  es_low_wage       = if (!is.null(res_low))           res_low$tidy           else NULL,
  es_high_wage      = if (!is.null(res_high))          res_high$tidy          else NULL,
  es_mid_wage       = if (!is.null(res_mid))           res_mid$tidy           else NULL,
  es_low_emp        = if (!is.null(res_low_emp))       res_low_emp$tidy       else NULL,
  es_high_emp       = if (!is.null(res_high_emp))      res_high_emp$tidy      else NULL,
  es_male_low       = if (!is.null(res_male_low))      res_male_low$tidy      else NULL,
  es_female_low     = if (!is.null(res_female_low))    res_female_low$tidy    else NULL,
  es_male_high      = if (!is.null(res_male_high))     res_male_high$tidy     else NULL,
  es_female_high    = if (!is.null(res_female_high))   res_female_high$tidy   else NULL,
  es_white_low      = if (!is.null(res_white_low))     res_white_low$tidy     else NULL,
  es_nonwhite_low   = if (!is.null(res_nonwhite_low))  res_nonwhite_low$tidy  else NULL,
  es_white_high     = if (!is.null(res_white_high))    res_white_high$tidy    else NULL,
  es_nonwhite_high  = if (!is.null(res_nonwhite_high)) res_nonwhite_high$tidy else NULL,
  es_low_lag3       = if (!is.null(res_low_lag3))      res_low_lag3$tidy      else NULL,
  es_high_lag3      = if (!is.null(res_high_lag3))     res_high_lag3$tidy     else NULL,
  es_low_dr         = if (!is.null(res_low_dr))        res_low_dr$tidy        else NULL,
  es_high_dr        = if (!is.null(res_high_dr))       res_high_dr$tidy       else NULL,
  es_low_dr_censo   = if (!is.null(res_low_dr_censo))  res_low_dr_censo$tidy  else NULL,
  es_high_dr_censo  = if (!is.null(res_high_dr_censo)) res_high_dr_censo$tidy else NULL,
  es_low_buffer     = if (!is.null(res_low_buffer))    res_low_buffer$tidy    else NULL,
  es_placebo_low    = if (!is.null(res_placebo_low))   res_placebo_low$tidy   else NULL,
  es_gap_low        = if (!is.null(res_gap_low))       res_gap_low$tidy       else NULL,
  es_gap_high       = if (!is.null(res_gap_high))      res_gap_high$tidy      else NULL,
  es_young_low      = if (!is.null(res_young_low))     res_young_low$tidy     else NULL,
  es_old_low        = if (!is.null(res_old_low))       res_old_low$tidy       else NULL,
  es_triple_poor    = if (!is.null(res_triple_young_poor)) res_triple_young_poor$tidy else NULL,
  es_college_low    = if (!is.null(res_college_low))   res_college_low$tidy   else NULL,
  es_bc_gap_low     = if (!is.null(res_bc_gap_low))    res_bc_gap_low$tidy    else NULL,
  es_bc_exposed_low = if (!is.null(res_bc_exposed_low)) res_bc_exposed_low$tidy else NULL,
  es_bc_control_low = if (!is.null(res_bc_control_low)) res_bc_control_low$tidy else NULL
))

for (nm in names(es_exports_v5)) {
  write_csv(es_exports_v5[[nm]],
            sprintf("Tables_Final/es_%s_v5.csv", nm))
}
message(sprintf("   %d event-study CSVs exported.", length(es_exports_v5)))



# ==============================================================================
# SECTION 9T: WILD (MULTIPLIER) BOOTSTRAP -- cluster-robust inference
# ==============================================================================
# Re-runs the two headline wage models with multiplier bootstrap SEs.
# This allows for arbitrary cross-unit correlation in the error process
# and is more conservative than the analytical SEs when spatial spillovers
# are present.
# bstrap = TRUE, boot_type = "multiplier" activates the wild bootstrap
# in the did package (Callaway & Sant'Anna 2021, Section 4).
# ==============================================================================
message(" SECTION 9T: Wild Bootstrap SEs (multiplier)...")

res_low_boot  <- NULL
res_high_boot <- NULL

tryCatch({
  message("   -> Low Dose / Wild Bootstrap...")
  # Build pre-aggregated data (same as run_cs() internals, wmean)
  d_low_agg <- data_low %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop")

  att_low_boot <- did::att_gt(
    yname = "y", tname = "ano", idname = "id_num", gname = "g",
    data  = d_low_agg,
    control_group = "notyettreated",
    bstrap = TRUE, boot_type = "multiplier", nboot = 999,
    allow_unbalanced_panel = TRUE, print_details = FALSE
  )
  es_low_boot    <- did::aggte(att_low_boot, type = "dynamic", na.rm = TRUE)
  att_low_simple <- did::aggte(att_low_boot, type = "simple",  na.rm = TRUE)
  pretrend_low_boot <- pretrend_ftest(att_low_boot)

  res_low_boot <- list(
    att_gt = att_low_boot, es = es_low_boot,
    tidy = broom::tidy(es_low_boot), simple = att_low_simple,
    pretrend = pretrend_low_boot, label = "Low/Boot",
    n_obs = nrow(d_low_agg),
    n_treated = n_distinct(d_low_agg$id_num[d_low_agg$g > 0]),
    n_notyet  = n_distinct(d_low_agg$id_num[d_low_agg$g == 0])
  )
  message(sprintf("   Low/Boot: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                  att_low_simple$overall.att, att_low_simple$overall.se,
                  pretrend_low_boot$p.value))

}, error = function(e) message(sprintf("   WARNING: Low bootstrap failed: %s", safe_err_msg(e))))

tryCatch({
  message("   -> High Dose / Wild Bootstrap...")
  d_high_agg <- data_high %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop")

  att_high_boot <- did::att_gt(
    yname = "y", tname = "ano", idname = "id_num", gname = "g",
    data  = d_high_agg,
    control_group = "notyettreated",
    bstrap = TRUE, boot_type = "multiplier", nboot = 999,
    allow_unbalanced_panel = TRUE, print_details = FALSE
  )
  es_high_boot    <- did::aggte(att_high_boot, type = "dynamic", na.rm = TRUE)
  att_high_simple <- did::aggte(att_high_boot, type = "simple",  na.rm = TRUE)
  pretrend_high_boot <- pretrend_ftest(att_high_boot)

  res_high_boot <- list(
    att_gt = att_high_boot, es = es_high_boot,
    tidy = broom::tidy(es_high_boot), simple = att_high_simple,
    pretrend = pretrend_high_boot, label = "High/Boot",
    n_obs = nrow(d_high_agg),
    n_treated = n_distinct(d_high_agg$id_num[d_high_agg$g > 0]),
    n_notyet  = n_distinct(d_high_agg$id_num[d_high_agg$g == 0])
  )
  message(sprintf("   High/Boot: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                  att_high_simple$overall.att, att_high_simple$overall.se,
                  pretrend_high_boot$p.value))

}, error = function(e) message(sprintf("   WARNING: High bootstrap failed: %s", safe_err_msg(e))))

# Compare analytical vs bootstrap SEs
if (!is.null(res_low_boot)) {
  message(sprintf("   Bootstrap comparison Low:  Analytical SE=%.4f | Bootstrap SE=%.4f",
                  res_low$simple$overall.se, res_low_boot$simple$overall.se))
}
if (!is.null(res_high_boot)) {
  message(sprintf("   Bootstrap comparison High: Analytical SE=%.4f | Bootstrap SE=%.4f",
                  res_high$simple$overall.se, res_high_boot$simple$overall.se))
}

# Bootstrap event-study figures
if (!is.null(res_low_boot) && !is.null(res_high_boot)) {
  p_boot <- plot_event_study(
    res_low_boot$tidy, res_high_boot$tidy,
    "Low Dose (bootstrap)", "High Dose (bootstrap)",
    "#E69F00", "#0072B2",
    "Dose-Response: Wild Bootstrap SEs (999 replications)",
    ylab = "ATT on log(wages)"
  )
  ggsave("Figures_Final/fig_bootstrap_v5.pdf", p_boot, width = 9, height = 5.5, dpi = 300)
  message("   Bootstrap figure saved.")
}



# ==============================================================================
# SECTION 9U: POPULATION-WEIGHTED CS ESTIMATOR
# ==============================================================================
# Weight each microregion by its 2010 youth population (pop_18_24) so that
# the ATT reflects the average effect for a randomly drawn young worker,
# rather than for a randomly drawn microregion. Larger microregions (e.g.
# Sao Paulo, Rio) contribute proportionally more to the aggregate estimate.
#
# Implemented via run_cs() with weightsname by augmenting the aggregated
# data with per-microregion pop_18_24 weights.
# ==============================================================================
message(" SECTION 9U: Population-Weighted CS Estimator...")

res_low_wt  <- NULL
res_high_wt <- NULL

tryCatch({
  # Build population weight map (one row per microregion)
  pop_weights <- df_pop_micro %>%
    mutate(id_num = as.integer(id_microrregiao)) %>%
    select(id_num, pop_wt = pop)

  # Low Dose -- population-weighted
  message("   -> Low Dose / Population-Weighted...")
  d_low_wt <- data_low %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(pop_weights, by = "id_num") %>%
    mutate(pop_wt = replace_na(pop_wt, 1))

  att_low_wt <- did::att_gt(
    yname = "y", tname = "ano", idname = "id_num", gname = "g",
    weightsname = "pop_wt",
    data  = d_low_wt,
    control_group = "notyettreated",
    allow_unbalanced_panel = TRUE, print_details = FALSE
  )
  es_low_wt    <- did::aggte(att_low_wt, type = "dynamic", na.rm = TRUE)
  att_low_wt_s <- did::aggte(att_low_wt, type = "simple",  na.rm = TRUE)
  pretrend_low_wt <- pretrend_ftest(att_low_wt)

  res_low_wt <- list(
    att_gt = att_low_wt, es = es_low_wt,
    tidy = broom::tidy(es_low_wt), simple = att_low_wt_s,
    pretrend = pretrend_low_wt, label = "Low/PopWt",
    n_obs = nrow(d_low_wt),
    n_treated = n_distinct(d_low_wt$id_num[d_low_wt$g > 0]),
    n_notyet  = n_distinct(d_low_wt$id_num[d_low_wt$g == 0])
  )
  message(sprintf("   Low/PopWt: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                  att_low_wt_s$overall.att, att_low_wt_s$overall.se,
                  pretrend_low_wt$p.value))

}, error = function(e) message(sprintf("   WARNING: Low/PopWt failed: %s", safe_err_msg(e))))

tryCatch({
  message("   -> High Dose / Population-Weighted...")
  d_high_wt <- data_high %>%
    group_by(id_num, ano, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(pop_weights, by = "id_num") %>%
    mutate(pop_wt = replace_na(pop_wt, 1))

  att_high_wt <- did::att_gt(
    yname = "y", tname = "ano", idname = "id_num", gname = "g",
    weightsname = "pop_wt",
    data  = d_high_wt,
    control_group = "notyettreated",
    allow_unbalanced_panel = TRUE, print_details = FALSE
  )
  es_high_wt    <- did::aggte(att_high_wt, type = "dynamic", na.rm = TRUE)
  att_high_wt_s <- did::aggte(att_high_wt, type = "simple",  na.rm = TRUE)
  pretrend_high_wt <- pretrend_ftest(att_high_wt)

  res_high_wt <- list(
    att_gt = att_high_wt, es = es_high_wt,
    tidy = broom::tidy(es_high_wt), simple = att_high_wt_s,
    pretrend = pretrend_high_wt, label = "High/PopWt",
    n_obs = nrow(d_high_wt),
    n_treated = n_distinct(d_high_wt$id_num[d_high_wt$g > 0]),
    n_notyet  = n_distinct(d_high_wt$id_num[d_high_wt$g == 0])
  )
  message(sprintf("   High/PopWt: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                  att_high_wt_s$overall.att, att_high_wt_s$overall.se,
                  pretrend_high_wt$p.value))

}, error = function(e) message(sprintf("   WARNING: High/PopWt failed: %s", safe_err_msg(e))))

# Population-weighted vs unweighted comparison figure
if (!is.null(res_low_wt) && !is.null(res_high_wt)) {
  p_wt <- plot_event_study(
    res_low_wt$tidy, res_high_wt$tidy,
    "Low Dose (pop-weighted)", "High Dose (pop-weighted)",
    "#E69F00", "#0072B2",
    "Population-Weighted CS: ATT for a Randomly Drawn Worker"
  )
  ggsave("Figures_Final/fig_popwt_v5.pdf", p_wt, width = 9, height = 5.5, dpi = 300)
  message("   Population-weighted figure saved.")
}

# ==============================================================================
# SECTION 9V: BIENNIAL PANEL (2-year collapsed windows)
# ==============================================================================
# Collapse the annual panel into 2-year periods (2005-06, 2007-08, ..., 2017-18)
# to reduce within-unit serial autocorrelation and improve precision of
# long-run estimates. Treatment timing g_biennial = ceiling(g / 2) * 2.
# This follows the "period-collapsing" approach discussed in Borusyak et al.
# (2024) as a check on temporal autocorrelation in panel SEs.
# ==============================================================================
message(" SECTION 9V: Biennial Panel Robustness...")

res_low_biennal  <- NULL
res_high_biennal <- NULL

tryCatch({
  # Create biennial year variable: 2005->2006, 2006->2006, 2007->2008, ...
  biennial_year <- function(y) as.integer(ceiling(y / 2) * 2)

  d_biennial_low <- data_low %>%
    mutate(ano_bi = biennial_year(ano)) %>%
    group_by(id_num, ano_bi, g) %>%
    summarise(
      y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Recode g to biennial: earliest biennial period >= g
    mutate(
      g_bi = ifelse(g == 0, 0L, biennial_year(g))
    ) %>%
    rename(ano = ano_bi, g = g_bi) %>%
    # Keep only distinct (id_num, ano) -- average within 2-year window
    group_by(id_num, ano, g) %>%
    summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

  n_treated_bi <- n_distinct(d_biennial_low$id_num[d_biennial_low$g > 0])
  n_notyet_bi  <- n_distinct(d_biennial_low$id_num[d_biennial_low$g == 0])
  message(sprintf("   Biennial Low: %d obs | %d treated | %d not-yet",
                  nrow(d_biennial_low), n_treated_bi, n_notyet_bi))

  if (n_treated_bi >= 20 && n_notyet_bi >= 20) {
    att_bi_low <- did::att_gt(
      yname = "y", tname = "ano", idname = "id_num", gname = "g",
      data  = d_biennial_low,
      control_group = "notyettreated",
      allow_unbalanced_panel = TRUE, print_details = FALSE
    )
    es_bi_low    <- did::aggte(att_bi_low, type = "dynamic", na.rm = TRUE)
    att_bi_low_s <- did::aggte(att_bi_low, type = "simple",  na.rm = TRUE)
    pretrend_bi_low <- pretrend_ftest(att_bi_low)

    res_low_biennal <- list(
      att_gt = att_bi_low, es = es_bi_low,
      tidy = broom::tidy(es_bi_low), simple = att_bi_low_s,
      pretrend = pretrend_bi_low, label = "Low/Biennial",
      n_obs = nrow(d_biennial_low), n_treated = n_treated_bi, n_notyet = n_notyet_bi
    )
    message(sprintf("   Low/Biennial: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                    att_bi_low_s$overall.att, att_bi_low_s$overall.se,
                    pretrend_bi_low$p.value))
  }

  # High dose biennial
  d_biennial_high <- data_high %>%
    mutate(ano_bi = biennial_year(ano)) %>%
    group_by(id_num, ano_bi, g) %>%
    summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(g_bi = ifelse(g == 0, 0L, biennial_year(g))) %>%
    rename(ano = ano_bi, g = g_bi) %>%
    group_by(id_num, ano, g) %>%
    summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

  n_treated_bih <- n_distinct(d_biennial_high$id_num[d_biennial_high$g > 0])
  n_notyet_bih  <- n_distinct(d_biennial_high$id_num[d_biennial_high$g == 0])

  if (n_treated_bih >= 20 && n_notyet_bih >= 20) {
    att_bi_high <- did::att_gt(
      yname = "y", tname = "ano", idname = "id_num", gname = "g",
      data  = d_biennial_high,
      control_group = "notyettreated",
      allow_unbalanced_panel = TRUE, print_details = FALSE
    )
    es_bi_high    <- did::aggte(att_bi_high, type = "dynamic", na.rm = TRUE)
    att_bi_high_s <- did::aggte(att_bi_high, type = "simple",  na.rm = TRUE)
    pretrend_bi_high <- pretrend_ftest(att_bi_high)

    res_high_biennal <- list(
      att_gt = att_bi_high, es = es_bi_high,
      tidy = broom::tidy(es_bi_high), simple = att_bi_high_s,
      pretrend = pretrend_bi_high, label = "High/Biennial",
      n_obs = nrow(d_biennial_high), n_treated = n_treated_bih, n_notyet = n_notyet_bih
    )
    message(sprintf("   High/Biennial: ATT = %.4f (SE = %.4f) | Pre-trend p = %.4f",
                    att_bi_high_s$overall.att, att_bi_high_s$overall.se,
                    pretrend_bi_high$p.value))
  }

  # Biennial robustness figure
  if (!is.null(res_low_biennal) && !is.null(res_high_biennal)) {
    p_bi <- plot_event_study(
      res_low_biennal$tidy, res_high_biennal$tidy,
      "Low Dose (biennial)", "High Dose (biennial)",
      "#E69F00", "#0072B2",
      "Biennial Panel Robustness (2-year collapsed windows)"
    )
    ggsave("Figures_Final/fig_biennial_v5.pdf", p_bi, width = 9, height = 5.5, dpi = 300)
    message("   Biennial figure saved.")
  }

}, error = function(e) message(sprintf("   WARNING: Biennial panel failed: %s", safe_err_msg(e))))

# Add bootstrap, pop-weighted, biennial to summary table
extra_rob_results <- Filter(Negate(is.null), list(
  res_low_boot, res_high_boot,
  res_low_wt, res_high_wt,
  res_low_biennal, res_high_biennal
))

if (length(extra_rob_results) > 0) {
  tab_extra <- bind_rows(lapply(extra_rob_results, build_att_row))
  tab_att   <- bind_rows(tab_att, tab_extra)
  write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
  message(sprintf("   ATT table updated: now %d models (added bootstrap/popwt/biennial).",
                  nrow(tab_att)))
}



# ==============================================================================
# SECTION 11: LaTeX TABLE EXPORTS
# ==============================================================================
message(" [11B/12] Exporting LaTeX tables...")

# Table A1: Continuous dose results (fixest)
if (!is.null(cont_results)) {
  print(xtable(cont_results,
               caption = "Continuous Dose-Response: fixest TWFE Specifications",
               label   = "tab:continuous_dose",
               digits  = c(0, 0, 6, 6, 3, 4, 0)),
        file = "Tables_Final/tab_continuous_dose_v5.tex",
        include.rownames = FALSE)
  message("   Table A1 (continuous dose) saved.")
}

# Table A2: Bootstrap vs analytical SE comparison
boot_compare <- bind_rows(
  Filter(Negate(is.null), list(
    if (!is.null(res_low))       build_att_row(res_low)       else NULL,
    if (!is.null(res_low_boot))  build_att_row(res_low_boot)  else NULL,
    if (!is.null(res_high))      build_att_row(res_high)      else NULL,
    if (!is.null(res_high_boot)) build_att_row(res_high_boot) else NULL
  ))
) %>% select(Model, ATT, SE, CI_low, CI_high, Pretrend_p)

if (nrow(boot_compare) > 0) {
  write_csv(boot_compare, "Tables_Final/tab_bootstrap_compare_v5.csv")
  print(xtable(boot_compare,
               caption = "Analytical vs Wild Bootstrap Standard Errors",
               label   = "tab:bootstrap",
               digits  = c(0, 0, 4, 4, 4, 4, 4)),
        file = "Tables_Final/tab_bootstrap_compare_v5.tex",
        include.rownames = FALSE)
  message("   Table A2 (bootstrap comparison) saved.")
}

# Table A3: Population-weighted vs unweighted
wt_compare <- bind_rows(
  Filter(Negate(is.null), list(
    if (!is.null(res_low))    build_att_row(res_low)    else NULL,
    if (!is.null(res_low_wt)) build_att_row(res_low_wt) else NULL,
    if (!is.null(res_high))   build_att_row(res_high)   else NULL,
    if (!is.null(res_high_wt)) build_att_row(res_high_wt) else NULL
  ))
) %>% select(Model, ATT, SE, CI_low, CI_high, Pretrend_p)

if (nrow(wt_compare) > 0) {
  write_csv(wt_compare, "Tables_Final/tab_popwt_compare_v5.csv")
  print(xtable(wt_compare,
               caption = "Unweighted vs Population-Weighted CS Estimates",
               label   = "tab:popwt",
               digits  = c(0, 0, 4, 4, 4, 4, 4)),
        file = "Tables_Final/tab_popwt_compare_v5.tex",
        include.rownames = FALSE)
  message("   Table A3 (population-weighted comparison) saved.")
}

# Table A4: DR comparison (base vs census-enriched)
dr_compare <- bind_rows(
  Filter(Negate(is.null), list(
    if (!is.null(res_low))          build_att_row(res_low)          else NULL,
    if (!is.null(res_low_dr))       build_att_row(res_low_dr)       else NULL,
    if (!is.null(res_low_dr_censo)) build_att_row(res_low_dr_censo) else NULL,
    if (!is.null(res_high))         build_att_row(res_high)         else NULL,
    if (!is.null(res_high_dr))      build_att_row(res_high_dr)      else NULL,
    if (!is.null(res_high_dr_censo)) build_att_row(res_high_dr_censo) else NULL
  ))
) %>% select(Model, ATT, SE, CI_low, CI_high, t_val, p_val, Pretrend_p)

if (nrow(dr_compare) > 0) {
  write_csv(dr_compare, "Tables_Final/tab_dr_compare_v5.csv")
  print(xtable(dr_compare,
               caption = paste0("Doubly-Robust Specifications: Unconditional, Base DR,",
                                " and DR Enriched with Education Census Covariates"),
               label   = "tab:dr_compare",
               digits  = c(0, 0, 4, 4, 4, 4, 3, 4, 4)),
        file = "Tables_Final/tab_dr_compare_v5.tex",
        include.rownames = FALSE)
  message("   Table A4 (DR comparison) saved.")
}

# Table A5: Cohort models summary (Duflo-style)
cohort_summary_models <- Filter(Negate(is.null), list(
  res_gap_low, res_gap_high,
  res_young_low, res_young_high,
  res_old_low, res_old_high,
  res_triple_young_poor,
  res_bc_gap_low, res_bc_gap_high,
  res_bc_exposed_low, res_bc_control_low
))

if (length(cohort_summary_models) > 0) {
  tab_cohort <- bind_rows(lapply(cohort_summary_models, build_att_row)) %>%
    select(Model, ATT, SE, t_val, p_val, Pretrend_p, N_obs)
  write_csv(tab_cohort, "Tables_Final/tab_cohort_summary_v5.csv")
  print(xtable(tab_cohort,
               caption = "Cohort-Based Estimates: Duflo (2001) Age-Cohort Design",
               label   = "tab:cohort",
               digits  = c(0, 0, 4, 4, 3, 4, 4, 0)),
        file = "Tables_Final/tab_cohort_summary_v5.tex",
        include.rownames = FALSE)
  message("   Table A5 (cohort summary) saved.")
}

# Table A6: Pre-trend p-values (all models)
tab_pretrend <- tab_att %>%
  select(Model, ATT, SE, t_val, p_val, Pretrend_p) %>%
  filter(!is.na(Pretrend_p)) %>%
  arrange(Pretrend_p)
write_csv(tab_pretrend, "Tables_Final/tab_pretrend_v5.csv")
print(xtable(tab_pretrend,
             caption = "Pre-Treatment Parallel Trends Tests (Wald chi-squared, all models)",
             label   = "tab:pretrend",
             digits  = c(0, 0, 4, 4, 3, 4, 4)),
      file = "Tables_Final/tab_pretrend_v5.tex",
      include.rownames = FALSE)
message("   Table A6 (pre-trend summary) saved.")

# Table A7: Education Census covariate summary
if (exists("df_censo_micro") && nrow(df_censo_micro) > 0) {
  censo_summary <- df_censo_micro %>%
    summarise(
      N_micro              = n(),
      mean_college_density = round(mean(college_density, na.rm = TRUE), 3),
      sd_college_density   = round(sd(college_density,   na.rm = TRUE), 3),
      mean_private_share   = round(mean(private_enroll_share, na.rm = TRUE), 3),
      pct_zero_college     = round(mean(college_density == 0, na.rm = TRUE) * 100, 1)
    )
  write_csv(censo_summary, "Tables_Final/tab_censo_covars_v5.csv")
  message("   Table A7 (Education Census covariates) saved.")
}

message("   All LaTeX tables exported.")



# ==============================================================================
# SECTION 12: COMPREHENSIVE AUDIT LOG
# ==============================================================================
message(" [12/12] Saving Audit Log...")

audit_log_v5 <- list(
  script_version = "untitled5.R (v5)",
  timestamp      = Sys.time(),
  params         = params,

  # Sample
  n_microregions = n_distinct(df_final$id_microrregiao),
  n_years        = length(params$years),
  year_range     = c(min(params$years), max(params$years)),

  # Treatment construction
  bin_method   = "2019 cross-section (D_rt_2019), lagged 4 years",
  thresholds   = list(tau_20 = tau_20, tau_50 = tau_50, tau_75 = tau_75),
  thresholds_lag3 = list(tau_20 = tau_20_alt, tau_50 = tau_50_alt, tau_75 = tau_75_alt),
  bin_sizes    = as.list(setNames(bin_sizes$n_regions, bin_sizes$dose_bin)),
  g_range      = list(
    min = min(df_final$g[df_final$g > 0], na.rm = TRUE),
    max = max(df_final$g[df_final$g > 0], na.rm = TRUE)
  ),

  # v5 fixes applied
  fixes_applied = list(
    fix1_contdid     = "target_parameter passed as scalar string",
    fix2_honestdid   = "exhaustive VCov search + diagonal fallback",
    fix3_dr_covars   = "Low: 4 base + 2 census covars; High: 2+1 covars",
    fix4_triple      = "young x low x below-median income microregions",
    fix5_wald_late   = "dropped (weak first stage, F < 1)",
    fix6_placebo     = "placebo cohort 36-50 activated"
  ),

  # New data sources
  data_sources = list(
    rais        = "br_me_rais.microdados_vinculos (2005-2019)",
    prouni      = "br_mec_prouni.microdados (2005-2019)",
    census      = "br_ibge_censo_demografico.microdados_pessoa_2010",
    censo_educ  = "br_inep_censo_educacao_superior.curso (2009 snapshot)"
  ),

  # Headline ATTs
  att_headline = list(
    low_wages  = list(att = res_low$simple$overall.att,  se = res_low$simple$overall.se),
    high_wages = list(att = res_high$simple$overall.att, se = res_high$simple$overall.se),
    low_dr     = if (!is.null(res_low_dr))
                   list(att = res_low_dr$simple$overall.att, se = res_low_dr$simple$overall.se)
                 else NA,
    low_dr_censo = if (!is.null(res_low_dr_censo))
                     list(att = res_low_dr_censo$simple$overall.att,
                          se  = res_low_dr_censo$simple$overall.se)
                   else NA
  ),

  # Pre-trend p-values (all non-NULL models)
  pretrend_pvals = setNames(
    sapply(all_results_v5, function(r) r$pretrend$p.value),
    sapply(all_results_v5, function(r) r$label)
  ),

  # HonestDiD breakdown Mbar
  honest_did_breakdown = lapply(honest_results, function(x)
                                  list(breakdown = x$breakdown,
                                       used_full_vcov = x$used_full_vcov)),

  # Robustness checks
  bootstrap = list(
    low_analytical = res_low$simple$overall.se,
    low_bootstrap  = if (!is.null(res_low_boot))  res_low_boot$simple$overall.se  else NA,
    high_analytical= res_high$simple$overall.se,
    high_bootstrap = if (!is.null(res_high_boot)) res_high_boot$simple$overall.se else NA
  ),

  pop_weighted = list(
    low  = if (!is.null(res_low_wt))  list(att = res_low_wt$simple$overall.att,  se = res_low_wt$simple$overall.se)  else NA,
    high = if (!is.null(res_high_wt)) list(att = res_high_wt$simple$overall.att, se = res_high_wt$simple$overall.se) else NA
  ),

  biennial = list(
    low  = if (!is.null(res_low_biennal))  list(att = res_low_biennal$simple$overall.att,  se = res_low_biennal$simple$overall.se)  else NA,
    high = if (!is.null(res_high_biennal)) list(att = res_high_biennal$simple$overall.att, se = res_high_biennal$simple$overall.se) else NA
  ),

  # TWFE bias benchmark
  twfe = list(
    wages = if (!is.null(twfe_fit)) {
      tc <- tryCatch(fixest::coeftable(twfe_fit), error = function(e) NULL)
      if (!is.null(tc)) list(coef = tc[1,1], se = tc[1,2], p = tc[1,4]) else NA
    } else NA,
    employment = if (!is.null(twfe_emp)) {
      tc <- tryCatch(fixest::coeftable(twfe_emp), error = function(e) NULL)
      if (!is.null(tc)) list(coef = tc[1,1], se = tc[1,2], p = tc[1,4]) else NA
    } else NA
  ),

  # Contdid
  contdid_available = list(
    level = !is.null(res_contdid_level),
    slope = !is.null(res_contdid_slope)
  ),

  # Buffer zone
  buffer_zone = list(
    n_excluded  = length(border_ids),
    low_att     = if (!is.null(res_low_buffer)) res_low_buffer$simple$overall.att else NA,
    high_att    = if (!is.null(res_high_buffer)) res_high_buffer$simple$overall.att else NA
  ),

  # Case studies
  case_studies = case_micro,

  # Full ATT table
  att_summary = tab_att,

  # Continuous dose table
  cont_dose_results = cont_results,

  # Alt cohort cutoffs
  alt_cohort_pretrends = lapply(alt_cohort_results, function(r)
                                  list(label = r$label, att = r$simple$overall.att,
                                       pretrend_p = r$pretrend$p.value)),

  dr_divergence_flag = dr_divergence_flag
)

saveRDS(audit_log_v5, "Documentation_Final/audit_log_v5.rds")
message("   Audit log v5 saved.")

# Also write a human-readable audit summary as plain text
audit_txt <- c(
  "================================================================",
  "  UNTITLED5.R -- v5 AUDIT SUMMARY",
  sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "================================================================",
  sprintf("  Microregions : %d", as.integer(audit_log_v5$n_microregions)),
  sprintf("  Years        : %d-%d",
          as.integer(audit_log_v5$year_range[1]),
          as.integer(audit_log_v5$year_range[2])),
  sprintf("  tau_20=%.4f  tau_50=%.4f  tau_75=%.4f", tau_20, tau_50, tau_75),
  sprintf("  Bin sizes: Low=%s  High=%s",
          if ("Low"  %in% names(audit_log_v5$bin_sizes)) audit_log_v5$bin_sizes$Low  else "NA",
          if ("High" %in% names(audit_log_v5$bin_sizes)) audit_log_v5$bin_sizes$High else "NA"),
  sprintf("  Treatment g range: %d-%d",
          as.integer(audit_log_v5$g_range$min),
          as.integer(audit_log_v5$g_range$max)),
  "  --- Headline ATTs ---",
  sprintf("  Low/Wages : %+.4f (SE %.4f)",
          res_low$simple$overall.att, res_low$simple$overall.se),
  sprintf("  High/Wages: %+.4f (SE %.4f)",
          res_high$simple$overall.att, res_high$simple$overall.se),
  if (!is.null(res_low_dr))
    sprintf("  Low/DR    : %+.4f (SE %.4f)",
            res_low_dr$simple$overall.att, res_low_dr$simple$overall.se)
  else
    "  Low/DR    : skipped",
  if (!is.null(res_low_dr_censo))
    sprintf("  Low/DR/Censo: %+.4f (SE %.4f)",
            res_low_dr_censo$simple$overall.att, res_low_dr_censo$simple$overall.se)
  else
    "  Low/DR/Censo: skipped",
  sprintf("  Total models estimated: %d", as.integer(nrow(tab_att))),
  sprintf("  All pre-trend p > 0.10: %s",
          all(tab_att$Pretrend_p > 0.10, na.rm = TRUE)),
  "================================================================"
)
writeLines(audit_txt, "Documentation_Final/audit_summary_v5.txt")
message("   Audit summary text saved.")



# ==============================================================================
# SECTION 12B: print_all_results() -- call at any time after script completes
# ==============================================================================
print_all_results <- function() {

  divider <- function(title) {
    cat("\n", strrep("=", 72), "\n  ", title, "\n", strrep("=", 72), "\n", sep = "")
  }

  divider("SAMPLE & TREATMENT CONSTRUCTION")
  cat(sprintf("  Microregions : %d | Years : %d-%d\n",
              as.integer(n_distinct(df_final$id_microrregiao)),
              as.integer(min(params$years)),
              as.integer(max(params$years))))
  cat(sprintf("  tau_20=%.4f  tau_50=%.4f  tau_75=%.4f\n", tau_20, tau_50, tau_75))
  print(bin_sizes)

  divider("TABLE 1: BIN SUMMARY")
  print(tab1, n = Inf)

  divider("TABLE 2: COVARIATE BALANCE")
  print(tab2, n = Inf)

  divider("TABLE 3: CASE STUDY MICROREGIONS")
  print(case_micro, n = Inf)

  divider("PRIMARY SPECIFICATION -- Doubly-Robust Wages")
  for (r in Filter(Negate(is.null),
                   list(res_low, res_low_dr, res_low_dr_censo,
                        res_high, res_high_dr, res_high_dr_censo))) {
    p <- 2 * pnorm(-abs(r$simple$overall.att / r$simple$overall.se))
    cat(sprintf("  %-22s ATT = %+.4f  SE = %.4f  t = %+.3f  p = %.3f\n",
                r$label, r$simple$overall.att, r$simple$overall.se,
                r$simple$overall.att / r$simple$overall.se, p))
  }

  divider("CONTINUOUS DOSE (fixest)")
  if (!is.null(cont_results)) print(as.data.frame(cont_results))
  for (nm_fit in c("fit_young_log", "fit_bc_log", "fit_bc_quad")) {
    obj <- tryCatch(get(nm_fit), error = function(e) NULL)
    if (!is.null(obj)) {
      ct <- coeftable(obj)
      cat(sprintf("  %-20s  coef=%.6f  t=%.3f  p=%.4f\n",
                  nm_fit, ct[1,1], ct[1,3], ct[1,4]))
    }
  }

  divider("ROBUSTNESS CHECKS")
  for (r in Filter(Negate(is.null),
                   list(res_low_lag3, res_high_lag3,
                        res_low_25_75, res_high_25_75,
                        res_low_33_67, res_high_33_67,
                        res_low_buffer, res_high_buffer,
                        res_low_boot, res_high_boot,
                        res_low_wt, res_high_wt,
                        res_low_biennal, res_high_biennal))) {
    cat(sprintf("  %-22s ATT = %+.4f  SE = %.4f  Pre-trend p = %.4f\n",
                r$label, r$simple$overall.att, r$simple$overall.se,
                r$pretrend$p.value))
  }

  divider("COHORT MODELS (Duflo-style)")
  for (r in Filter(Negate(is.null),
                   list(res_gap_low, res_gap_high,
                        res_young_low, res_young_high,
                        res_old_low, res_old_high,
                        res_triple_young_poor,
                        res_college_low, res_college_high,
                        res_gap_male_low, res_gap_male_high,
                        res_gap_female_low, res_gap_female_high,
                        res_gap_white_low, res_gap_white_high,
                        res_gap_nonwhite_low, res_gap_nonwhite_high))) {
    cat(sprintf("  %-28s ATT = %+.4f  SE = %.4f  Pre-trend p = %.4f\n",
                r$label, r$simple$overall.att, r$simple$overall.se,
                r$pretrend$p.value))
  }

  divider("BIRTH-COHORT (Duflo-strict: 1987-93 vs 1975-82)")
  for (r in Filter(Negate(is.null),
                   list(res_bc_gap_low, res_bc_gap_high,
                        res_bc_exposed_low, res_bc_control_low))) {
    cat(sprintf("  %-22s ATT = %+.4f  SE = %.4f  Pre-trend p = %.4f\n",
                r$label, r$simple$overall.att, r$simple$overall.se,
                r$pretrend$p.value))
  }

  divider("PLACEBO COHORT (ages 36-50) -- Fix 6")
  for (r in Filter(Negate(is.null), list(res_placebo_low, res_placebo_high))) {
    cat(sprintf("  %-22s ATT = %+.4f  SE = %.4f  Pre-trend p = %.4f\n",
                r$label, r$simple$overall.att, r$simple$overall.se,
                r$pretrend$p.value))
  }

  divider("TWFE BIAS BENCHMARK")
  for (tw in list(list(fit = twfe_fit, lbl = "TWFE/Wages"),
                  list(fit = twfe_emp, lbl = "TWFE/Employment"))) {
    if (!is.null(tw$fit)) {
      tc <- tryCatch(fixest::coeftable(tw$fit), error = function(e) NULL)
      if (!is.null(tc))
        cat(sprintf("  %-22s coef = %+.4f  SE = %.4f  p = %.4f\n",
                    tw$lbl, tc[1,1], tc[1,2], tc[1,4]))
    }
  }

  divider("HONESTDID BREAKDOWN Mbar")
  for (nm in names(honest_results)) {
    hr <- honest_results[[nm]]
    cat(sprintf("  %-25s Mbar = %.2f  (VCov: %s)\n",
                nm, hr$breakdown,
                if (isTRUE(hr$used_full_vcov)) "full matrix" else "diagonal fallback"))
  }

  divider("TABLE 4: ATT SUMMARY (ALL MODELS)")
  print(tab_att, n = Inf)

  divider("PRE-TREND P-VALUES (all models, sorted)")
  pt <- tab_att %>% select(Model, ATT, SE, p_val, Pretrend_p) %>%
    filter(!is.na(Pretrend_p)) %>% arrange(Pretrend_p)
  print(as.data.frame(pt), row.names = FALSE)

  divider("OUTPUT FILES")
  for (d in c("Tables_Final", "Figures_Final", "Documentation_Final")) {
    if (dir.exists(d)) {
      files <- list.files(d)
      cat(sprintf("\n  %s/ (%d files)\n", d, length(files)))
      cat(paste0("    ", files, "\n"), sep = "")
    }
  }

  cat("\n", strrep("=", 72), "\n", sep = "")
  invisible(NULL)
}

# ==============================================================================
# v5 AUDIT SUMMARY BLOCK (auto-prints at end of script execution)
# ==============================================================================
message("\n\n")
message(strrep("=", 72))
message("  UNTITLED5.R -- v5 EXECUTION COMPLETE")
message(strrep("=", 72))
message(sprintf("  Timestamp     : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
message(sprintf("  Microregions  : %d | Years : %d-%d",
                n_distinct(df_final$id_microrregiao),
                min(params$years), max(params$years)))
message(sprintf("  Thresholds    : tau_20=%.4f  tau_50=%.4f  tau_75=%.4f",
                tau_20, tau_50, tau_75))
message(sprintf("  Models run    : %d total", nrow(tab_att)))
message(sprintf("  Pre-trends OK : %s",
                if (all(tab_att$Pretrend_p > 0.10, na.rm = TRUE)) "YES (all p > 0.10)"
                else sprintf("WARN -- %d model(s) with p <= 0.10",
                             sum(tab_att$Pretrend_p <= 0.10, na.rm = TRUE))))
message("  --- Headline ---")
message(sprintf("  Low/Wages     : ATT = %+.4f (SE = %.4f, p = %.3f)",
                res_low$simple$overall.att, res_low$simple$overall.se,
                2*pnorm(-abs(res_low$simple$overall.att/res_low$simple$overall.se))))
message(sprintf("  High/Wages    : ATT = %+.4f (SE = %.4f, p = %.3f)",
                res_high$simple$overall.att, res_high$simple$overall.se,
                2*pnorm(-abs(res_high$simple$overall.att/res_high$simple$overall.se))))
if (!is.null(res_low_dr_censo)) {
  message(sprintf("  Low/DR/Censo  : ATT = %+.4f (SE = %.4f, p = %.3f)",
                  res_low_dr_censo$simple$overall.att,
                  res_low_dr_censo$simple$overall.se,
                  2*pnorm(-abs(res_low_dr_censo$simple$overall.att /
                               res_low_dr_censo$simple$overall.se))))
}
message("  --- v5 Fixes Applied ---")
message("  Fix 1: contdid  -- scalar target_parameter ('level' + 'slope')")
message("  Fix 2: HonestDiD -- exhaustive VCov search + diagonal fallback")
message("  Fix 3: DR covars -- 4/6 (Low) | 2/3 (High) + Education Census")
message("  Fix 4: Triple   -- young x low-dose x below-median income")
message("  Fix 5: Wald LATE -- dropped (weak first stage F < 1)")
message("  Fix 6: Placebo  -- cohort 36-50 activated")
message("  --- New additions ---")
message("  Education Census (INEP 2009) integrated as DR covariate")
message("  Wild bootstrap SEs (999 replications)")
message("  Population-weighted CS estimator")
message("  Biennial panel robustness check")
message("  Mid dose bin (p50-p75) added")
message("  Type: call print_all_results() for full output")
message(strrep("=", 72))
message("DONE. All outputs saved to Figures_Final/ Tables_Final/ Documentation_Final/")



# ==============================================================================
# SECTION 3I: IBGE PIB -- MUNICIPAL GDP (time-varying, 2005-2019)
# ==============================================================================
message(" SECTION 3I: Downloading IBGE PIB (Municipal GDP 2005-2019)...")

df_pib_raw <- load_or_fetch(
  "Documentation_Final/cache_pib.rds",
  function() {
    basedosdados::read_sql("
      SELECT
        ano,
        id_municipio,
        pib                AS pib_nominal,
        impostos_liquidos  AS impostos,
        va                 AS valor_adicionado,
        va_agropecuaria,
        va_industria,
        va_servicos,
        va_adespam         AS va_admin_publica
      FROM `basedosdados.br_ibge_pib.municipio`
      WHERE ano BETWEEN 2005 AND 2019
    ")
  },
  "ibge_pib_2005_2019"
)
message(sprintf("   PIB raw: %d rows", nrow(df_pib_raw)))

# Annual population for per-capita GDP
df_pop_annual <- load_or_fetch(
  "Documentation_Final/cache_pop_annual.rds",
  function() {
    basedosdados::read_sql("
      SELECT
        ano,
        id_municipio,
        populacao
      FROM `basedosdados.br_ibge_populacao.municipio`
      WHERE ano BETWEEN 2005 AND 2019
    ")
  },
  "ibge_populacao_2005_2019"
)
message(sprintf("   Population annual: %d rows", nrow(df_pop_annual)))
gc()

# ==============================================================================
# SECTION 3K: INEP CENSO EDUCACAO SUPERIOR -- HEI COUNTS (IES level)
# ==============================================================================
# Note: basedosdados only has 2009+. We use 2009 as the earliest pre-treatment
# proxy. The IES table has one row per institution per year.
# ==============================================================================
message(" SECTION 3K: Downloading INEP HEI counts (IES-level)...")

df_hei_raw <- load_or_fetch(
  "Documentation_Final/cache_hei_ies.rds",
  function() {
    # Primary: IES-level table
    tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_municipio,
          id_ies,
          tipo_categoria_administrativa,
          tipo_organizacao_academica
        FROM `basedosdados.br_inep_censo_educacao_superior.ies`
        WHERE ano BETWEEN 2009 AND 2019
          AND id_municipio IS NOT NULL
      "),
      error = function(e) {
        message(sprintf("   IES table error: %s -- falling back to curso counts", e$message))
        # Fallback: use curso table, count distinct id_ies per municipality-year
        tryCatch(
          basedosdados::read_sql("
            SELECT
              ano,
              id_municipio,
              id_ies,
              tipo_organizacao_administrativa AS tipo_categoria_administrativa,
              tipo_organizacao_academica
            FROM `basedosdados.br_inep_censo_educacao_superior.curso`
            WHERE ano BETWEEN 2009 AND 2019
              AND id_municipio IS NOT NULL
            GROUP BY ano, id_municipio, id_ies,
                     tipo_organizacao_administrativa, tipo_organizacao_academica
          "),
          error = function(e2) {
            message(sprintf("   Fallback also failed: %s", e2$message))
            NULL
          }
        )
      }
    )
  },
  "inep_hei_counts"
)

if (!is.null(df_hei_raw)) {
  message(sprintf("   HEI raw: %d rows", nrow(df_hei_raw)))
} else {
  message("   SKIPPED: HEI data not available.")
}
gc()

# ==============================================================================
# SECTION 3K2: INEP HEI QUALITY SCORES -- IGC (institutional) + CPC (course)
# ------------------------------------------------------------------------------
# IGC (Indice Geral de Cursos) and CPC (Conceito Preliminar de Curso) are the
# two official MEC quality scales (values 1-5).  They are published at the HEI
# level by INEP.  We query the two canonical basedosdados tables and fall
# back gracefully if either is unavailable.
# ==============================================================================
message(" SECTION 3K2: Downloading INEP HEI quality scores (IGC / CPC)...")

df_hei_quality <- load_or_fetch(
  "Documentation_Final/cache_hei_quality.rds",
  function() {
    out <- list()
    # IGC table
    out$igc <- tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_ies,
          indice_geral_cursos          AS igc,
          indice_geral_cursos_continuo AS igc_continuo
        FROM `basedosdados.br_inep_indicadores_educacao_superior.ies`
        WHERE ano BETWEEN 2007 AND 2019
      "),
      error = function(e) {
        message(sprintf("   IGC primary table failed: %s", e$message))
        tryCatch(
          basedosdados::read_sql("
            SELECT
              ano,
              id_ies,
              AVG(CAST(conceito_igc AS FLOAT64))          AS igc,
              AVG(CAST(conceito_igc_continuo AS FLOAT64)) AS igc_continuo
            FROM `basedosdados.br_inep_censo_educacao_superior.ies`
            WHERE ano BETWEEN 2007 AND 2019
            GROUP BY ano, id_ies
          "),
          error = function(e2) {
            message(sprintf("   IGC fallback also failed: %s", e2$message))
            NULL
          })
      })
    # CPC table
    out$cpc <- tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_ies,
          AVG(CAST(conceito_preliminar_curso AS FLOAT64))          AS cpc,
          AVG(CAST(conceito_preliminar_curso_continuo AS FLOAT64)) AS cpc_continuo
        FROM `basedosdados.br_inep_indicadores_educacao_superior.curso`
        WHERE ano BETWEEN 2007 AND 2019
        GROUP BY ano, id_ies
      "),
      error = function(e) {
        message(sprintf("   CPC primary table failed: %s", e$message))
        NULL
      })
    out
  },
  "INEP HEI quality (IGC / CPC)"
)

# HEI-year wide quality panel
hei_quality_panel <- NULL
if (!is.null(df_hei_quality)) {
  parts <- list()
  if (!is.null(df_hei_quality$igc)) parts$igc <- df_hei_quality$igc %>% mutate(ano = as.integer(ano))
  if (!is.null(df_hei_quality$cpc)) parts$cpc <- df_hei_quality$cpc %>% mutate(ano = as.integer(ano))
  if (length(parts) == 2) {
    hei_quality_panel <- parts$igc %>% full_join(parts$cpc, by = c("ano", "id_ies"))
  } else if (length(parts) == 1) {
    hei_quality_panel <- parts[[1]]
  }
  if (!is.null(hei_quality_panel)) {
    # Pre-treatment (earliest-available) snapshot per HEI
    hei_quality_pre <- hei_quality_panel %>%
      group_by(id_ies) %>%
      summarise(
        igc_pre = first(na.omit(igc)),
        cpc_pre = first(na.omit(cpc)),
        .groups = "drop"
      )
    message(sprintf("   HEI quality panel: %d HEI-year obs | %d HEIs with pre-treatment IGC/CPC",
                    nrow(hei_quality_panel), nrow(hei_quality_pre)))
  }
} else {
  hei_quality_pre <- NULL
  message("   HEI quality panel unavailable; quality-drop robustness will be skipped.")
}

# ==============================================================================
# SECTION 3L: SIM HOMICIDE RATES (optional, best-effort)
# ==============================================================================
message(" SECTION 3L: Downloading SIM homicide rates (best-effort)...")

df_homicide_raw <- load_or_fetch(
  "Documentation_Final/cache_sim_homicide.rds",
  function() {
    # Try aggregated SIM table first (ICD-10 X85-Y09: assault homicides)
    tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_municipio,
          SUM(CAST(numero_obitos AS INT64)) AS homicidios
        FROM `basedosdados.br_ms_sim.municipio_causa`
        WHERE ano BETWEEN 2005 AND 2019
          AND (
            (causa_basica >= 'X85' AND causa_basica <= 'X99')
            OR (causa_basica >= 'Y00' AND causa_basica <= 'Y09')
          )
        GROUP BY ano, id_municipio
      "),
      error = function(e1) {
        message(sprintf("   SIM aggregated query failed: %s", e1$message))
        # Try alternative table name
        tryCatch(
          basedosdados::read_sql("
            SELECT
              ano,
              id_municipio,
              COUNT(*) AS homicidios
            FROM `basedosdados.br_ms_sim.microdados`
            WHERE ano BETWEEN 2005 AND 2019
              AND causabas >= 'X85' AND causabas <= 'Y09'
            GROUP BY ano, id_municipio
          "),
          error = function(e2) {
            message(sprintf("   SIM microdados also failed: %s", e2$message))
            message("   Homicide data unavailable -- will skip in models.")
            NULL
          }
        )
      }
    )
  },
  "sim_homicide_2005_2019"
)

if (!is.null(df_homicide_raw)) {
  message(sprintf("   Homicide raw: %d rows", nrow(df_homicide_raw)))
} else {
  message("   SKIPPED: SIM homicide data not available.")
}
gc()

# ==============================================================================
# SECTION 3M: POLICY CONFOUNDERS -- FIES, LEI 12.711 / 2012, REUNI
# ------------------------------------------------------------------------------
# Three federal policies co-expanded with ProUni and could confound the
# reduced-form wage estimates:
#
#   FIES (Fundo de Financiamento Estudantil): public student loans, strongly
#     expanded 2010-2014.  Queried from br_mec_fies if available; otherwise
#     we fall back to the INEP Censo ES flag on FIES-enrolled students.
#
#   Lei 12.711 / 2012 (federal quota law): affirmative action at federal HEIs
#     binding from 2013 for all federal institutions; we build an indicator
#     post_quota_t x federal_HEI_density_r from the INEP Censo ES.
#
#   REUNI (Decreto 6.096 / 2007): federal HEI expansion.  We use the count
#     of federal HEIs per microregion in 2009 (earliest Censo ES year) as a
#     REUNI exposure proxy.
# ==============================================================================
message(" SECTION 3M: Downloading FIES and building policy-confounder panels...")

df_fies_raw <- load_or_fetch(
  "Documentation_Final/cache_fies.rds",
  function() {
    tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_municipio,
          COUNT(*) AS n_fies
        FROM `basedosdados.br_mec_fies.microdados`
        WHERE ano BETWEEN 2005 AND 2019
        GROUP BY ano, id_municipio
      "),
      error = function(e) {
        message(sprintf("   FIES table not in basedosdados (%s) -- trying INEP fallback",
                        e$message))
        tryCatch(
          basedosdados::read_sql("
            SELECT
              ano,
              id_municipio,
              SUM(CAST(in_fies AS INT64)) AS n_fies
            FROM `basedosdados.br_inep_censo_educacao_superior.aluno`
            WHERE ano BETWEEN 2009 AND 2019
              AND in_fies = '1'
            GROUP BY ano, id_municipio
          "),
          error = function(e2) {
            message(sprintf("   FIES fallback failed: %s", e2$message))
            NULL
          })
      })
  },
  "FIES scholarships"
)

if (!is.null(df_fies_raw)) {
  message(sprintf("   FIES raw: %d rows", nrow(df_fies_raw)))
} else {
  message("   SKIPPED: FIES data not available; confounder control will be NA.")
}
gc()


# ==============================================================================
# SECTION 6C: MICROREGION-LEVEL COVARIATES FROM NEW SOURCES
# ==============================================================================
message(" SECTION 6C: Building time-varying covariates (PIB, HEI, homicide)...")

# --- 6C.1: PIB microregion-year panel ---
df_pib_micro <- NULL
tryCatch({
  df_pib_micro <- df_pib_raw %>%
    mutate(ano = as.integer(ano)) %>%
    inner_join(
      df_geo %>% distinct(id_municipio, id_microrregiao),
      by = "id_municipio"
    ) %>%
    # Join annual population for per-capita GDP
    left_join(
      df_pop_annual %>%
        mutate(ano = as.integer(ano)) %>%
        inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                   by = "id_municipio") %>%
        group_by(id_microrregiao, ano) %>%
        summarise(pop_total = sum(as.numeric(populacao), na.rm = TRUE),
                  .groups = "drop"),
      by = c("id_microrregiao", "ano")
    ) %>%
    group_by(id_microrregiao, ano) %>%
    summarise(
      pib_total      = sum(as.numeric(pib_nominal),        na.rm = TRUE),
      pop_total      = sum(as.numeric(pop_total),          na.rm = TRUE),
      va_total       = sum(as.numeric(valor_adicionado),   na.rm = TRUE),
      va_industry    = sum(as.numeric(va_industria),       na.rm = TRUE),
      va_services    = sum(as.numeric(va_servicos),        na.rm = TRUE),
      va_agro        = sum(as.numeric(va_agropecuaria),    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      id_num         = as.integer(id_microrregiao),
      log_pib        = log(pmax(pib_total, 1)),
      pib_pc         = ifelse(pop_total > 0, pib_total / pop_total, NA_real_),
      log_pib_pc     = log(pmax(pib_pc, 0.001)),
      share_industry = ifelse(va_total > 0, va_industry / va_total, NA_real_),
      share_services = ifelse(va_total > 0, va_services / va_total, NA_real_),
      share_agro     = ifelse(va_total > 0, va_agro     / va_total, NA_real_)
    ) %>%
    group_by(id_microrregiao) %>%
    arrange(ano) %>%
    mutate(
      pib_growth = (pib_total - lag(pib_total)) / pmax(lag(pib_total), 1)
    ) %>%
    ungroup()

  message(sprintf("   df_pib_micro: %d rows, %d microregions, years %d-%d",
                  nrow(df_pib_micro),
                  n_distinct(df_pib_micro$id_microrregiao),
                  as.integer(min(df_pib_micro$ano)), as.integer(max(df_pib_micro$ano))))
}, error = function(e) {
  message(sprintf("   WARNING: PIB microregion construction failed: %s", safe_err_msg(e)))
})

# --- 6C.2: Pre-treatment HEI snapshot (2009 -- earliest available) ---
df_hei_pre <- NULL
if (!is.null(df_hei_raw)) {
  tryCatch({
    df_hei_micro_full <- df_hei_raw %>%
      mutate(ano = as.integer(ano)) %>%
      inner_join(
        df_geo %>% distinct(id_municipio, id_microrregiao),
        by = "id_municipio"
      ) %>%
      group_by(id_microrregiao, ano) %>%
      summarise(
        n_heis_total  = n_distinct(id_ies),
        n_heis_private = sum(tipo_categoria_administrativa %in%
                               c("4", "5", "6", "7",
                                 "Privada com fins lucrativos",
                                 "Privada sem fins lucrativos",
                                 "Privada - com fins lucrativos",
                                 "Privada - sem fins lucrativos"),
                             na.rm = TRUE),
        n_universities = sum(tipo_organizacao_academica %in%
                               c("1", "Universidade"),
                             na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        log_n_priv_heis    = log(1 + n_heis_private),
        share_universities = n_universities / pmax(n_heis_total, 1)
      )

    # Use 2009 as pre-treatment proxy (earliest available)
    df_hei_pre <- df_hei_micro_full %>%
      filter(ano == 2009) %>%
      mutate(id_num = as.integer(id_microrregiao)) %>%
      select(id_num, id_microrregiao,
             n_private_heis_2009  = n_heis_private,
             log_n_priv_heis_2009 = log_n_priv_heis,
             n_heis_total_2009    = n_heis_total,
             share_univ_2009      = share_universities)

    message(sprintf("   df_hei_pre: %d microregions with 2009 HEI data",
                    nrow(df_hei_pre)))
  }, error = function(e) {
    message(sprintf("   WARNING: HEI aggregation failed: %s", safe_err_msg(e)))
  })
}

# --- 6C.3: Homicide rate panel ---
df_homicide_micro <- NULL
if (!is.null(df_homicide_raw)) {
  tryCatch({
    df_homicide_micro <- df_homicide_raw %>%
      mutate(ano = as.integer(ano)) %>%
      inner_join(
        df_geo %>% distinct(id_municipio, id_microrregiao),
        by = "id_municipio"
      ) %>%
      group_by(id_microrregiao, ano) %>%
      summarise(
        total_homicides = sum(as.numeric(homicidios), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # Join annual population
      left_join(
        df_pop_annual %>%
          mutate(ano = as.integer(ano)) %>%
          inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                     by = "id_municipio") %>%
          group_by(id_microrregiao, ano) %>%
          summarise(pop_total = sum(as.numeric(populacao), na.rm = TRUE),
                    .groups = "drop"),
        by = c("id_microrregiao", "ano")
      ) %>%
      mutate(
        id_num            = as.integer(id_microrregiao),
        homicide_rate     = (total_homicides / pmax(pop_total, 1)) * 100000,
        log_homicide_rate = log(1 + homicide_rate)
      )

    message(sprintf("   df_homicide_micro: %d microregion-years", nrow(df_homicide_micro)))
  }, error = function(e) {
    message(sprintf("   WARNING: Homicide microregion construction failed: %s", safe_err_msg(e)))
  })
}



# ==============================================================================
# SECTION 6D: UPDATE DR COVARIATE FRAME WITH PRE-TREATMENT HEI COUNTS
# ==============================================================================
# Pre-2005 HEI infrastructure directly addresses treatment endogeneity:
# microregions with more private universities were more likely to adopt ProUni
# earlier AND more intensively. Adding this to the DR propensity score controls
# for the supply-side selection mechanism.
#
# NOTE: We use 2009 as proxy for pre-2005 because basedosdados only has 2009+.
# Limitation noted in comments. For Low-Dose we can afford an extra covariate;
# for High-Dose we still cap at 3 to avoid singular matrix.
# ==============================================================================
message(" SECTION 6D: Updating DR covariate frames with HEI infrastructure...")

# Enrich df_covars_dr_enrich with HEI pre-treatment count (if available)
if (!is.null(df_hei_pre)) {
  df_covars_dr_enrich <- df_covars_dr_enrich %>%
    left_join(df_hei_pre %>% select(id_num, log_n_priv_heis_2009,
                                    n_private_heis_2009),
              by = "id_num") %>%
    mutate(
      log_n_priv_heis_2009 = replace_na(log_n_priv_heis_2009, 0),
      n_private_heis_2009  = replace_na(n_private_heis_2009, 0L)
    )

  # HEI-enriched DR formulas
  # Low-Dose: now 7 covariates (4 base + 2 census + 1 HEI)
  dr_xformla_low_hei  <- ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate +
                            share_in_school + college_density +
                            private_enroll_share + log_n_priv_heis_2009

  # High-Dose: 3 covariates (keep conservative -- HEI replaces private_enroll_share)
  dr_xformla_high_hei <- ~ log_pop_18_24 + avg_hh_income_pc_mw + log_n_priv_heis_2009

  message(sprintf("   df_covars_dr_enrich updated: %d rows, %d cols",
                  nrow(df_covars_dr_enrich), ncol(df_covars_dr_enrich)))
  message("   HEI DR formulas defined (Low: 7 covars | High: 3 covars)")

  # Re-run DR models with HEI covariate
  res_low_dr_hei  <- tryCatch({
    message("   -> Low Dose / DR + Censo + HEI...")
    run_cs(data_low, "log_wage_sm", label = "Low/DR/HEI",
           xformla   = dr_xformla_low_hei,
           covars_df = df_covars_dr_enrich)
  }, error = function(e) {
    message(sprintf("   WARNING: Low/DR/HEI failed: %s", safe_err_msg(e)))
    NULL
  })

  res_high_dr_hei <- tryCatch({
    message("   -> High Dose / DR + HEI...")
    run_cs(data_high, "log_wage_sm", label = "High/DR/HEI",
           xformla   = dr_xformla_high_hei,
           covars_df = df_covars_dr_enrich)
  }, error = function(e) {
    message(sprintf("   WARNING: High/DR/HEI failed: %s", safe_err_msg(e)))
    NULL
  })

  if (!is.null(res_low_dr_hei)) {
    message(sprintf("   Low/DR/HEI: ATT = %.4f (SE = %.4f, t = %.2f)",
                    res_low_dr_hei$simple$overall.att,
                    res_low_dr_hei$simple$overall.se,
                    res_low_dr_hei$simple$overall.att / res_low_dr_hei$simple$overall.se))
  }

  # Append to tab_att if it exists
  if (exists("tab_att") && !is.null(res_low_dr_hei)) {
    tab_att <- bind_rows(tab_att,
                         build_att_row(res_low_dr_hei),
                         if (!is.null(res_high_dr_hei)) build_att_row(res_high_dr_hei) else NULL)
    write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
    message("   ATT table updated with HEI DR models.")
  }
} else {
  message("   SKIPPED: HEI data not available, keeping existing DR formulas.")
  res_low_dr_hei  <- NULL
  res_high_dr_hei <- NULL
}

# ==============================================================================
# SECTION 6E: POLICY CONFOUNDER COVARIATES -- FIES + REUNI + LEI 12.711
# ------------------------------------------------------------------------------
# Builds three microregion-level confounder measures that enter either as
# pre-treatment propensity-score controls (Sant'Anna-Zhao DR) or as
# time-varying covariates in the TWFE robustness check:
#   - fies_density_pre : FIES grants per 1000 youth in the first year observed
#                        (proxy for pre-ProUni FIES exposure intensity);
#   - reuni_exposure   : log(1 + n_federal_HEIs_2009), a REUNI absorption
#                        proxy using the earliest Censo ES year;
#   - federal_share_2012 : share of federal HEIs in 2012; Lei 12.711 (2012)
#                        binds federal HEIs only from 2013, so a non-zero
#                        share identifies microregions exposed to the quota
#                        law.
# ==============================================================================
message(" SECTION 6E: Building FIES / REUNI / Lei 12.711 confounder covariates...")

df_fies_micro <- NULL
df_fies_panel <- NULL
tryCatch({
  if (!is.null(df_fies_raw)) {
    df_fies_micro <- df_fies_raw %>%
      mutate(ano = as.integer(ano)) %>%
      inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                 by = "id_municipio") %>%
      group_by(id_microrregiao, ano) %>%
      summarise(n_fies = sum(as.numeric(n_fies), na.rm = TRUE), .groups = "drop")

    fies_first_year <- min(df_fies_micro$ano, na.rm = TRUE)
    fies_snapshot <- df_fies_micro %>%
      filter(ano == fies_first_year) %>%
      inner_join(df_pop_micro, by = "id_microrregiao") %>%
      transmute(id_microrregiao,
                fies_density_pre = 1000 * n_fies / pmax(pop, 1))

    if (exists("df_covars_dr_enrich")) {
      df_covars_dr_enrich <- df_covars_dr_enrich %>%
        left_join(
          fies_snapshot %>%
            mutate(id_num = match(as.character(id_microrregiao),
                                   as.character(sort(unique(df_final$id_microrregiao))))) %>%
            select(id_num, fies_density_pre),
          by = "id_num"
        ) %>%
        mutate(fies_density_pre = replace_na(fies_density_pre, 0))
      message(sprintf("   df_covars_dr_enrich gained fies_density_pre (first FIES year = %d)",
                      as.integer(fies_first_year)))
    }

    # Time-varying FIES density panel for TWFE robustness
    df_fies_panel <- expand_grid(
      id_microrregiao = unique(df_pop_micro$id_microrregiao),
      ano             = params$years
    ) %>%
      left_join(df_fies_micro, by = c("id_microrregiao", "ano")) %>%
      left_join(df_pop_micro,  by = "id_microrregiao") %>%
      mutate(n_fies = replace_na(n_fies, 0),
             fies_density_rt = 1000 * n_fies / pmax(pop, 1))
  }
}, error = function(e) message(sprintf("   FIES covariate build failed: %s", safe_err_msg(e))))

# REUNI exposure via federal HEI count (if available in df_hei_pre)
if (exists("df_hei_pre") && !is.null(df_hei_pre) &&
    exists("df_covars_dr_enrich")) {
  if ("n_federal_heis_2009" %in% names(df_hei_pre)) {
    df_covars_dr_enrich <- df_covars_dr_enrich %>%
      left_join(df_hei_pre %>% select(id_num, n_federal_heis_2009),
                by = "id_num") %>%
      mutate(reuni_exposure = log1p(replace_na(n_federal_heis_2009, 0L)))
    message("   reuni_exposure = log1p(n_federal_heis_2009) added.")
  } else {
    message("   REUNI proxy skipped (no federal-HEI count field in df_hei_pre).")
  }
}

# Lei 12.711: pre-treatment federal-HEI share at microregion level (2012)
tryCatch({
  if (!is.null(df_hei_raw)) {
    is_federal <- function(x) grepl("federal", tolower(as.character(x)))
    fed_share_2012 <- df_hei_raw %>%
      mutate(ano = as.integer(ano)) %>%
      filter(ano == 2012) %>%
      inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                 by = "id_municipio") %>%
      group_by(id_microrregiao) %>%
      summarise(
        n_ies_total   = n_distinct(id_ies),
        n_ies_federal = n_distinct(id_ies[is_federal(tipo_categoria_administrativa)]),
        .groups = "drop"
      ) %>%
      mutate(federal_share_2012 = n_ies_federal / pmax(n_ies_total, 1))

    if (exists("df_covars_dr_enrich")) {
      df_covars_dr_enrich <- df_covars_dr_enrich %>%
        left_join(
          fed_share_2012 %>%
            mutate(id_num = match(as.character(id_microrregiao),
                                   as.character(sort(unique(df_final$id_microrregiao))))) %>%
            select(id_num, federal_share_2012),
          by = "id_num"
        ) %>%
        mutate(federal_share_2012 = replace_na(federal_share_2012, 0))
      message("   federal_share_2012 added (Lei 12.711 exposure proxy).")
    }
  }
}, error = function(e) message(sprintf("   Lei 12.711 covariate build failed: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9Z: POLICY-CONFOUNDER-CONTROLLED DR + TWFE RE-RUNS
# ==============================================================================
message(" SECTION 9Z: Policy-confounder-controlled DR and TWFE re-runs...")

res_low_dr_policy  <- NULL
res_high_dr_policy <- NULL
twfe_fies          <- NULL

tryCatch({
  if (exists("df_covars_dr_enrich")) {
    extra <- c()
    if ("fies_density_pre"   %in% names(df_covars_dr_enrich)) extra <- c(extra, "fies_density_pre")
    if ("reuni_exposure"     %in% names(df_covars_dr_enrich)) extra <- c(extra, "reuni_exposure")
    if ("federal_share_2012" %in% names(df_covars_dr_enrich)) extra <- c(extra, "federal_share_2012")

    if (length(extra) > 0) {
      base_low_terms  <- c("pop_18_24", "avg_hh_income_pc_mw", "literacy_rate",
                           "share_in_school", "college_density", "private_enroll_share")
      base_high_terms <- c("log_pop_18_24", "avg_hh_income_pc_mw", "college_density")
      f_low  <- reformulate(c(base_low_terms,  extra))
      f_high <- reformulate(c(base_high_terms, extra))
      message(sprintf("   Policy-augmented xformla (Low):  %s",
                      paste(deparse(f_low),  collapse = " ")))
      message(sprintf("   Policy-augmented xformla (High): %s",
                      paste(deparse(f_high), collapse = " ")))

      res_low_dr_policy  <- tryCatch(
        run_cs(data_low, "log_wage_sm", label = "Low/DR/Policy",
               xformla = f_low,  covars_df = df_covars_dr_enrich,
               est_method = "dr"),
        error = function(e) { message(sprintf("   Low/DR/Policy failed: %s",
                                              safe_err_msg(e))); NULL })
      res_high_dr_policy <- tryCatch(
        run_cs(data_high, "log_wage_sm", label = "High/DR/Policy",
               xformla = f_high, covars_df = df_covars_dr_enrich,
               est_method = "dr"),
        error = function(e) { message(sprintf("   High/DR/Policy failed: %s",
                                              safe_err_msg(e))); NULL })

      for (r in Filter(Negate(is.null), list(res_low_dr_policy, res_high_dr_policy))) {
        tab_att <- bind_rows(tab_att, build_att_row(r))
      }
      write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
    } else {
      message("   SKIP: no policy-confounder columns were built.")
    }
  }
}, error = function(e) message(sprintf("   Policy DR re-run failed: %s", safe_err_msg(e))))

# Time-varying TWFE with fies_density_rt (direct robustness check)
tryCatch({
  if (!is.null(df_fies_panel)) {
    twfe_fies_data <- df_final %>%
      mutate(treated = as.integer(g > 0 & ano >= g)) %>%
      left_join(df_fies_panel %>% select(id_microrregiao, ano, fies_density_rt),
                by = c("id_microrregiao", "ano")) %>%
      mutate(fies_density_rt = replace_na(fies_density_rt, 0))
    twfe_fies <- fixest::feols(log_wage_sm ~ treated + fies_density_rt |
                                 id_microrregiao + ano,
                               data    = twfe_fies_data,
                               cluster = ~id_microrregiao)
    message("   TWFE with fies_density_rt covariate:")
    print(fixest::coeftable(twfe_fies))

    tc <- tryCatch(fixest::coeftable(twfe_fies), error = function(e) NULL)
    if (!is.null(tc) && "treated" %in% rownames(tc)) {
      tab_att <- bind_rows(tab_att, tibble(
        Model      = "TWFE/Wages+FIES",
        ATT        = round(tc["treated", "Estimate"],   4),
        SE         = round(tc["treated", "Std. Error"], 4),
        CI_low     = round(tc["treated", "Estimate"] - 1.96 * tc["treated", "Std. Error"], 4),
        CI_high    = round(tc["treated", "Estimate"] + 1.96 * tc["treated", "Std. Error"], 4),
        t_val      = round(tc["treated", "t value"],    3),
        p_val      = round(tc["treated", "Pr(>|t|)"],   4),
        Pretrend_p = NA_real_,
        N_obs      = as.integer(stats::nobs(twfe_fies)),
        N_treated  = NA_integer_,
        N_notyet   = NA_integer_
      ))
      write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
    }
  }
}, error = function(e) message(sprintf("   TWFE+FIES failed: %s", safe_err_msg(e))))

# Balance table -- pre-treatment GDP and HEI by dose bin
tryCatch({
  if (!is.null(df_pib_micro) && !is.null(df_hei_pre)) {
    df_balance_enriched <- df_final %>%
      filter(ano == max(params$years), !is.na(dose_bin)) %>%
      distinct(id_microrregiao, dose_bin) %>%
      left_join(df_hei_pre %>% select(id_microrregiao,
                                       n_private_heis_2009,
                                       log_n_priv_heis_2009),
                by = "id_microrregiao") %>%
      left_join(
        df_pib_micro %>%
          filter(ano == 2008) %>%
          select(id_microrregiao, pib_pc_2008 = pib_pc,
                 share_industry_2008 = share_industry),
        by = "id_microrregiao"
      )

    balance_hei_summary <- df_balance_enriched %>%
      filter(!is.na(dose_bin)) %>%
      group_by(dose_bin) %>%
      summarise(
        mean_private_heis = round(mean(n_private_heis_2009, na.rm = TRUE), 1),
        sd_private_heis   = round(sd(n_private_heis_2009,   na.rm = TRUE), 1),
        mean_pib_pc_2008  = round(mean(pib_pc_2008,  na.rm = TRUE), 0),
        sd_pib_pc_2008    = round(sd(pib_pc_2008,    na.rm = TRUE), 0),
        mean_ind_share    = round(mean(share_industry_2008, na.rm = TRUE), 3),
        .groups = "drop"
      )
    write_csv(balance_hei_summary, "Tables_Final/tab_balance_hei_v5.csv")
    message("   Balance table (HEI + GDP) saved.")
  }
}, error = function(e) message(sprintf("   WARNING: Balance table enrichment failed: %s", safe_err_msg(e))))



# ==============================================================================
# SECTION 9W: BIRTH-COHORT QUADRATIC BATTERY WITH PROGRESSIVE GDP CONTROLS
# ==============================================================================
# Purpose: Verify that the concave dose-response in the birth-cohort design
# survives inclusion of time-varying macroeconomic controls. If the quadratic
# (D_rt:post + D_rt^2:post) remains stable across columns (1)-(4), the finding
# is not driven by differential economic shocks across microregions.
# Reference: Addendum B, model battery specification.
# ==============================================================================
message(" SECTION 9W: Birth-cohort quadratic battery with GDP controls...")

tryCatch({

  # 1. Defensively re-create cont_bc_data if Section 9Q result not in scope
  if (!exists("cont_bc_data") || is.null(cont_bc_data)) {
    if (!is.null(df_bc_wide) && nrow(df_bc_wide) > 0) {
      cont_bc_data <- df_bc_wide %>%
        mutate(post = as.integer(g > 0 & ano >= g))
      message("   cont_bc_data re-created from df_bc_wide")
    } else {
      stop("df_bc_wide not available -- Section 9W cannot proceed")
    }
  }

  bc_data_w <- cont_bc_data %>%
    filter(!is.na(wage_gap_bc), !is.na(D_rt), !is.na(post))

  message(sprintf("   Base birth-cohort sample: %d microregion-year obs", nrow(bc_data_w)))

  # 2. Merge time-varying PIB controls (microregion x year)
  if (!is.null(df_pib_micro)) {
    bc_data_w <- bc_data_w %>%
      left_join(
        df_pib_micro %>%
          mutate(id_num = as.integer(id_microrregiao)) %>%
          select(id_num, ano, log_pib, pib_growth, share_industry),
        by = c("id_num", "ano")
      )
    message(sprintf("   PIB joined: %d/%d obs with log_pib",
                    sum(!is.na(bc_data_w$log_pib)), nrow(bc_data_w)))
  } else {
    bc_data_w <- bc_data_w %>%
      mutate(log_pib = NA_real_, pib_growth = NA_real_, share_industry = NA_real_)
    message("   WARNING: df_pib_micro NULL -- GDP controls unavailable")
  }

  # 3. Enrollment: re-aggregate Education Census 2009 to microregion
  #    (pre-treatment snapshot applied as constant cross-section; logs absorb level diffs)
  if (exists("df_censo_educ") && !is.null(df_censo_educ) &&
      exists("df_geo")        && !is.null(df_geo)) {
    log_enroll_2009 <- tryCatch({
      df_censo_educ %>%
        inner_join(
          df_geo %>% distinct(id_municipio, id_microrregiao),
          by = "id_municipio"
        ) %>%
        group_by(id_microrregiao) %>%
        summarise(
          total_ingressantes_micro = sum(total_ingressantes, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        transmute(
          id_num         = as.integer(id_microrregiao),
          log_enrollment = log(1 + total_ingressantes_micro)
        )
    }, error = function(e) {
      message(sprintf("   WARNING: enrollment aggregation failed: %s", e$message))
      NULL
    })
    if (!is.null(log_enroll_2009)) {
      bc_data_w <- bc_data_w %>%
        left_join(log_enroll_2009, by = "id_num") %>%
        mutate(log_enrollment = replace_na(log_enrollment, 0))
      message(sprintf("   Enrollment (2009) joined: %d obs non-zero",
                      sum(bc_data_w$log_enrollment > 0, na.rm = TRUE)))
    } else {
      bc_data_w <- bc_data_w %>% mutate(log_enrollment = NA_real_)
    }
  } else {
    bc_data_w <- bc_data_w %>% mutate(log_enrollment = NA_real_)
    message("   Enrollment control unavailable (df_censo_educ not in scope)")
  }

  # 4. Merge time-varying homicide rate (microregion x year)
  if (!is.null(df_homicide_micro)) {
    bc_data_w <- bc_data_w %>%
      left_join(
        df_homicide_micro %>%
          select(id_num, ano, log_homicide_rate),
        by = c("id_num", "ano")
      )
    message(sprintf("   Homicide joined: %d/%d obs non-NA",
                    sum(!is.na(bc_data_w$log_homicide_rate)), nrow(bc_data_w)))
  } else {
    bc_data_w <- bc_data_w %>% mutate(log_homicide_rate = NA_real_)
    message("   WARNING: df_homicide_micro NULL -- homicide control unavailable")
  }

  # Helper: run feols safely and report key coefficients
  safe_feols_bc <- function(fml, data, label) {
    tryCatch({
      if (nrow(data) < 50) {
        message(sprintf("   %-28s | SKIPPED (n=%d < 50)", label, nrow(data)))
        return(NULL)
      }
      fit <- feols(fml, data = data, cluster = ~id_num)
      ct  <- coeftable(fit)
      if (nrow(ct) >= 2) {
        message(sprintf("   %-28s | D=%+.5f (t=%+.2f) | D^2=%+.7f (t=%+.2f)",
                        label, ct[1, 1], ct[1, 3], ct[2, 1], ct[2, 3]))
      }
      fit
    }, error = function(e) {
      message(sprintf("   %-28s | FAILED: %s", label, e$message))
      NULL
    })
  }

  # Pre-filter subsets to avoid repeated inline filtering
  bc_full    <- bc_data_w
  bc_pib     <- bc_data_w %>% filter(!is.na(log_pib))
  bc_pib_gr  <- bc_data_w %>% filter(!is.na(log_pib), !is.na(pib_growth))
  bc_pib_ind <- bc_data_w %>% filter(!is.na(log_pib), !is.na(share_industry))
  bc_pib_enr <- bc_data_w %>% filter(!is.na(log_pib), !is.na(pib_growth),
                                      log_enrollment > 0)
  bc_pib_hom <- bc_data_w %>% filter(!is.na(log_pib), !is.na(log_homicide_rate))

  # Model 1: Baseline quadratic (replicates fit_bc_quad from Section 9Q)
  fit_bc_1 <- safe_feols_bc(
    wage_gap_bc ~ D_rt:post + I(D_rt^2):post | id_num + ano,
    bc_full, "(1) Baseline"
  )

  # Model 2: + log(GDP level) x post -- absorbs differential economic levels
  fit_bc_2 <- safe_feols_bc(
    wage_gap_bc ~ D_rt:post + I(D_rt^2):post + log_pib:post | id_num + ano,
    bc_pib, "(2) + log GDP"
  )

  # Model 3: + GDP growth x post -- absorbs differential economic momentum
  fit_bc_3 <- safe_feols_bc(
    wage_gap_bc ~ D_rt:post + I(D_rt^2):post + log_pib:post + pib_growth:post | id_num + ano,
    bc_pib_gr, "(3) + GDP Growth"
  )

  # Model 4: + industry share x post -- structural composition of local economy
  fit_bc_4 <- safe_feols_bc(
    wage_gap_bc ~ D_rt:post + I(D_rt^2):post + log_pib:post + share_industry:post | id_num + ano,
    bc_pib_ind, "(4) + Industry Share"
  )

  # Model 5: Kitchen sink -- adds HE enrollment (only if obs sufficient)
  fit_bc_5 <- if (nrow(bc_pib_enr) >= 50 && !all(is.na(bc_pib_enr$log_enrollment))) {
    safe_feols_bc(
      wage_gap_bc ~ D_rt:post + I(D_rt^2):post + log_pib:post + pib_growth:post +
        log_enrollment:post | id_num + ano,
      bc_pib_enr, "(5) Kitchen Sink"
    )
  } else {
    message("   (5) Kitchen Sink SKIPPED (enrollment unavailable or n < 50)")
    NULL
  }

  # Model 6: + homicide rate x post -- institutional quality / security environment
  fit_bc_6 <- if (nrow(bc_pib_hom) >= 50 && !all(is.na(bc_pib_hom$log_homicide_rate))) {
    safe_feols_bc(
      wage_gap_bc ~ D_rt:post + I(D_rt^2):post + log_pib:post +
        log_homicide_rate:post | id_num + ano,
      bc_pib_hom, "(6) + Homicide"
    )
  } else {
    message("   (6) + Homicide SKIPPED (homicide data unavailable or n < 50)")
    NULL
  }

  # --- Delta method: peak dose D* = -beta1 / (2 * beta2) ---
  # Only meaningful if beta2 < 0 (concave response)
  compute_dstar_bc <- function(fit, label) {
    if (is.null(fit)) return(NA_real_)
    tryCatch({
      ct <- coeftable(fit)
      rn <- rownames(ct)
      # Identify the D_rt:post coefficient (linear) and I(D_rt^2):post (quadratic)
      idx_b1 <- which(rn == "D_rt:post")
      idx_b2 <- grep("D_rt.2..post", rn)   # regex: matches "I(D_rt^2):post"
      if (length(idx_b1) == 0 || length(idx_b2) == 0) return(NA_real_)
      b1 <- ct[idx_b1[1], "Estimate"]
      b2 <- ct[idx_b2[1], "Estimate"]
      if (is.na(b1) || is.na(b2) || b2 >= 0) return(NA_real_)  # concavity required
      d_star <- -b1 / (2 * b2)
      message(sprintf("   D* (%s) = %.2f scholarships per 1,000 youth", label, d_star))
      d_star
    }, error = function(e) NA_real_)
  }

  d_star_9w <- list(
    bc_1 = compute_dstar_bc(fit_bc_1, "(1)"),
    bc_2 = compute_dstar_bc(fit_bc_2, "(2)"),
    bc_3 = compute_dstar_bc(fit_bc_3, "(3)"),
    bc_4 = compute_dstar_bc(fit_bc_4, "(4)")
  )

  # --- modelsummary robustness table ---
  bc_model_list_9w <- Filter(Negate(is.null), list(
    "(1) Baseline"        = fit_bc_1,
    "(2) + Log GDP"       = fit_bc_2,
    "(3) + GDP Growth"    = fit_bc_3,
    "(4) + Indust. Share" = fit_bc_4,
    "(5) Kitchen Sink"    = fit_bc_5,
    "(6) + Homicide"      = fit_bc_6
  ))

  if (length(bc_model_list_9w) >= 1) {

    bc_coef_map_9w <- c(
      "D_rt:post"              = "Dose $\\times$ Post",
      "I(D_rt^2):post"         = "Dose$^2$ $\\times$ Post",
      "log_pib:post"           = "Log GDP $\\times$ Post",
      "pib_growth:post"        = "GDP Growth $\\times$ Post",
      "share_industry:post"    = "Industry Share $\\times$ Post",
      "log_enrollment:post"    = "Log Enrollment $\\times$ Post",
      "log_homicide_rate:post" = "Log Homicide $\\times$ Post"
    )

    bc_gof_9w <- tribble(
      ~raw,            ~clean,           ~fmt,
      "nobs",          "Observations",   0,
      "r.squared",     "$R^2$",          3,
      "adj.r.squared", "Adj.\\ $R^2$",   3
    )

    bc_note_9w <- paste0(
      "Dependent variable: cohort wage gap (log wage, exposed cohort 1987--93 ",
      "minus control cohort 1975--82) aggregated to microregion $\\times$ year. ",
      "Dose = ProUni scholarships per 1,000 youth (18--24). ",
      "All controls interacted with Post indicator. ",
      "Microregion and year fixed effects included in all columns. ",
      "Standard errors clustered by microregion. ",
      "GDP: IBGE PIB municipal (2005--2019), aggregated to microregion. ",
      "Enrollment: INEP Education Census 2009 (pre-treatment snapshot). ",
      "Homicide: SIM/DATASUS ICD-10 X85--Y09."
    )

    # LaTeX export
    tryCatch({
      modelsummary(
        bc_model_list_9w,
        coef_map = bc_coef_map_9w,
        gof_map  = bc_gof_9w,
        stars    = c("*" = .1, "**" = .05, "***" = .01),
        title    = paste0("Cohort Wage Gap: Quadratic Dose--Response Robustness",
                          " to Macroeconomic Controls"),
        notes    = bc_note_9w,
        output   = "Tables_Final/tab_bc_quad_robustness_v5.tex",
        escape   = FALSE
      )
      message("   Saved: Tables_Final/tab_bc_quad_robustness_v5.tex")
    }, error = function(e)
      message(sprintf("   WARNING: LaTeX export failed: %s", e$message)))

    # CSV export (more robust fallback)
    tryCatch({
      modelsummary(
        bc_model_list_9w,
        coef_map = bc_coef_map_9w,
        gof_map  = bc_gof_9w,
        stars    = c("*" = .1, "**" = .05, "***" = .01),
        output   = "Tables_Final/tab_bc_quad_robustness_v5.csv"
      )
      message("   Saved: Tables_Final/tab_bc_quad_robustness_v5.csv")
    }, error = function(e)
      message(sprintf("   WARNING: CSV export failed: %s", e$message)))

    # Save D* robustness summary
    d_star_df_9w <- tibble(
      spec     = names(d_star_9w),
      d_star   = unlist(d_star_9w),
      controls = c("Baseline", "+ Log GDP", "+ GDP Growth", "+ Indust. Share")
    ) %>% filter(!is.na(d_star))

    if (nrow(d_star_df_9w) > 0) {
      write_csv(d_star_df_9w, "Tables_Final/tab_bc_dstar_v5.csv")
      message(sprintf("   D* range: [%.2f, %.2f] scholarships/1,000 youth",
                      min(d_star_df_9w$d_star), max(d_star_df_9w$d_star)))
    } else {
      message("   NOTE: No concave quadratic -- D* not defined (b2 >= 0 in all models)")
    }

    message(sprintf("   Section 9W complete: %d/%d models estimated",
                    length(bc_model_list_9w), 6L))

  } else {
    message("   WARNING: Zero models estimated in Section 9W -- table skipped")
  }

}, error = function(e) {
  message(sprintf("   SECTION 9W FAILED: %s", safe_err_msg(e)))
})


# ==============================================================================
# SECTION 9X: COMMITTEE DEFENSE BATTERY
# ==============================================================================
# Responses to examiner probes on the birth-cohort quadratic, inference, and
# policy magnitudes. Produces four artefacts:
#   (1) Effective-cluster count in the treatment arm (Q6)
#   (2) HonestDiD breakdown M-bar for Low/Wages (Q7)
#   (3) Fiscal magnitude: BRL per scholarship vs. programme cost (Q12)
#   (4) Attenuation diagnostic for fit_bc_1 vs fit_bc_2 (Q10)
# ==============================================================================
message(" SECTION 9X: Committee defense battery...")

tryCatch({

  # ---------------------------------------------------------------------------
  # (1) Fieller-method CI for D* = -beta1 / (2*beta2)  --  REMOVED IN FINAL DRAFT
  # ---------------------------------------------------------------------------
  # The Fieller subblock and the writing of tab_fieller_dstar_v5.csv have been
  # commented out: the peak-dose point estimate is no longer reported as a
  # primary object in the final manuscript.
  #
  # fieller_dstar <- function(fit, alpha = 0.05, label = "") { ... }
  # fieller_tab   <- bind_rows( fieller_dstar(fit_bc_1, ...), ... )
  # write_csv(fieller_tab, "Tables_Final/tab_fieller_dstar_v5.csv")

  # ---------------------------------------------------------------------------
  # (2) Effective-cluster count in treated arm (Low and High)
  # ---------------------------------------------------------------------------
  #     MacKinnon & Webb (2018): wild bootstrap inference requires a sufficient
  #     number of *treated* clusters. Report n_treated as well as the share of
  #     clusters with any post-period treatment variation.
  effective_clusters <- function(dat, label) {
    tryCatch({
      clus <- dat %>%
        group_by(id_num) %>%
        summarise(has_treat = any(g > 0, na.rm = TRUE),
                  n_obs     = n(), .groups = "drop")
      tibble(
        spec            = label,
        n_clusters      = nrow(clus),
        n_treated       = sum(clus$has_treat),
        n_notyet        = sum(!clus$has_treat),
        share_treated   = mean(clus$has_treat)
      )
    }, error = function(e)
      tibble(spec = label, n_clusters = NA, n_treated = NA,
             n_notyet = NA, share_treated = NA))
  }

  eff_tab <- bind_rows(
    if (exists("data_low"))  effective_clusters(data_low,  "Low/Wages")  else NULL,
    if (exists("data_high")) effective_clusters(data_high, "High/Wages") else NULL,
    if (exists("data_mid"))  effective_clusters(data_mid,  "Mid/Wages")  else NULL
  )
  if (!is.null(eff_tab) && nrow(eff_tab) > 0) {
    write_csv(eff_tab, "Tables_Final/tab_effective_clusters_v5.csv")
    message("   Effective clusters:")
    for (i in seq_len(nrow(eff_tab))) {
      r <- eff_tab[i, ]
      message(sprintf("      %-12s n_clusters=%d  n_treated=%d  n_notyet=%d  share=%.2f",
                      r$spec, r$n_clusters, r$n_treated, r$n_notyet, r$share_treated))
    }
  }

  # ---------------------------------------------------------------------------
  # (3) HonestDiD breakdown M-bar for Low/Wages
  # ---------------------------------------------------------------------------
  #     Breakdown M-bar = smallest M for which 0 enters the robust 95% CI.
  #     Interpreted as the degree of violation of parallel trends required to
  #     overturn a statistically significant result.
  breakdown_mbar <- NA_real_
  if (exists("honest_low_wages") && !is.null(honest_low_wages)) {
    tryCatch({
      hld <- honest_low_wages
      # Expect columns: Mbar, lb, ub
      if (all(c("Mbar", "lb", "ub") %in% names(hld))) {
        crosses <- hld %>%
          mutate(includes_zero = (lb <= 0 & ub >= 0)) %>%
          filter(includes_zero) %>%
          arrange(Mbar)
        if (nrow(crosses) > 0) {
          breakdown_mbar <- crosses$Mbar[1]
          message(sprintf("   HonestDiD breakdown M-bar (Low/Wages) = %.3f",
                          breakdown_mbar))
        } else {
          message("   HonestDiD (Low/Wages): CI excludes zero for all M-bar in grid")
          breakdown_mbar <- max(hld$Mbar, na.rm = TRUE) + 0.01   # report ">" lower bound
        }
        write_csv(tibble(spec = "Low/Wages", breakdown_mbar = breakdown_mbar),
                  "Tables_Final/tab_honestdid_breakdown_v5.csv")
      } else {
        message("   HonestDiD table missing expected columns (Mbar/lb/ub)")
      }
    }, error = function(e)
      message(sprintf("   HonestDiD breakdown calc failed: %s", e$message)))
  } else {
    message("   honest_low_wages not in scope -- breakdown M-bar skipped")
  }

  # ---------------------------------------------------------------------------
  # (4) Fiscal magnitude: BRL per scholarship vs programme cost
  # ---------------------------------------------------------------------------
  # ATT is on avg log wage of formal workers aged 22-35 in microregion-year.
  # Conversion: delta_wage_per_worker = exp(att) - 1; total_wage_gain_per_year
  # = delta_wage_per_worker * mean_monthly_wage * 12 * n_formal_workers.
  # Per-scholarship gain = total / n_scholarships_in_microregion.
  # Benchmark: ProUni tax expenditure approximately R$5,000-7,000/scholarship/yr
  # (Ministry of Finance estimates, CEPAL 2015).
  fiscal_tab <- tryCatch({
    if (exists("res_low_dr") && !is.null(res_low_dr) &&
        exists("df_final")   && !is.null(df_final)) {

      att_low <- res_low_dr$simple$overall.att
      se_low  <- res_low_dr$simple$overall.se

      # Minimum wage in BRL per month -- use 2019 real value as numeraire
      mw_brl <- 998     # 2019 minimum wage, BRL

      # From df_final: mean log wage (in min-wage units) among Low-Dose microregions
      if ("log_wage_sm" %in% names(df_final)) {
        mean_wage_mw <- df_final %>%
          filter(dose_bin == "Low", ano >= 2009) %>%
          summarise(m = mean(exp(log_wage_sm), na.rm = TRUE)) %>%
          pull(m)
      } else {
        mean_wage_mw <- 2.0   # conservative fallback (2x min wage)
      }
      mean_monthly_brl <- mean_wage_mw * mw_brl

      # Avg formal workers aged 22-35 per Low-Dose microregion-year
      if ("n_formal" %in% names(df_final)) {
        avg_n_formal <- df_final %>%
          filter(dose_bin == "Low", ano >= 2009) %>%
          summarise(m = mean(n_formal, na.rm = TRUE)) %>% pull(m)
      } else {
        avg_n_formal <- 5000   # order-of-magnitude placeholder
      }

      # Avg scholarships per Low-Dose microregion-year (post-treatment)
      avg_n_scholarships <- df_final %>%
        filter(dose_bin == "Low", ano >= 2009) %>%
        summarise(m = mean(D_rt * avg_n_formal / 1000, na.rm = TRUE)) %>% pull(m)
      if (!is.finite(avg_n_scholarships) || avg_n_scholarships < 1) avg_n_scholarships <- 50

      # Wage gain per worker per month and per year
      pct_gain           <- exp(att_low) - 1
      gain_per_worker_mo <- pct_gain * mean_monthly_brl
      gain_total_yr      <- gain_per_worker_mo * 12 * avg_n_formal
      gain_per_schol_yr  <- gain_total_yr / pmax(avg_n_scholarships, 1)

      # ProUni tax-expenditure benchmark (2019 BRL, per scholarship per year)
      prouni_cost_low  <- 5000
      prouni_cost_high <- 7000

      out <- tibble(
        spec                    = "Low/DR",
        att                     = att_low,
        se                      = se_low,
        pct_gain                = pct_gain,
        mean_monthly_wage_brl   = mean_monthly_brl,
        avg_n_formal_workers    = avg_n_formal,
        avg_n_scholarships_yr   = avg_n_scholarships,
        gain_per_worker_month   = gain_per_worker_mo,
        gain_total_yr_brl       = gain_total_yr,
        gain_per_scholarship_yr = gain_per_schol_yr,
        cost_per_scholarship_lo = prouni_cost_low,
        cost_per_scholarship_hi = prouni_cost_high,
        bcr_low                 = gain_per_schol_yr / prouni_cost_high,
        bcr_high                = gain_per_schol_yr / prouni_cost_low
      )
      write_csv(out, "Tables_Final/tab_fiscal_magnitude_v5.csv")
      message(sprintf("   Fiscal magnitude (Low/DR):"))
      message(sprintf("      ATT = %+.4f  =>  %+.2f%% of mean wage", att_low, 100 * pct_gain))
      message(sprintf("      Gain per worker per month    : R$ %.2f", gain_per_worker_mo))
      message(sprintf("      Gain per microregion per year: R$ %.0f (%.0f workers)",
                      gain_total_yr, avg_n_formal))
      message(sprintf("      Gain per scholarship per year: R$ %.0f", gain_per_schol_yr))
      message(sprintf("      ProUni cost per scholarship  : R$ %d-%d",
                      prouni_cost_low, prouni_cost_high))
      message(sprintf("      Implied benefit-cost ratio   : [%.2f, %.2f]",
                      out$bcr_low, out$bcr_high))
      out
    } else NULL
  }, error = function(e) {
    message(sprintf("   Fiscal magnitude calc failed: %s", e$message))
    NULL
  })

  # ---------------------------------------------------------------------------
  # (5) Attenuation diagnostic: fit_bc_1 vs fit_bc_2 coefficient stability
  # ---------------------------------------------------------------------------
  # If adding log_pib:post shrinks the D_rt:post coefficient substantially,
  # the quadratic is partly capturing a GDP story. If the coefficient is
  # stable, the concave response is orthogonal to local GDP.
  if (exists("fit_bc_1") && !is.null(fit_bc_1) &&
      exists("fit_bc_2") && !is.null(fit_bc_2)) {
    tryCatch({
      b1_m1 <- coef(fit_bc_1)["D_rt:post"]
      b1_m2 <- coef(fit_bc_2)["D_rt:post"]
      pct_change <- 100 * (b1_m2 - b1_m1) / abs(b1_m1)
      message(sprintf("   Attenuation from GDP control: %+.1f%% (%.5f -> %.5f)",
                      pct_change, b1_m1, b1_m2))
      write_csv(tibble(
        baseline_coef    = b1_m1,
        gdp_adjusted     = b1_m2,
        pct_change       = pct_change,
        orthogonal_to_gdp = abs(pct_change) < 20
      ), "Tables_Final/tab_attenuation_gdp_v5.csv")
    }, error = function(e)
      message(sprintf("   Attenuation diag failed: %s", e$message)))
  }

  message("   Section 9X complete.")

}, error = function(e) {
  message(sprintf("   SECTION 9X FAILED: %s", safe_err_msg(e)))
})

# ==============================================================================
# SECTION 9Y: HEI-LEVEL ROBUSTNESS -- IGC/CPC QUALITY FILTERS
# ------------------------------------------------------------------------------
# With id_ies now retained in the ProUni microdata, we can:
#   (A) compute an HEI-level first stage -- total ProUni scholarships per HEI
#       per year divided by HEI enrolment;
#   (B) tag each scholarship with the issuing HEI's pre-treatment IGC/CPC
#       quality tier, then re-run the Low-Dose wage model after dropping the
#       microregions whose exposure is dominated by bottom-quartile HEIs.
#
# If either id_ies or IGC/CPC is unavailable (legacy cache, table missing),
# the section degrades gracefully to a NULL result.
# ==============================================================================
message(" SECTION 9Y: HEI-level quality-drop robustness...")

res_low_hq  <- NULL
res_low_noq <- NULL
tab_hei_fs  <- NULL

tryCatch({
  has_ies <- "id_ies" %in% names(df_prouni_raw) &&
             any(!is.na(df_prouni_raw$id_ies))
  has_q   <- !is.null(hei_quality_pre)

  if (!has_ies) {
    message("   SKIP: ProUni cache lacks id_ies (rerun after clearing cache).")
  } else if (!has_q) {
    message("   SKIP: HEI quality panel unavailable (IGC/CPC missing).")
  } else {
    # --- (A) HEI-level first stage ---
    # scholarships per HEI-year and (best-effort) HEI enrolment from Censo ES
    prouni_by_ies_year <- df_prouni_raw %>%
      mutate(ano = as.integer(ano)) %>%
      filter(!is.na(id_ies)) %>%
      group_by(id_ies, ano) %>%
      summarise(n_bolsas = sum(n_bolsas, na.rm = TRUE), .groups = "drop")

    hei_enrol <- tryCatch(
      basedosdados::read_sql("
        SELECT
          ano,
          id_ies,
          SUM(CAST(quantidade_matriculas AS INT64)) AS n_matriculas
        FROM `basedosdados.br_inep_censo_educacao_superior.ies`
        WHERE ano BETWEEN 2007 AND 2019
        GROUP BY ano, id_ies
      "),
      error = function(e) NULL)

    if (!is.null(hei_enrol)) {
      hei_fs_panel <- prouni_by_ies_year %>%
        inner_join(hei_enrol %>% mutate(ano = as.integer(ano)),
                   by = c("id_ies", "ano")) %>%
        left_join(hei_quality_pre, by = "id_ies") %>%
        mutate(prouni_share = n_bolsas / pmax(n_matriculas, 1))

      tab_hei_fs <- hei_fs_panel %>%
        summarise(
          n_hei_years      = n(),
          n_ies            = n_distinct(id_ies),
          mean_prouni_share = mean(prouni_share, na.rm = TRUE),
          med_prouni_share  = median(prouni_share, na.rm = TRUE),
          p90_prouni_share  = quantile(prouni_share, 0.90, na.rm = TRUE)
        )
      write_csv(hei_fs_panel,
                "Tables_Final/tab_hei_first_stage_panel_v5.csv")
      write_csv(tab_hei_fs,
                "Tables_Final/tab_hei_first_stage_summary_v5.csv")
      message(sprintf("   HEI-level first stage: %d HEI-years across %d HEIs",
                      tab_hei_fs$n_hei_years, tab_hei_fs$n_ies))
    }

    # --- (B) drop microregions dominated by bottom-quartile-CPC HEIs ---
    cpc_q1 <- quantile(hei_quality_pre$cpc_pre, 0.25, na.rm = TRUE)
    bottom_ies <- hei_quality_pre %>%
      filter(!is.na(cpc_pre) & cpc_pre <= cpc_q1) %>%
      pull(id_ies)
    if (length(bottom_ies) > 0) {
      bottom_schol_by_micro <- df_prouni_raw %>%
        mutate(ano = as.integer(ano)) %>%
        inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                   by = "id_municipio") %>%
        group_by(id_microrregiao) %>%
        summarise(
          total_bol       = sum(n_bolsas, na.rm = TRUE),
          bottom_bol      = sum(n_bolsas[id_ies %in% bottom_ies], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(share_bottom = bottom_bol / pmax(total_bol, 1))

      dominated_micro <- bottom_schol_by_micro %>%
        filter(share_bottom >= 0.50) %>%
        pull(id_microrregiao)

      message(sprintf("   Dropping %d microregions dominated by bottom-CPC HEIs (share >= 50%%)",
                      length(dominated_micro)))

      data_low_hq <- data_low %>% filter(!(id_microrregiao %in% dominated_micro))
      res_low_hq  <- tryCatch(
        run_cs(data_low_hq, "log_wage_sm", label = "Low/HighQualityHEI"),
        error = function(e) { message(sprintf("   Low/HQ failed: %s", safe_err_msg(e))); NULL })

      # Complement: the dropped set by itself (expected to be near-zero / noisy)
      data_low_noq <- data_low %>% filter(id_microrregiao %in% dominated_micro)
      res_low_noq  <- tryCatch(
        run_cs(data_low_noq, "log_wage_sm", label = "Low/LowQualityHEI"),
        error = function(e) { message(sprintf("   Low/LQ failed: %s", safe_err_msg(e))); NULL })

      if (!is.null(res_low_hq)) {
        tab_att <- bind_rows(tab_att, build_att_row(res_low_hq))
        write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
      }
      if (!is.null(res_low_noq)) {
        tab_att <- bind_rows(tab_att, build_att_row(res_low_noq))
        write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
      }
    }
  }
}, error = function(e) {
  message(sprintf("   SECTION 9Y FAILED: %s", safe_err_msg(e)))
})


# ==============================================================================
# SECTION 9AA: DOSE-RESPONSE CURVE FIGURE
# ------------------------------------------------------------------------------
# The birth-cohort quadratic (Section 9W) identifies a concave dose-response
# relationship but there is no direct visualisation.  We build one here:
# x-axis = mean D_rt within each dose bin; y-axis = within-stratum ATT with
# 95% CI; overlay = fitted quadratic curve from fit_bc_1 (interpreted as
# the marginal dose-response on the wage gap).
# ==============================================================================
message(" SECTION 9AA: Dose-response curve figure...")

tryCatch({
  if (exists("res_low") && exists("res_high") &&
      !is.null(res_low) && !is.null(res_high)) {

    mean_d_bin <- df_final %>%
      filter(ano >= 2015, !is.na(dose_bin)) %>%
      group_by(dose_bin) %>%
      summarise(mean_D = mean(D_rt, na.rm = TRUE), .groups = "drop")
    md_lookup <- setNames(mean_d_bin$mean_D, as.character(mean_d_bin$dose_bin))

    bin_row <- function(bin_label, res_obj) {
      if (is.null(res_obj) || !(bin_label %in% names(md_lookup))) return(NULL)
      tibble(bin    = bin_label,
             mean_D = unname(md_lookup[bin_label]),
             att    = res_obj$simple$overall.att,
             se     = res_obj$simple$overall.se)
    }
    res_mid_obj <- if (exists("res_mid")) res_mid else NULL
    dose_resp_pts <- bind_rows(
      bin_row("Low",  res_low),
      bin_row("Mid",  res_mid_obj),
      bin_row("High", res_high)
    ) %>% mutate(lwr = att - Z_95 * se, upr = att + Z_95 * se)

    q_grid <- seq(0, max(df_final$D_rt, na.rm = TRUE), length.out = 200)
    q_fit  <- NULL
    if (exists("fit_bc_1") && !is.null(fit_bc_1)) {
      cf <- coef(fit_bc_1)
      b1 <- if ("D_rt"       %in% names(cf)) cf[["D_rt"]]       else NA_real_
      b2 <- if ("I(D_rt^2)"  %in% names(cf)) cf[["I(D_rt^2)"]]  else NA_real_
      if (!is.na(b1) && !is.na(b2))
        q_fit <- tibble(D = q_grid, y = b1*q_grid + b2*q_grid^2)
    }

    p_dr <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbar(data = dose_resp_pts,
                    aes(x = mean_D, ymin = lwr, ymax = upr, colour = bin),
                    width = 1, linewidth = 0.8) +
      geom_point(data = dose_resp_pts,
                 aes(x = mean_D, y = att, colour = bin), size = 4) +
      scale_colour_manual(name = "", values = DOSE_COLORS) +
      labs(title    = "Dose-Response: Wage ATT as a Function of ProUni Exposure Intensity",
           subtitle = "Bin-level CS ATT (points, 95% CI) with birth-cohort quadratic overlay",
           x = expression(paste("Mean lagged dose ", D[rt], " in 2015-2019")),
           y = "ATT on log wage")
    if (!is.null(q_fit)) {
      p_dr <- p_dr + geom_line(data = q_fit, aes(x = D, y = y),
                               linetype = "solid", linewidth = 0.8, colour = "black")
    }
    print(p_dr)
    ggsave("Figures_Final/fig_dose_response_v5.pdf", p_dr, width = 8, height = 5.5, dpi = 300)
    message("   Saved: Figures_Final/fig_dose_response_v5.pdf")
  }
}, error = function(e) message(sprintf("   9AA failed: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9AB: FISCAL BCR WITH SENSITIVITY BANDS
# ------------------------------------------------------------------------------
# Refines the point BCR from Section 9X by computing LOW and HIGH case bounds:
#   LOW case : ATT - 1.96 SE (upper end of Low-Dose SE), cost R$7,000/year
#   HIGH case: ATT + 1.96 SE,                         cost R$5,000/year
# ==============================================================================
message(" SECTION 9AB: Fiscal BCR with sensitivity bands...")

tryCatch({
  if (exists("res_low_dr") && !is.null(res_low_dr)) {
    att <- res_low_dr$simple$overall.att
    se  <- res_low_dr$simple$overall.se
    lo_att <- att - Z_95 * se
    hi_att <- att + Z_95 * se

    # Reuse data-driven wage / worker / scholarship magnitudes from Section 9X
    # (fiscal_tab) when available; fall back to conservative order-of-magnitude
    # placeholders otherwise, so 9X and 9AB report mutually consistent BCRs.
    if (exists("fiscal_tab") && is.data.frame(fiscal_tab) && nrow(fiscal_tab) >= 1) {
      mean_wage_brl_mo          <- fiscal_tab$mean_monthly_wage_brl[1]
      workers_per_micro         <- fiscal_tab$avg_n_formal_workers[1]
      scholarships_per_micro_yr <- fiscal_tab$avg_n_scholarships_yr[1]
      cost_lo <- fiscal_tab$cost_per_scholarship_lo[1]
      cost_hi <- fiscal_tab$cost_per_scholarship_hi[1]
      cost_pt <- 0.5 * (cost_lo + cost_hi)
    } else {
      mean_wage_brl_mo          <- 3000
      workers_per_micro         <- 5000
      scholarships_per_micro_yr <- 50
      cost_lo <- 5000; cost_hi <- 7000; cost_pt <- 6000
    }
    annual_wage_gain <- function(a) (exp(a) - 1) * mean_wage_brl_mo * 12

    gain_per_schol <- function(a) annual_wage_gain(a) * workers_per_micro /
                                   max(scholarships_per_micro_yr, 1)
    gps_low  <- gain_per_schol(lo_att)
    gps_pt   <- gain_per_schol(att)
    gps_high <- gain_per_schol(hi_att)

    # LOW-case pairs conservative ATT with highest cost; HIGH-case pairs
    # optimistic ATT with lowest cost -- joint worst/best-case framing.
    bcr_tab <- tibble(
      Case                 = c("LOW-case", "Point", "HIGH-case"),
      ATT                  = c(lo_att, att, hi_att),
      WageGainPerMicroYr   = c(gps_low, gps_pt, gps_high) * scholarships_per_micro_yr,
      GainPerScholarshipYr = c(gps_low, gps_pt, gps_high),
      CostPerScholarshipYr = c(cost_hi, cost_pt, cost_lo),
      BCR                  = round(c(gps_low / cost_hi,
                                     gps_pt  / cost_pt,
                                     gps_high/ cost_lo), 2)
    )
    write_csv(bcr_tab, "Tables_Final/tab_fiscal_bcr_ci_v5.csv")
    bcr_caption <- sprintf(
      "Benefit-cost ratio for the Low-Dose ATT with 95\\%% CI sensitivity bands.  Wage-gain denominator uses mean monthly wage R\\$%s, %s formal workers per Low-Dose microregion, and %s scholarships per microregion per year (from Section 9X fiscal\\_tab when available).",
      formatC(mean_wage_brl_mo,          format = "d", big.mark = ","),
      formatC(round(workers_per_micro),  format = "d", big.mark = ","),
      formatC(round(scholarships_per_micro_yr), format = "d", big.mark = ","))
    tryCatch(print(xtable(bcr_tab,
                          caption = bcr_caption,
                          label   = "tab:bcr_ci",
                          digits  = c(0,0,4,0,0,0,2)),
                   file = "Tables_Final/tab_fiscal_bcr_ci_v5.tex",
                   include.rownames = FALSE),
             error = function(e) NULL)
    print(bcr_tab)
  }
}, error = function(e) message(sprintf("   9AB failed: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9AC: COURSE-FIELD AND MODALITY HETEROGENEITY FROM PROUNI
# ------------------------------------------------------------------------------
# Uses the modalidade_ensino (presencial vs EAD) and nivel_bolsa (integral vs
# parcial) fields added to the ProUni SQL.  Aggregates ProUni scholarships
# within each stratum and recomputes the Low-Dose wage ATT separately when
# the microregion's scholarship share is dominated by one modality.  A
# divergence in ATT between presencial-dominated and EAD-dominated Low-Dose
# markets would indicate that the wage channel is specific to traditional
# classroom instruction.
# ==============================================================================
message(" SECTION 9AC: Modality (presencial/EAD) heterogeneity...")

res_low_presencial <- NULL
res_low_ead        <- NULL

tryCatch({
  if ("modalidade_ensino" %in% names(df_prouni_raw) &&
      any(!is.na(df_prouni_raw$modalidade_ensino))) {
    # EAD detector: basedosdados encodes modalidade_ensino as a numeric code
    # (1 = Presencial, 2 = EAD). Convert to numeric first; fall back to text
    # regex for any legacy/decoded string representations.
    is_ead <- function(x) {
      v <- suppressWarnings(as.numeric(as.character(x)))
      ifelse(!is.na(v), v == 2,
             grepl("\\bead\\b|dist.ncia|distancia|\\bremoto\\b|semipresencial",
                   tolower(as.character(x))))
    }
    mod_share <- df_prouni_raw %>%
      mutate(ead_flag = is_ead(modalidade_ensino)) %>%
      inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                 by = "id_municipio") %>%
      group_by(id_microrregiao) %>%
      summarise(
        total = sum(n_bolsas, na.rm = TRUE),
        ead   = sum(n_bolsas[ead_flag], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(total > 0) %>%
      mutate(ead_share = ead / total)

    presencial_micro <- mod_share %>% filter(ead_share <  0.25) %>% pull(id_microrregiao)
    ead_micro        <- mod_share %>% filter(ead_share >= 0.50) %>% pull(id_microrregiao)

    data_low_presencial <- data_low %>% filter(id_microrregiao %in% presencial_micro)
    data_low_ead        <- data_low %>% filter(id_microrregiao %in% ead_micro)

    message(sprintf("   Presencial-dominated Low microregions: %d | EAD-dominated: %d",
                    length(presencial_micro), length(ead_micro)))

    res_low_presencial <- tryCatch(
      run_cs(data_low_presencial, "log_wage_sm", label = "Low/Presencial"),
      error = function(e) { message(sprintf("   Low/Pres failed: %s", safe_err_msg(e))); NULL })
    res_low_ead        <- tryCatch(
      run_cs(data_low_ead,        "log_wage_sm", label = "Low/EAD"),
      error = function(e) { message(sprintf("   Low/EAD failed: %s", safe_err_msg(e))); NULL })

    for (r in Filter(Negate(is.null), list(res_low_presencial, res_low_ead))) {
      tab_att <- bind_rows(tab_att, build_att_row(r))
    }
    write_csv(tab_att, "Tables_Final/tab_att_v5.csv")
  } else {
    message("   SKIP: modalidade_ensino not present (legacy cache or missing field).")
  }
}, error = function(e) message(sprintf("   9AC failed: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9AD: PROUNI SCHOLARSHIP DEMOGRAPHIC BREAKDOWN
# ------------------------------------------------------------------------------
# The thesis claims (a) majority-female scholarship uptake and (b) significant
# non-white uptake consistent with the affirmative-action design.  These claims
# must be backed by actual shares computed from df_prouni_raw, not asserted.
#
# IMPORTANT: basedosdados ProUni codes differ from RAIS codes:
#   sexo:     1 = Feminino, 2 = Masculino  (OPPOSITE of RAIS)
#   raca_cor: 1=Amarela, 2=Branca, 3=Indígena, 4=Parda, 5=Preta
#   tipo_bolsa: 1=Integral, 2=Parcial (50%), 3=Complementar (25%)
#
# Outputs:
#   Tables_Final/tab_prouni_desc_v5.tex  -- LaTeX descriptive table
#   Tables_Final/tab_prouni_desc_v5.csv  -- machine-readable version
# ==============================================================================
message(" SECTION 9AD: ProUni scholarship demographic breakdown...")

tryCatch({
  req_cols <- c("sexo", "raca_cor", "tipo_bolsa", "n_bolsas", "ano")
  if (!all(req_cols %in% names(df_prouni_raw))) {
    missing_c <- setdiff(req_cols, names(df_prouni_raw))
    message(sprintf("   SKIP: missing columns in df_prouni_raw: %s",
                    paste(missing_c, collapse = ", ")))
  } else {

    # ------------------------------------------------------------------
    # 9AD.1  Overall totals and shares (pooled 2005-2019)
    # ------------------------------------------------------------------
    total_schol <- sum(df_prouni_raw$n_bolsas, na.rm = TRUE)

    # Female share (ProUni sexo 1 = Feminino)
    female_schol <- sum(df_prouni_raw$n_bolsas[
      suppressWarnings(as.numeric(as.character(df_prouni_raw$sexo))) == 1], na.rm = TRUE)

    # Non-white share (raca_cor 3=Indígena, 4=Parda, 5=Preta)
    nonwhite_schol <- sum(df_prouni_raw$n_bolsas[
      suppressWarnings(as.numeric(as.character(df_prouni_raw$raca_cor))) %in% c(3, 4, 5)],
      na.rm = TRUE)

    # Integral share (tipo_bolsa 1 = Integral)
    integral_schol <- sum(df_prouni_raw$n_bolsas[
      suppressWarnings(as.numeric(as.character(df_prouni_raw$tipo_bolsa))) == 1],
      na.rm = TRUE)

    pct_female    <- 100 * female_schol   / max(total_schol, 1)
    pct_nonwhite  <- 100 * nonwhite_schol / max(total_schol, 1)
    pct_integral  <- 100 * integral_schol / max(total_schol, 1)

    message(sprintf(
      "   ProUni 2005-2019: total=%s | female=%.1f%% | non-white=%.1f%% | integral=%.1f%%",
      format(total_schol, big.mark = ","), pct_female, pct_nonwhite, pct_integral))

    # ------------------------------------------------------------------
    # 9AD.2  Year-by-year female and non-white trend
    # ------------------------------------------------------------------
    trend_tab <- df_prouni_raw %>%
      mutate(
        sexo_n    = suppressWarnings(as.numeric(as.character(sexo))),
        raca_n    = suppressWarnings(as.numeric(as.character(raca_cor))),
        bolsa_n   = suppressWarnings(as.numeric(as.character(tipo_bolsa))),
        female    = sexo_n == 1,
        nonwhite  = raca_n %in% c(3, 4, 5),
        integral  = bolsa_n == 1
      ) %>%
      group_by(ano) %>%
      summarise(
        n_schol      = sum(n_bolsas, na.rm = TRUE),
        n_female     = sum(n_bolsas[female],   na.rm = TRUE),
        n_nonwhite   = sum(n_bolsas[nonwhite], na.rm = TRUE),
        n_integral   = sum(n_bolsas[integral], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pct_female   = 100 * n_female   / pmax(n_schol, 1),
        pct_nonwhite = 100 * n_nonwhite / pmax(n_schol, 1),
        pct_integral = 100 * n_integral / pmax(n_schol, 1)
      ) %>%
      arrange(ano)

    write_csv(trend_tab, "Tables_Final/tab_prouni_trend_v5.csv")

    # ------------------------------------------------------------------
    # 9AD.3  Summary table for LaTeX
    # ------------------------------------------------------------------
    overall_row <- tibble(
      Category      = "Overall (2005--2019)",
      `N scholarships` = format(total_schol, big.mark = ","),
      `Female (\\%)`   = sprintf("%.1f", pct_female),
      `Non-white (\\%)` = sprintf("%.1f", pct_nonwhite),
      `Integral (\\%)`  = sprintf("%.1f", pct_integral)
    )

    # Early (2005-2009) vs Late (2015-2019) split for trend evidence
    period_tab <- df_prouni_raw %>%
      mutate(
        period    = case_when(ano <= 2009 ~ "2005--2009",
                              ano >= 2015 ~ "2015--2019",
                              TRUE ~ NA_character_),
        sexo_n    = suppressWarnings(as.numeric(as.character(sexo))),
        raca_n    = suppressWarnings(as.numeric(as.character(raca_cor))),
        bolsa_n   = suppressWarnings(as.numeric(as.character(tipo_bolsa)))
      ) %>%
      filter(!is.na(period)) %>%
      group_by(period) %>%
      summarise(
        n_schol    = sum(n_bolsas, na.rm = TRUE),
        n_female   = sum(n_bolsas[sexo_n == 1],            na.rm = TRUE),
        n_nonwhite = sum(n_bolsas[raca_n %in% c(3,4,5)],  na.rm = TRUE),
        n_integral = sum(n_bolsas[bolsa_n == 1],           na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        `N scholarships`  = format(n_schol, big.mark = ","),
        `Female (\\%)`    = sprintf("%.1f", 100 * n_female   / pmax(n_schol, 1)),
        `Non-white (\\%)` = sprintf("%.1f", 100 * n_nonwhite / pmax(n_schol, 1)),
        `Integral (\\%)`  = sprintf("%.1f", 100 * n_integral / pmax(n_schol, 1))
      ) %>%
      rename(Category = period) %>%
      select(Category, `N scholarships`, `Female (\\%)`,
             `Non-white (\\%)`, `Integral (\\%)`)

    prouni_desc_tab <- bind_rows(period_tab, overall_row)
    write_csv(prouni_desc_tab, "Tables_Final/tab_prouni_desc_v5.csv")

    # LaTeX output
    lat_lines <- c(
      "\\begin{table}[htbp]",
      "\\centering",
      "\\caption{ProUni scholarship recipient characteristics, 2005--2019}",
      "\\label{tab:prouni_desc}",
      "\\begin{tabular}{lrrrr}",
      "\\toprule",
      "Period & N scholarships & Female (\\%) & Non-white (\\%) & Integral (\\%) \\\\",
      "\\midrule"
    )
    for (i in seq_len(nrow(prouni_desc_tab))) {
      r <- prouni_desc_tab[i, ]
      lat_lines <- c(lat_lines,
        sprintf("%s & %s & %s & %s & %s \\\\",
                r$Category, r$`N scholarships`,
                r$`Female (\\%)`, r$`Non-white (\\%)`, r$`Integral (\\%)`))
    }
    lat_lines <- c(lat_lines,
      "\\bottomrule",
      paste0("\\multicolumn{5}{l}{\\footnotesize{\\textit{Notes:} Non-white = Ind\\'igena + Parda + Preta ",
             "(ProUni raca\\_cor codes 3--5). Integral = full tuition; Parcial = 50\\%.",
             " Source: \\texttt{basedosdados.br\\_mec\\_prouni.microdados}.}}"),
      "\\end{tabular}",
      "\\end{table}"
    )
    writeLines(lat_lines, "Tables_Final/tab_prouni_desc_v5.tex")
    message("   -> Tables_Final/tab_prouni_desc_v5.tex written.")

    # ------------------------------------------------------------------
    # 9AD.4  Time-trend figure: female and non-white share over time
    # ------------------------------------------------------------------
    p_demog_trend <- trend_tab %>%
      select(ano, pct_female, pct_nonwhite) %>%
      pivot_longer(cols = c(pct_female, pct_nonwhite),
                   names_to = "group", values_to = "pct") %>%
      mutate(group = recode(group,
               pct_female   = "Female",
               pct_nonwhite = "Non-white (Ind\u00edgena + Parda + Preta)")) %>%
      ggplot(aes(x = ano, y = pct, colour = group, linetype = group)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      scale_colour_manual(values = c("Female" = "#1f77b4",
                                     "Non-white (Ind\u00edgena + Parda + Preta)" = "#d62728")) +
      labs(
        title    = "ProUni scholarship recipient composition, 2005--2019",
        subtitle = "Share of total annual scholarship grants",
        x = "Year", y = "Share (%)",
        colour = NULL, linetype = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave("Figures_Final/fig_prouni_demog_trend_v5.pdf",
           p_demog_trend, width = 7, height = 4.5)
    message("   -> Figures_Final/fig_prouni_demog_trend_v5.pdf written.")
  }
}, error = function(e) message(sprintf("   9AD failed: %s", safe_err_msg(e))))

# ==============================================================================
# SECTION 9AE: WEIGHTED SCHOLARSHIP EQUIVALENTS -- DOSE ROBUSTNESS
# ------------------------------------------------------------------------------
# The main dose D_rt counts all scholarships equally (integral + parcial).
# An integral scholarship covers 100% of tuition; parcial covers 50%.
# A "full-scholarship equivalent" (FSE) dose weights them by value:
#   FSE = integral * 1.0 + parcial * 0.5 + complementar * 0.25
#
# If the dose bin classification is insensitive to FSE weighting, the main
# results do not depend on this aggregation choice.
# Output: Tables_Final/tab_dose_robustness_v5.csv
#         Figures_Final/fig_dose_weighted_v5.pdf
# ==============================================================================
message(" SECTION 9AE: Weighted scholarship equivalents (dose robustness)...")

tryCatch({
  if (!("tipo_bolsa" %in% names(df_prouni_raw))) {
    message("   SKIP: tipo_bolsa not in df_prouni_raw.")
  } else {

    # Build FSE-weighted aggregation at microregion-year level
    df_prouni_fse <- df_prouni_raw %>%
      mutate(
        ano    = as.integer(ano),
        bolsa_n = suppressWarnings(as.numeric(as.character(tipo_bolsa))),
        fse_weight = case_when(
          bolsa_n == 1 ~ 1.00,   # Integral
          bolsa_n == 2 ~ 0.50,   # Parcial 50%
          bolsa_n == 3 ~ 0.25,   # Complementar 25%
          TRUE         ~ 1.00    # Unknown: treat as full
        ),
        n_fse = n_bolsas * fse_weight
      ) %>%
      inner_join(df_geo %>% distinct(id_municipio, id_microrregiao),
                 by = "id_municipio") %>%
      group_by(id_microrregiao, ano) %>%
      summarise(
        new_schol_raw = sum(n_bolsas, na.rm = TRUE),
        new_schol_fse = sum(n_fse,    na.rm = TRUE),
        .groups = "drop"
      )

    # Reconstruct lagged FSE density alongside original D_rt
    df_dose_compare <- df_panel %>%
      left_join(df_prouni_fse %>% select(id_microrregiao, ano, new_schol_fse),
                by = c("id_microrregiao", "ano")) %>%
      mutate(new_schol_fse = replace_na(new_schol_fse, 0)) %>%
      group_by(id_microrregiao) %>%
      arrange(ano) %>%
      mutate(
        cum_fse   = cumsum(new_schol_fse),
        D_rt_fse  = (dplyr::lag(cum_fse, params$lag_years, default = 0) / pop) * 1000
      ) %>%
      ungroup()

    # Compare 2019 cross-section: D_rt vs D_rt_fse
    compare_2019 <- df_dose_compare %>%
      filter(ano == max(params$years)) %>%
      select(id_microrregiao, D_rt, D_rt_fse)

    dose_cor <- cor(compare_2019$D_rt, compare_2019$D_rt_fse,
                    use = "complete.obs")

    # Bin classification under FSE dose using same quantile thresholds
    fse_thresholds <- quantile(compare_2019$D_rt_fse,
                               c(params$q_low, params$q_mid, params$q_high),
                               na.rm = TRUE)

    compare_2019 <- compare_2019 %>%
      left_join(df_final %>% distinct(id_microrregiao, dose_bin), by = "id_microrregiao") %>%
      mutate(
        dose_bin_fse = case_when(
          D_rt_fse >  fse_thresholds[1] & D_rt_fse <= fse_thresholds[2] ~ "Low",
          D_rt_fse >  fse_thresholds[2] & D_rt_fse <= fse_thresholds[3] ~ "Mid",
          D_rt_fse >  fse_thresholds[3]                                   ~ "High",
          TRUE ~ NA_character_
        )
      )

    # Reclassification rate (only among microregions classified in both)
    reclass_tab <- compare_2019 %>%
      filter(!is.na(dose_bin), !is.na(dose_bin_fse)) %>%
      summarise(
        n_total      = n(),
        n_unchanged  = sum(dose_bin == dose_bin_fse),
        n_changed    = sum(dose_bin != dose_bin_fse),
        pct_unchanged = 100 * mean(dose_bin == dose_bin_fse)
      )

    message(sprintf(
      "   FSE vs raw D_rt: correlation=%.4f | %d of %d microregions unchanged (%.1f%%)",
      dose_cor, reclass_tab$n_unchanged, reclass_tab$n_total, reclass_tab$pct_unchanged))

    # Save cross-tabulation
    cross_tab <- compare_2019 %>%
      filter(!is.na(dose_bin), !is.na(dose_bin_fse)) %>%
      count(dose_bin, dose_bin_fse) %>%
      arrange(dose_bin, dose_bin_fse)
    write_csv(cross_tab, "Tables_Final/tab_dose_robustness_v5.csv")

    # Scatter: D_rt vs D_rt_fse
    p_dose_wt <- compare_2019 %>%
      filter(!is.na(dose_bin)) %>%
      ggplot(aes(x = D_rt, y = D_rt_fse, colour = dose_bin)) +
      geom_point(alpha = 0.55, size = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = DOSE_COLORS, name = "Raw dose bin") +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
               label = sprintf("r = %.4f", dose_cor), size = 3.5) +
      labs(
        title    = "Dose measure robustness: raw vs FSE-weighted scholarships",
        subtitle = sprintf("FSE weights: Integral=1.0, Parcial=0.5 | %.1f%% of microregions bin-stable",
                           reclass_tab$pct_unchanged),
        x = "D_rt  (raw scholarship count, per 1,000 youth)",
        y = "D_rt_fse  (full-scholarship equivalents, per 1,000 youth)"
      ) +
      theme_bw(base_size = 11)

    ggsave("Figures_Final/fig_dose_weighted_v5.pdf",
           p_dose_wt, width = 6.5, height = 5)
    message("   -> Figures_Final/fig_dose_weighted_v5.pdf written.")
  }
}, error = function(e) message(sprintf("   9AE failed: %s", safe_err_msg(e))))
