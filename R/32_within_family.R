# =============================================================================
# 32_within_family.R : sibling and co-twin designs for the onset-anchored model
# =============================================================================
# Shared familial confounders (genes, home, neighbourhood) can make bingers
# differ from non-bingers in brain development. Comparing relatives removes
# everything constant within family.
#  W1 sibling design: participants in families with >= 2 members in the
#     imaging panel. Within-between (Mundlak) mixed model: each onset term
#     (eventual, post, tso) enters as its family mean plus its deviation from
#     the family mean; the deviation coefficients are within-family effects,
#     equivalent to family fixed effects.
#  W2 co-twin design: the same model restricted to twin/triplet sets (same
#     family and same birth id); within-pair effects also hold age and much of
#     the genetic background constant.
# Families discordant for onset carry the within-family information; their
# counts are reported (suppressed below 10).
lp_log("32_within_family")
if (!need("lme4")) stop("lme4 required")
d <- img_data()
fam_n <- d |> distinct(id, fam, birth_id) |> count(fam, name = "n_members")
d <- d |> left_join(fam_n, by = "fam") |> filter(n_members >= 2)
twin_sets <- d |> distinct(id, fam, birth_id) |> count(fam, birth_id, name = "n_birth") |> filter(n_birth >= 2)
d$twin <- paste(d$fam, d$birth_id) %in% paste(twin_sets$fam, twin_sets$birth_id)
fx <- list(); disc <- list()
mundlak <- function(dd, grp) dd |> group_by(.data[[grp]]) |>
  mutate(across(c(eventual, post, tso), list(fm = ~ mean(.x), w = ~ .x - mean(.x)))) |> ungroup()
for (des in c("W1_sibling", "W2_cotwin")) {
  dd <- if (des == "W2_cotwin") d |> filter(twin) |> mutate(set = paste(fam, birth_id)) else d |> mutate(set = fam)
  if (n_distinct(dd$set) < 10) { lp_log("  ", des, ": fewer than 10 sets, skipped"); next }
  pers <- dd |> distinct(set, id, eventual)
  disc[[des]] <- tibble(module = "32", model = des, n_sets = n_distinct(pers$set), n_id = n_distinct(pers$id),
                        n_discordant_sets = sum(tapply(pers$eventual, pers$set, function(x) length(unique(x)) > 1)))
  dd <- mundlak(dd, "set")
  for (m in names(cfg$brain)) {
    dm <- dd |> filter(!is.na(.data[[m]]))
    f <- paste0(m, " ~ time + eventual_w + eventual_fm + time:eventual + post_w + post_fm + tso_w + tso_fm + age_bl_c + female",
                icv_term(m), " + (1 | set) + (1 + time | id)", if (cfg$scanner_re) " + (1 | scanner)" else "")
    fit <- try_fit(paste(des, m), lme4::lmer(as.formula(f), data = dm, control = lme4::lmerControl(calc.derivs = FALSE)))
    if (is.null(fit)) next
    if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
    fx[[paste(des, m)]] <- tidy_fixed(fit, "32", des, m) |> filter(grepl("_w$|_fm$", term)) |>
      mutate(component = ifelse(grepl("_w$", term), "within_family", "between_family"),
             n_id = n_distinct(dm$id), n_sets = n_distinct(dm$set))
  }
}
safe_write_csv(bind_rows(disc), "32_family_sets", count_cols = c("n_sets", "n_id", "n_discordant_sets"))
out <- bind_rows(fx)
safe_write_csv(out, "32_within_family", count_cols = c("n_id", "n_sets"))
if (nrow(out)) key_write(out |> filter(outcome == cfg$primary, component == "within_family") |>
                           select(model, outcome, term, estimate, conf.low, conf.high, p.value, n_id), "32")
