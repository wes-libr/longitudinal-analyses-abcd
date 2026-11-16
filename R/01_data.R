# =============================================================================
# 01_data.R : build analysis panels from ABCD 7.0 with NBDCtools
# =============================================================================
# Produces (in cfg$dir_derived, never in the repo):
#   panel_img     one row per participant x imaging wave (QC-passed T1)
#   panel_core    one row per participant x core assessment (annual + mid-year)
#   events        binge-onset interval data, one row per participant
#   baseline      baseline covariates, one row per participant
# and (in the repo) aggregate count tables.
#
# NBDCtools functions used: get_latest_release_abcd(), get_dd_abcd(),
# get_sessions_abcd(), create_dataset_abcd(), filter_events_abcd(),
# filter_empty_rows(). Nonresponse codes are converted with value_to_na = TRUE,
# which calls transf_value_to_na() and keeps the code as a tagged NA; the
# codes are tabulated below before the factors are converted to numbers.

lp_log("01_data: building panels")

## ---- 1. Requested variables -------------------------------------------------
brain_vars <- unlist(cfg$brain, use.names = FALSE)
design     <- c(cfg$site, cfg$family, cfg$age, cfg$sex, cfg$edu, cfg$income, cfg$weight)
imaging    <- c(brain_vars, cfg$icv, cfg$qc_t1, cfg$scanner)
substance  <- c(cfg$binge_vars, cfg$alc_lifetime)
other      <- c(unname(cfg$pds), unname(cfg$pds_mean), cfg$ext, cfg$visit_dt, cfg$visit_type, cfg$birth_id)
# Regional thickness: every Desikan ROI the 7.0 dictionary lists for cfg$roi_table
roi_vars   <- tryCatch(grep(cfg$roi_pattern, get_dd_abcd(release = cfg$release, tables = cfg$roi_table)$name, value = TRUE),
                       error = function(e) character(0))
all_vars   <- unique(c(design, imaging, substance, other, roi_vars))
essential  <- unique(c(cfg$site, cfg$family, cfg$age, cfg$sex,
                       unlist(cfg$brain[cfg$primary]), cfg$qc_t1, cfg$binge_vars[1]))

## ---- 2. Release and dictionary checks (NBDCtools metadata) ------------------
latest <- tryCatch(get_latest_release_abcd(), error = function(e) NA)
if (!identical(latest, cfg$release))
  lp_log("  WARNING: NBDCtoolsData latest ABCD release is ", latest, ", config expects ",
         cfg$release, ". Update NBDCtoolsData (>= 5.0.0 for 7.0) or cfg$release.")

dd <- get_dd_abcd(release = cfg$release, vars = all_vars)
absent <- setdiff(all_vars, dd$name)
if (length(intersect(absent, essential)))
  stop("Essential variable(s) not in the ABCD ", cfg$release, " dictionary: ",
       paste(intersect(absent, essential), collapse = ", "), call. = FALSE)
if (length(absent)) {
  lp_log("  dropping ", length(absent), " variable(s) not in the ", cfg$release,
         " dictionary: ", paste(absent, collapse = ", "))
  all_vars <- setdiff(all_vars, absent)
}
safe_write_csv(dd |> filter(name %in% all_vars) |>
                 select(any_of(c("name", "table_name", "label", "type_var", "unit"))),
               "00_variables_used")

sess <- tryCatch(get_sessions_abcd(release = cfg$release), error = function(e) NULL)
if (!is.null(sess)) safe_write_csv(sess, "00_sessions_table")

## ---- 3. Load with NBDCtools -------------------------------------------------
# `release` deliberately left at its default: see config.R.
raw <- create_dataset_abcd(
  dir_data        = cfg$dir_data,
  vars            = all_vars,
  categ_to_factor = TRUE,
  add_labels      = TRUE,
  value_to_na     = TRUE            # nonresponse codes -> tagged NA
) |>
  filter_empty_rows()
lp_log("  create_dataset_abcd(): ", nrow(raw), " rows, ", ncol(raw), " columns")

## ---- 4. Nonresponse reasons (from the tagged NAs) ---------------------------
na_reason_table <- function(d) {
  if (!requireNamespace("haven", quietly = TRUE)) return(NULL)
  vars <- setdiff(names(d), c("participant_id", "session_id"))
  bind_rows(lapply(vars, function(v) {
    x <- d[[v]]
    tags <- tryCatch(haven::na_tag(unclass(x)), error = function(e) NULL)
    if (is.null(tags) || all(is.na(tags))) return(NULL)
    tibble(variable = v, code = tags[!is.na(tags)]) |> count(variable, code, name = "n")
  }))
}
reasons <- na_reason_table(raw)
if (!is.null(reasons) && nrow(reasons)) safe_write_csv(reasons, "00_nonresponse_codes", count_cols = "n")

## ---- 5. Rename, convert, derive ---------------------------------------------
pick <- function(v) if (length(v) && all(v %in% names(raw))) rowSums(sapply(v, function(z) num(raw[[z]])), na.rm = FALSE) else NA_real_

d <- tibble(
  id         = as.character(raw$participant_id),
  session_id = as.character(raw$session_id),
  site       = as.character(raw[[cfg$site]]),
  fam        = as.character(raw[[cfg$family]]),
  age        = num(raw[[cfg$age]]),
  sex        = num(raw[[cfg$sex]]),
  edu        = if (cfg$edu %in% names(raw)) num(raw[[cfg$edu]]) else NA_real_,
  income     = if (cfg$income %in% names(raw)) num(raw[[cfg$income]]) else NA_real_,
  wt         = if (cfg$weight %in% names(raw)) num(raw[[cfg$weight]]) else NA_real_,
  icv        = pick(cfg$icv),
  qc_t1      = num(raw[[cfg$qc_t1]]),
  scanner    = if (cfg$scanner %in% names(raw)) as.character(raw[[cfg$scanner]]) else NA_character_,
  alc_lt     = if (cfg$alc_lifetime %in% names(raw)) num(raw[[cfg$alc_lifetime]]) else NA_real_,
  ext        = if (cfg$ext %in% names(raw)) num(raw[[cfg$ext]]) else NA_real_,
  pds_f      = if (cfg$pds[["f"]] %in% names(raw)) num(raw[[cfg$pds[["f"]]]]) else NA_real_,
  pds_m      = if (cfg$pds[["m"]] %in% names(raw)) num(raw[[cfg$pds[["m"]]]]) else NA_real_,
  visit_date = if (cfg$visit_dt %in% names(raw)) as.Date(substr(as.character(raw[[cfg$visit_dt]]), 1, 10)) else as.Date(NA),
  remote     = if (cfg$visit_type %in% names(raw)) as.integer(num(raw[[cfg$visit_type]]) %in% c(2, 3)) else NA_integer_,
  birth_id   = if (cfg$birth_id %in% names(raw)) as.character(raw[[cfg$birth_id]]) else NA_character_,
  pdsm_f     = if (cfg$pds_mean[["f"]] %in% names(raw)) num(raw[[cfg$pds_mean[["f"]]]]) else NA_real_,
  pdsm_m     = if (cfg$pds_mean[["m"]] %in% names(raw)) num(raw[[cfg$pds_mean[["m"]]]]) else NA_real_
)
roi_vars <- intersect(roi_vars, names(raw))
roi_names <- paste0("roi_", sub("^mr_y_smri__thk__dsk__", "", sub("_mean$", "", roi_vars)))
for (k in seq_along(roi_vars)) d[[roi_names[k]]] <- num(raw[[roi_vars[k]]])
save_derived(tibble(roi = roi_names, variable = roi_vars), "roi_names")
for (b in names(cfg$brain)) d[[b]] <- pick(cfg$brain[[b]])
binge_cols <- intersect(cfg$binge_vars, names(raw))
bmat <- sapply(binge_cols, function(z) num(raw[[z]]))
if (is.null(dim(bmat))) bmat <- matrix(bmat, ncol = 1)
d$binge_days <- apply(bmat, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
d$binge_any  <- as.integer(d$binge_days > 0)

# Static fields arrive on the Static session (session_id = NA) in 7.0; carry
# them to every visit BEFORE any event filtering. Site and SES are taken from
# the baseline visit so that family is nested in one site for the whole panel.
d <- d |>
  arrange(id, is.na(session_id), session_id) |>
  group_by(id) |>
  mutate(across(c(fam, sex, wt, birth_id), ~ dplyr::first(na.omit(.x))),
         site_bl = dplyr::first(na.omit(site[session_id %in% "ses-00A"])),
         site_bl = ifelse(is.na(site_bl), dplyr::first(na.omit(site)), site_bl),
         edu_bl  = dplyr::first(na.omit(edu[session_id %in% "ses-00A"])),
         inc_bl  = dplyr::first(na.omit(income[session_id %in% "ses-00A"])),
         age_bl  = dplyr::first(na.omit(age[session_id %in% "ses-00A"]))) |>
  ungroup() |>
  mutate(site     = site_bl,
         female   = as.integer(sex == 2),
         pds      = ifelse(female == 1, pds_f, pds_m),
         pds_mean = ifelse(female == 1, pdsm_f, pdsm_m),
         # 7.0 codes highest caregiver education in 5 levels (1 = no diploma ...
         # 5 = graduate/professional), not years of schooling
         edu_cat  = factor(edu_bl, levels = 1:5,
                           labels = c("lt_HS", "HS_GED", "some_college", "bachelor", "graduate")),
         wave_yr  = years_from_session(session_id),
         age      = ifelse(is.na(age) & !is.na(age_bl) & !is.na(wave_yr), age_bl + wave_yr, age),
         time     = age - age_bl)            # years since baseline, from actual visit age

# Calendar period and birth cohort (age-period-cohort). Period is the decimal
# calendar year of the visit; cohort is the decimal birth year (period - age);
# pandemic exposure = years lived since cfg$covid_start at the visit.
dec_year <- function(x) as.numeric(format(x, "%Y")) + (as.numeric(format(x, "%j")) - 1) / 365.25
c0 <- dec_year(as.Date(cfg$covid_start)); c1 <- dec_year(as.Date(cfg$covid_end))
d <- d |>
  mutate(cal = dec_year(visit_date)) |>
  group_by(id) |>
  mutate(cal_bl = dplyr::first(na.omit(cal[session_id %in% "ses-00A"])),
         cal = ifelse(is.na(cal) & !is.na(cal_bl) & !is.na(time), cal_bl + time, cal)) |>
  ungroup() |>
  mutate(birth_year = cal - age,
         covid_era  = factor(ifelse(is.na(cal), NA, ifelse(cal < c0, "pre", ifelse(cal < c1, "restricted", "post"))),
                             levels = c("pre", "restricted", "post")),
         pand_yrs   = pmax(0, cal - c0))

## ---- 6. Panels via filter_events_abcd() ------------------------------------
core <- d |> filter(!is.na(session_id)) |>
  filter_events_abcd(c("annual", "mid_year"), connect = "or")
annual <- core |> filter_events_abcd("annual")

img <- d |> filter(!is.na(session_id)) |>
  filter_events_abcd(c("annual", "even")) |>          # imaging waves 00A/02A/04A/06A
  filter(qc_t1 == 1, !is.na(.data[[cfg$primary]])) |>
  mutate(is_bl = session_id == "ses-00A")
for (b in names(cfg$brain))  img[[paste0(b, "_z")]] <- std_by_baseline(img[[b]], img$is_bl)
img$icv_z <- std_by_baseline(img$icv, img$is_bl)
img$scanner[is.na(img$scanner)] <- "unknown"

## ---- 7. Binge onset (interval-censored at assessment ages) ------------------
events <- core |>
  filter(!is.na(binge_any), !is.na(time)) |>
  arrange(id, time) |>
  group_by(id) |>
  summarise(
    prevalent = binge_any[1] == 1,
    event     = as.integer(any(binge_any == 1)),
    t_right   = if (any(binge_any == 1)) time[which(binge_any == 1)[1]] else max(time),
    t_left    = if (any(binge_any == 1)) { k <- which(binge_any == 1)[1]; if (k > 1) time[k - 1] else 0 } else max(time),
    n_assess  = n(), .groups = "drop") |>
  filter(!prevalent) |>
  mutate(t_event = if (cfg$onset_time == "midpoint") ifelse(event == 1, (t_left + t_right) / 2, t_right) else t_right)

baseline <- d |> group_by(id) |>
  summarise(site = first(site), fam = first(fam), female = first(female),
            cal_bl = first(cal_bl), birth_year = first(na.omit(birth_year)), birth_id = first(birth_id),
            age_bl = first(age_bl), edu_bl = first(edu_bl), edu_cat = first(edu_cat),
            inc_bl = first(inc_bl), wt = first(wt), .groups = "drop")

events <- events |> inner_join(baseline, by = "id") |>
  filter(id %in% img$id, t_event > 0)

## ---- 8. Save and describe --------------------------------------------------
save_derived(img, "panel_img"); save_derived(core, "panel_core")
save_derived(annual, "panel_annual"); save_derived(events, "events"); save_derived(baseline, "baseline")

safe_write_csv(n_table(img, session_id), "01_n_imaging_by_wave", count_cols = c("n_obs", "n_id"))
safe_write_csv(core |> filter(!is.na(binge_any)) |> group_by(session_id) |>
                 summarise(n_assessed = n(), n_binge = sum(binge_any), .groups = "drop"),
               "01_binge_by_session", count_cols = c("n_assessed", "n_binge"))
safe_write_csv(tibble(n_participants = nrow(events), n_onsets = sum(events$event),
                      median_followup_yr = median(events$t_event)),
               "01_onset_summary", count_cols = c("n_participants", "n_onsets"))
lp_log("01_data done: ", n_distinct(img$id), " participants with imaging; ",
       sum(events$event), " binge onsets among ", nrow(events), " at risk")
