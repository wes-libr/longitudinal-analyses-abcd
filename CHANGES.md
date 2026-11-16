# longitudinal.dev tutorial fixes, 2026-11-14

Revised R scripts are in `longitudinal_dev_R/tutorials/`. The originals are in
`_superseded/longitudinal_dev_R_2026-11-13/tutorials/`.
`longitudinal-dev_fixes_2026-11-14.patch` applies the same code changes to the site's markdown sources
(`git apply` from the root of github.com/beginDevSci/longitudinal-dev, commit 6ccd81a).

## 1. Random effects (LMM / GLMM): 5 tutorials, 8 model calls
`(1 | site:family_id:participant_id)` is one grouping factor in lme4/glmmTMB, with one level per
participant. It fits no site or family variance, despite the inline comment saying it does.
Replaced with the nested structure:

    (1 | site) + (1 | site:family_id) + (1 | site:family_id:participant_id)
    (1 | site) + (1 | site:family_id) + (1 + time | site:family_id:participant_id)   # random-slope version

Files: lmm-random-intercept, lmm-time-invariant-covariates (both models, so the interaction model
matches the main model), glmm (with the variance-component table now looked up by group name rather
than `vcov[1]`, plus family and site rows), glmm-binary (m1-m3), glmm-interactions.
In lmm-random-intercept, `predict(re.form = ~ (1 | site:family_id:participant_id))` is now `re.form = NULL`.

Not changed: lmm-random-slopes, lmm-time-varying-covariates and glmm-count use `(1 + time | participant_id)`,
which is what their comments describe. Moving them to the nested structure would be a consistency change, not a bug fix.

## 2. GEE clustering: 3 tutorials
gee, gee-count and gee-time-varying-covariate clustered on `participant_id`, which treats siblings and
twins as independent. They now cluster on `family_id`, with the data sorted by family before
`geeglm()` (required for contiguous clusters). The diagnostics tables report the number of families
(clusters) as well as the number of participants. Exchangeable within family is only an approximate
working correlation; the sandwich SEs remain valid.

## 3. Time coding: 7 tutorials
`time = as.numeric(session_id) - 1` returns factor positions, not years. Its spacing depends on which
levels the session factor carries. The published random-intercept text reads its slope as "per biannual
assessment", for example, and glmm-binary labels time 0-3 as Baseline / Year 2 / 4 / 6.
Now `time` is years since baseline, parsed from the label (`ses-02A` → 2), in lmm-random-intercept,
lmm-random-slopes, lmm-time-varying-covariates, glmm-binary, glmm-count and gee-count.
In lmm-time-invariant-covariates (biennial waves relabeled 0-3) it is `2 * (wave index)`.
glmm-binary's prediction grid and axis breaks are now taken from the data instead of hard-coded 0:3.
gee, gee-time-varying, glmm and glmm-interactions already use explicit factors or consecutive annual waves,
so their time coding is unchanged.

## Testing
`tools/test_fixes.R` runs on simulated ABCD-like data: 22 sites, 3,000 families (80/17/3% of size 1/2/3),
3,736 participants, waves at years 0, 2, 4, 6, and a session factor that includes non-annual levels.
- Time: the old coding gave 1, 5, 7, 9 for years 0, 2, 4, 6 on that factor; the new coding gives 0, 2, 4, 6.
- Old LMM: one grouping factor with 3,736 levels (= participants). Its participant variance of 71.9 absorbs site and family.
- New LMM recovers site 9.4, family 18.8, participant 44.2 and residual 26.1 (true 8 / 20 / 40 / 25). The random-slope
  version converges without warnings.
- glmmTMB binary with the nested random slopes fits, and the `VarCorr()$cond` extraction still works. glmer's group-name lookup works
  (the site variance can be singular for binary outcomes, which is expected).
- GEE: family-clustered robust SEs are about 3-4% larger for intercept and family-level terms.

## Still to do (needs real data)
The results tables, figures and interpretation text in the 11 changed tutorials must be regenerated. In particular:
- Slopes are now per year, not per wave.
- Variance components and ICCs now split across site, family and participant.
- GEE α is now a within-family correlation.
All tutorials still load release 6.0.
