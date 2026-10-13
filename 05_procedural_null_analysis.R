# 05_procedural_null_analysis.R
#
# Compare empirical ROI-optimization gains with the procedural-null
# distributions.
#
# Empirical and simulated correlations are Fisher Z-transformed before
# optimization gains are calculated.
#
# The procedural-null comparison is descriptive:
#   exceedance count = number of simulated gains >= empirical mean gain



library(tidyverse)


# -------------------------------------------------------------------------
# Load empirical correlations
# -------------------------------------------------------------------------

full_run <- read_csv(
  "rtp_afni_correlations_full_run_revision_4mm.csv",
  show_col_types = FALSE
)

happy_rest <- read_csv(
  "rtp_afni_correlations_happyrest_4mm_rev_26_08_11.csv",
  show_col_types = FALSE
)


# -------------------------------------------------------------------------
# Put empirical Full Run, Happy, and Rest data in long format
# -------------------------------------------------------------------------

full_long <- full_run %>%

  transmute(
    ID,
    Target = ROI,
    Visit,
    Run,
    Condition = "FullRun",
    Original = Corr_Original,
    Static = Corr_Static,
    Dynamic = Corr_Dynamic
  ) %>%

  pivot_longer(
    cols = c(Original, Static, Dynamic),
    names_to = "Strategy",
    values_to = "Correlation"
  )


happy_long <- happy_rest %>%

  transmute(
    ID,
    Target = ROI,
    Visit,
    Run,
    Condition = "Happy",
    Original = Corr_Original_happy,
    Static = Corr_Static_happy,
    Dynamic = Corr_Dynamic_concat_happy
  ) %>%

  pivot_longer(
    cols = c(Original, Static, Dynamic),
    names_to = "Strategy",
    values_to = "Correlation"
  )


rest_long <- happy_rest %>%

  transmute(
    ID,
    Target = ROI,
    Visit,
    Run,
    Condition = "Rest",
    Original = Corr_Original_rest,
    Static = Corr_Static_rest,
    Dynamic = Corr_Dynamic_concat_rest
  ) %>%

  pivot_longer(
    cols = c(Original, Static, Dynamic),
    names_to = "Strategy",
    values_to = "Correlation"
  )


empirical <- bind_rows(
  full_long,
  happy_long,
  rest_long
)


# -------------------------------------------------------------------------
# Fisher Z-transform empirical run-level correlations
# -------------------------------------------------------------------------

empirical <- empirical %>%

  mutate(
    Correlation = pmin(
      pmax(Correlation, -0.9999999),
      0.9999999
    ),

    FisherZ = atanh(Correlation)
  )


# -------------------------------------------------------------------------
# Average runs within visit
# -------------------------------------------------------------------------

visit_level <- empirical %>%

  group_by(
    ID,
    Target,
    Visit,
    Condition,
    Strategy
  ) %>%

  summarise(
    MeanZ = mean(FisherZ, na.rm = TRUE),
    .groups = "drop"
  )


# -------------------------------------------------------------------------
# Average visits within participant
# -------------------------------------------------------------------------

participant_level <- visit_level %>%

  group_by(
    ID,
    Target,
    Condition,
    Strategy
  ) %>%

  summarise(
    MeanZ = mean(MeanZ, na.rm = TRUE),
    .groups = "drop"
  ) %>%

  mutate(
    ParticipantR = tanh(MeanZ)
  )


# -------------------------------------------------------------------------
# Calculate participant-level empirical gains
# -------------------------------------------------------------------------

empirical_gains <- participant_level %>%

  pivot_wider(
    names_from = Strategy,
    values_from = c(
      MeanZ,
      ParticipantR
    )
  ) %>%

  mutate(

    StaticGainZ =
      MeanZ_Static - MeanZ_Original,

    DynamicGainZ =
      MeanZ_Dynamic - MeanZ_Static,

    StaticGainR =
      ParticipantR_Static - ParticipantR_Original,

    DynamicGainR =
      ParticipantR_Dynamic - ParticipantR_Static
  )


# Mean empirical gain for each Target x Condition

empirical_mean <- empirical_gains %>%

  group_by(
    Target,
    Condition
  ) %>%

  summarise(

    N = n(),

    EmpiricalStaticGainZ =
      mean(StaticGainZ),

    EmpiricalDynamicGainZ =
      mean(DynamicGainZ),

    EmpiricalStaticGainR =
      mean(StaticGainR),

    EmpiricalDynamicGainR =
      mean(DynamicGainR),

    .groups = "drop"
  )


# =========================================================================
# Load procedural-null simulations
# =========================================================================

null_full <- read_csv(
  "procedural_null_fullrun.csv",
  show_col_types = FALSE
) %>%

  mutate(
    Condition = "FullRun"
  )


null_conditions <- read_csv(
  "procedural_null_conditions.csv",
  show_col_types = FALSE
)


null_data <- bind_rows(
  null_full,
  null_conditions
)


# -------------------------------------------------------------------------
# Fisher Z-transform simulated correlations
# -------------------------------------------------------------------------

null_data <- null_data %>%

  mutate(

    Z_Original =
      atanh(
        pmin(
          pmax(Corr_Original, -0.9999999),
          0.9999999
        )
      ),

    Z_Static =
      atanh(
        pmin(
          pmax(Corr_Static, -0.9999999),
          0.9999999
        )
      ),

    Z_Dynamic =
      atanh(
        pmin(
          pmax(Corr_Dynamic, -0.9999999),
          0.9999999
        )
      ),


    # Fisher-Z gains used for the null comparison

    StaticGainZ =
      Z_Static - Z_Original,

    DynamicGainZ =
      Z_Dynamic - Z_Static,


    # Raw-r gains retained for descriptive interpretation

    StaticGainR =
      Corr_Static - Corr_Original,

    DynamicGainR =
      Corr_Dynamic - Corr_Static
  )


# =========================================================================
# Static gain: Static - Original
# =========================================================================

static_results <- null_data %>%

  left_join(
    empirical_mean,
    by = c(
      "Target",
      "Condition"
    )
  ) %>%

  group_by(
    Target,
    Condition,
    N,
    EmpiricalStaticGainZ,
    EmpiricalStaticGainR
  ) %>%

  summarise(

    NNull = n(),

    NullMeanGainZ =
      mean(StaticGainZ),

    NullMeanGainR =
      mean(StaticGainR),

    ExceedanceCount =
      sum(
        StaticGainZ >= EmpiricalStaticGainZ
      ),

    ExceedanceProportion =
      ExceedanceCount / NNull,

    .groups = "drop"
  ) %>%

  mutate(
    Optimization = "Static - Original"
  )


# =========================================================================
# Dynamic gain: Dynamic - Static
# =========================================================================

dynamic_results <- null_data %>%

  left_join(
    empirical_mean,
    by = c(
      "Target",
      "Condition"
    )
  ) %>%

  group_by(
    Target,
    Condition,
    N,
    EmpiricalDynamicGainZ,
    EmpiricalDynamicGainR
  ) %>%

  summarise(

    NNull = n(),

    NullMeanGainZ =
      mean(DynamicGainZ),

    NullMeanGainR =
      mean(DynamicGainR),

    ExceedanceCount =
      sum(
        DynamicGainZ >= EmpiricalDynamicGainZ
      ),

    ExceedanceProportion =
      ExceedanceCount / NNull,

    .groups = "drop"
  ) %>%

  mutate(
    Optimization = "Dynamic - Static"
  )


# -------------------------------------------------------------------------
# Print results
# -------------------------------------------------------------------------

cat("\nSTATIC OPTIMIZATION\n")

print(
  static_results,
  n = Inf
)


cat("\nDYNAMIC OPTIMIZATION\n")

print(
  dynamic_results,
  n = Inf
)


# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

write_csv(
  static_results,
  "procedural_null_static_results.csv"
)

write_csv(
  dynamic_results,
  "procedural_null_dynamic_results.csv"
)
