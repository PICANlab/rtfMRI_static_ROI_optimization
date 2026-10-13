# 03_primary_statistical_analysis.R
#
# Primary statistical analysis of RTP-AFNI correspondence.
#
# R 4.5.0
# tidyverse 2.0.0
# afex 1.4.1
# emmeans 1.11.0


library(tidyverse)
library(afex)
library(emmeans)


# -------------------------------------------------------------------------
# Load run-level results
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
# Put Full Run, Happy, and Rest correlations in the same long-format table
# -------------------------------------------------------------------------

full_long <- full_run %>%

  select(
    ID, Group, Visit, Run, ROI,
    Corr_Original,
    Corr_Static,
    Corr_Dynamic
  ) %>%

  pivot_longer(
    cols = starts_with("Corr_"),
    names_to = "Strategy",
    values_to = "Correlation"
  ) %>%

  mutate(
    Strategy = str_remove(Strategy, "Corr_"),
    Condition = "FullRun"
  )


happy_long <- happy_rest %>%

  transmute(
    ID,
    Group,
    Visit,
    Run,
    ROI,
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
    Group,
    Visit,
    Run,
    ROI,
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


data_long <- bind_rows(
  full_long,
  happy_long,
  rest_long
) %>%

  mutate(

    Group = factor(
      Group,
      levels = c(0, 1),
      labels = c(
        "Parietal (HIPS)",
        "Amygdala"
      )
    ),

    Condition = factor(
      Condition,
      levels = c(
        "FullRun",
        "Happy",
        "Rest"
      )
    ),

    Strategy = factor(
      Strategy,
      levels = c(
        "Original",
        "Static",
        "Dynamic"
      )
    )
  )


# -------------------------------------------------------------------------
# Fisher Z-transform run-level correlations
# -------------------------------------------------------------------------

data_long <- data_long %>%

  mutate(

    Correlation = pmin(
      pmax(Correlation, -0.9999999),
      0.9999999
    ),

    FisherZ = atanh(Correlation)
  )


# -------------------------------------------------------------------------
# Average runs within each visit
# -------------------------------------------------------------------------

visit_level <- data_long %>%

  group_by(
    ID,
    Group,
    Visit,
    Condition,
    Strategy
  ) %>%

  summarise(
    VisitMeanZ = mean(FisherZ, na.rm = TRUE),
    .groups = "drop"
  )


# -------------------------------------------------------------------------
# Average visits within each participant
#
# Participant is the inferential unit.
# -------------------------------------------------------------------------

participant_level <- visit_level %>%

  group_by(
    ID,
    Group,
    Condition,
    Strategy
  ) %>%

  summarise(
    MeanZ = mean(VisitMeanZ, na.rm = TRUE),
    .groups = "drop"
  ) %>%

  mutate(
    ParticipantR = tanh(MeanZ)
  )


# -------------------------------------------------------------------------
# Descriptive Pearson-r values
#
# Fisher Z is used for inference.
# Back-transformed r is used for interpretation and figures.
# -------------------------------------------------------------------------

descriptives <- participant_level %>%

  group_by(
    Group,
    Condition,
    Strategy
  ) %>%

  summarise(
    N = n(),
    MeanR = mean(ParticipantR),
    SDR = sd(ParticipantR),
    .groups = "drop"
  )


print(descriptives)


# =========================================================================
# 1. Original ROI only
# =========================================================================

original_data <- participant_level %>%
  filter(Strategy == "Original")


anova_original <- aov_ez(
  id = "ID",
  dv = "MeanZ",
  data = original_data,
  between = "Group",
  within = "Condition",
  type = 3,
  anova_table = list(
    correction = "GG",
    es = "ges"
  )
)


cat("\nORIGINAL ROI: GROUP x CONDITION\n")

print(anova_original)


cat("\nGroup means:\n")

print(
  emmeans(
    anova_original,
    ~ Group
  )
)


cat("\nCondition contrasts:\n")

print(
  pairs(
    emmeans(
      anova_original,
      ~ Condition
    ),
    adjust = "tukey"
  )
)


# =========================================================================
# 2. Original versus Static
# =========================================================================

static_data <- participant_level %>%
  filter(Strategy %in% c("Original", "Static")) %>%
  droplevels()


anova_static <- aov_ez(
  id = "ID",
  dv = "MeanZ",
  data = static_data,
  between = "Group",
  within = c(
    "Strategy",
    "Condition"
  ),
  type = 3,
  anova_table = list(
    correction = "GG",
    es = "ges"
  )
)


cat("\nORIGINAL versus STATIC\n")

print(anova_static)


cat("\nOverall Original versus Static:\n")

print(
  pairs(
    emmeans(
      anova_static,
      ~ Strategy
    ),
    adjust = "tukey"
  )
)


cat("\nOriginal versus Static by condition:\n")

print(
  pairs(
    emmeans(
      anova_static,
      ~ Strategy | Condition
    ),
    adjust = "tukey"
  )
)


cat("\nOriginal versus Static by group:\n")

print(
  pairs(
    emmeans(
      anova_static,
      ~ Strategy | Group
    ),
    adjust = "tukey"
  )
)


# =========================================================================
# 3. Original versus Static versus Dynamic
# =========================================================================

anova_dynamic <- aov_ez(
  id = "ID",
  dv = "MeanZ",
  data = participant_level,
  between = "Group",
  within = c(
    "Strategy",
    "Condition"
  ),
  type = 3,
  anova_table = list(
    correction = "GG",
    es = "ges"
  )
)


cat("\nORIGINAL versus STATIC versus DYNAMIC\n")

print(anova_dynamic)


cat("\nOverall strategy contrasts:\n")

print(
  pairs(
    emmeans(
      anova_dynamic,
      ~ Strategy
    ),
    adjust = "tukey"
  )
)


cat("\nStrategy contrasts by condition:\n")

print(
  pairs(
    emmeans(
      anova_dynamic,
      ~ Strategy | Condition
    ),
    adjust = "tukey"
  )
)


cat("\nStrategy contrasts by group:\n")

print(
  pairs(
    emmeans(
      anova_dynamic,
      ~ Strategy | Group
    ),
    adjust = "tukey"
  )
)


cat("\nStrategy contrasts within Group x Condition:\n")

print(
  pairs(
    emmeans(
      anova_dynamic,
      ~ Strategy | Group * Condition
    ),
    adjust = "tukey"
  )
)


# =========================================================================
# 4. Participant-level optimization gains
# =========================================================================

gains <- participant_level %>%

  select(
    ID,
    Group,
    Condition,
    Strategy,
    MeanZ,
    ParticipantR
  ) %>%

  pivot_wider(
    names_from = Strategy,
    values_from = c(
      MeanZ,
      ParticipantR
    )
  ) %>%

  mutate(

    StaticMinusOriginal_Z =
      MeanZ_Static - MeanZ_Original,

    DynamicMinusStatic_Z =
      MeanZ_Dynamic - MeanZ_Static,

    StaticMinusOriginal_R =
      ParticipantR_Static - ParticipantR_Original,

    DynamicMinusStatic_R =
      ParticipantR_Dynamic - ParticipantR_Static
  )


gain_summary <- gains %>%

  group_by(
    Group,
    Condition
  ) %>%

  summarise(

    N = n(),

    MeanStaticGainZ =
      mean(StaticMinusOriginal_Z),

    MeanDynamicGainZ =
      mean(DynamicMinusStatic_Z),

    MeanStaticGainR =
      mean(StaticMinusOriginal_R),

    MeanDynamicGainR =
      mean(DynamicMinusStatic_R),

    .groups = "drop"
  )


cat("\nPARTICIPANT-LEVEL OPTIMIZATION GAINS\n")

print(gain_summary)
