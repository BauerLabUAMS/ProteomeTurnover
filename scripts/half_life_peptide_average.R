#   Copyright 2025 Michael A. Bauer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# Load required libraries
library(nlfitr)
library(broom)
library(qs)
library(modelr)
library(nplyr)

##################################################
# Preprocess input: calculate percent light, extract time, and group data
# Original input is from a peptide quantification file
df1 <- df2_only_quant %>% 
  # Mutate to calculate the percent of light-labeled peptide
  mutate(Percent_light = 1 - Percent_heavy) %>%
  mutate(time = as.numeric(timepoint)) %>% 
  group_by(PG.FastaHeaders, PG.Genes, condition, replicate) %>% 
  arrange(time) %>% 
  ungroup()

# Save intermediate dataframe
qsave(df1, "/results/percentlight_peptide.qs")

# Filter out timepoints outside the 0–48 hr window, and peptides that don't decrease in light labeling
df2 <- df1 %>%
  filter(time <= 48, time >= 0, !is.na(Percent_light)) %>%
  group_by(PG.Genes, peptide, PG.FastaHeaders, condition, replicate) %>%
  nest() %>%
  ungroup() %>%
  mutate(length = map_int(data, nrow)) %>%
  filter(length > 9) %>%  # Ensure enough timepoints (at least 10)
  mutate(
    t0 = data %>% map_dbl(~ slice(., which(time == 0))$Percent_light),
    t48 = data %>% map_dbl(~ slice(., which(time == 48))$Percent_light)
  ) %>%
  filter(t0 > t48) %>%  # Retain peptides with decreasing light over time
  select(-length, -t0, -t48)

# Fit decay models using a nonlinear least squares (nls) function
df3 <- df2 %>%  
  mutate(
    nls_model = map(data, safely(~ fitdecay1(x = time, y = Percent_light, data = .x, k = 0.2, ylo = 0, yhi = 1, weigh = TRUE))) %>% map("result"),
    glance_data = map(nls_model, safely(~ glance(.x))),
    augment_data = map(nls_model, safely(~ augment(.x))),
    tidy_data = map(nls_model, safely(~ tidy(.x))) %>% map("result")
  ) %>%
  filter(map_int(tidy_data, length) > 0)

df3 <- df2 %>%
  mutate(
    nls_model = future_pmap(data, ~ safely(fitdecay1(x = time, y = Percent_light, data = .x, k = 0.2, ylo = 0, yhi = 1, weigh = TRUE)) %>% map("result")),
    glance_data = future_map(nls_model, ~ safely(glance(.x))),
    augment_data = future_map(nls_model, ~ safely(augment(.x))),
    tidy_data = future_map(nls_model, ~ safely(tidy(.x))) %>% map("result")
  ) %>%
  filter(map_int(tidy_data, length) > 0)

# Calculate pseudo and adjusted R-squared for each model
df3 <- df3 %>% mutate(nlsR2 = map(nls_model, ~ nlsR2(nls_model = .x)))
df3 <- df3 %>%
  mutate(
    pseudo.R.squared = map_dbl(nlsR2, ~ .x[["pseudo.R.squared"]]),
    adj.R.squared = map_dbl(nlsR2, ~ .x[["adj.R.squared"]])
  )

###############################################################################################
# Filter peptides with high model fit (adjusted R² ≥ 0.9), unnest and summarize
df3_R9 <- df3 %>% 
  filter(adj.R.squared >= 0.9) %>% 
  select(PG.FastaHeaders, PG.Genes, peptide, condition, replicate, data) %>% 
  unnest(data)

# Calculate average Percent_light across replicates
df3_R9_avg <- df3_R9 %>% 
  group_by(PG.FastaHeaders, PG.Genes, condition, replicate, time) %>% 
  summarise(Percent_light_avg = mean(Percent_light, na.rm = TRUE)) %>%
  ungroup()

# Save average light ratios
df3_R9_avg %>%
  write_tsv("/results/average_Light_ratios_of_peptides.tsv")

# Re-nest averaged data for model fitting
df3_R9_avg_data <- df3_R9_avg %>%
  group_by(PG.FastaHeaders, PG.Genes, condition, replicate) %>%
  nest() %>%
  ungroup() %>%
  mutate(length = map_int(data, nrow)) %>%
  filter(length > 4) %>%
  select(-length)

# Fit decay models on averaged data
df3_R9_avg_data_model <- df3_R9_avg_data %>%
  mutate(
    nls_model = map(data, safely(~ fitdecay1(x = time, y = Percent_light_avg, data = .x, k = 0.2, ylo = 0, yhi = 1, weigh = TRUE))) %>% map("result"),
    augment_data = map(nls_model, augment),
    tidy_data = map(nls_model, safely(~ tidy(.x))) %>% map("result")
  ) %>%
  filter(map_int(tidy_data, length) > 0)

# Calculate model R-squared metrics
df3_R9_avg_data_model <- df3_R9_avg_data_model %>% mutate(nlsR2 = map(nls_model, ~ nlsR2(nls_model = .x)))
df3_R9_avg_data_model <- df3_R9_avg_data_model %>% 
  mutate(
    pseudo.R.squared = map_dbl(nlsR2, ~ .x[["pseudo.R.squared"]]),
    adj.R.squared = map_dbl(nlsR2, ~ .x[["adj.R.squared"]])
  )

# Extract estimated decay rates and compute half-lives
df4 <- df3_R9_avg_data_model %>%
  select(condition, replicate, PG.FastaHeaders, PG.Genes, tidy_data) %>%
  unnest(tidy_data) %>%
  select(condition, replicate, PG.FastaHeaders, PG.Genes, term, estimate) %>%
  spread(term, estimate)

df5 <- df4 %>%
  mutate(half_life = log(2) / k)

# Spread data by replicate-condition and filter for complete cases
df6_1 <- df5 %>%
  select(condition, replicate, PG.FastaHeaders, PG.Genes, half_life) %>%
  unite(Group, condition, replicate) %>%
  spread(Group, half_life)

df6_filter_1 <- df6_1 %>%
  filter(!(is.na(Acute_R1) | is.na(Acute_R2)) | !(is.na(Chronic_R1) | is.na(Chronic_R2)))

# Save half-life output
df6_filter_1 %>%
  write_tsv("../2023.10.17 SILAC AvE R1R2 2nd run/1st_run_protien_half_lifes_fitdecay_wieghted_refitted_average_ratios_v2_filter_all_timepoints.11.14.2023.tsv")

####################################################################################################
# Custom helper functions to calculate R² for nonlinear models

# Function to compute R² by comparing nls model vs null model
R2nls <- function(nls.obj) {
  if (class(nls.obj) != "nls") {
    stop("The object must be of class `nls`.")
  }
  da <- eval(nls.obj$data)
  resp.name <- all.vars(summary(nls.obj)$formula)[1]
  form <- paste(resp.name, "~1", sep = "")
  m0 <- stats::lm(form, da)
  an <- stats::anova(nls.obj, m0)
  sqn <- stats::deviance(nls.obj)
  sqe <- stats::deviance(m0)
  r2 <- 1 - (sqn / sqe)
  aov <- data.frame(
    " " = c("regression", "residuals"),
    Df = c(-an$Df[2], an$Res.Df[1]),
    `Sum Sq` = c(-an$Sum[2], an$Res.Sum[1]),
    `Mean Sq` = c(NA, NA),
    `F value` = c(NA, NA),
    `Pr(>F)` = c(NA, NA)
  )
  aov$`Mean Sq` <- aov$`Sum Sq` / aov$Df
  aov$`F value` <- c(aov$`Mean Sq`[1] / aov$`Mean Sq`[2], NA)
  aov$`Pr(>F)`[1] <- 1 - stats::pf(aov$`F value`[1], df1 = aov$Df[1], df2 = aov$Df[2])
  return(list(anova = aov, R2 = r2))
}

# Function to calculate pseudo and adjusted R² for nls models
nlsR2 <- function(nls_model = nls.obj){
  pred <- predict(nls_model)
  n <- length(pred)
  res <- resid(nls_model)
  w <- weights(nls_model)
  if (is.null(w)) w <- rep(1, n)
  rss <- sum(w * res ^ 2)
  resp <- pred + res
  center <- weighted.mean(resp, w)
  r.df <- summary(nls_model)$df[2]
  int.df <- 1
  tss <- sum(w * (resp - center)^2)
  r.sq <- 1 - rss / tss
  adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
  return(list(pseudo.R.squared = r.sq, adj.R.squared = adj.r.sq))
}

