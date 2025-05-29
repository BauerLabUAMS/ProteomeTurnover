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

# Load necessary libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())  # Apply clean publication-ready theme

# Read in the input CSV and clean column names by replacing spaces with underscores
df <- read_csv("chromatogram_library_grouped_proteins_ModifiedTransition_Results.csv", 
               guess_max = 100000) %>%
  {set_names(., gsub(" ", "_", names(.)))}

# Clean up protein identifiers
df$Protein <- gsub('sp\\|', '', df$Protein)  # Remove Swiss-Prot prefix
df$Protein <- gsub('\\|.*_HUMAN;', ';', df$Protein)  # Remove gene names from internal proteins
df$Protein <- gsub('\\|.*_HUMAN', '', df$Protein)  # Remove gene name from last protein in group
df$Protein <- gsub(';', '; ', df$Protein)  # Add space after semicolons for readability

# View cleaned data (optional interactive call)
df

# Convert Area to numeric and remove NA values
df1 <- df %>%
  mutate(Area = as.numeric(Area)) %>%
  filter(!is.na(Area))

#####################################################################################################
# Prepare fragment-level data for ratio analysis (SILAC light/heavy)
temp <- df1 %>%
  filter(Quantitative) %>%
  select(Protein, Peptide, Replicate, 
         Peptide_Modified_Sequence_Monoisotopic_Masses, 
         Isotope_Label_Type, Area) %>%
  mutate(Peptide_Modified_Sequence_Monoisotopic_Masses = 
           str_remove_all(Peptide_Modified_Sequence_Monoisotopic_Masses, "\\[\\+57]")) %>%  # Remove +57 mod (carbamidomethylation)
  group_by(Peptide, Protein, Replicate, Isotope_Label_Type) %>%
  summarize(Total_Area_Fragment = sum(Area, na.rm = TRUE)) %>%
  distinct(Replicate, Peptide, Isotope_Label_Type, Total_Area_Fragment, .keep_all = TRUE) %>%
  spread(Isotope_Label_Type, Total_Area_Fragment) %>%  # Spread into light and heavy columns
  mutate(
    total = ifelse(is.na(light) & !is.na(heavy), heavy, heavy + light),  # Default to heavy if light is missing
    fraction = ifelse(is.na(light), NA, light / total),  # Light over total
    has_fraction = !is.na(fraction)
  ) %>%
  ungroup()

#####################################################################################################
# Calculate percent heavy from total peptide intensity per replicate
df2 <- df1 %>%
  group_by(Peptide, Protein, Replicate, Isotope_Label_Type) %>%
  summarize(Intensity = sum(Area, na.rm = TRUE)) %>%
  group_by(Peptide, Protein, Replicate) %>%
  mutate(Percent_heavy = Intensity / sum(Intensity, na.rm = TRUE)) %>%
  filter(Isotope_Label_Type == "heavy") %>%
  ungroup()

# Same as df2, but only include peptides marked as Quantitative and split by Fragment_Ion
df2_only_quant <- df1 %>%
  filter(Quantitative) %>%
  group_by(Peptide, Protein, Replicate, Isotope_Label_Type, Fragment_Ion) %>%
  summarize(Intensity = sum(Area, na.rm = TRUE)) %>%
  group_by(Peptide, Protein, Replicate, Fragment_Ion) %>%
  mutate(Percent_heavy = Intensity / sum(Intensity, na.rm = TRUE)) %>%
  filter(Isotope_Label_Type == "heavy") %>%
  ungroup()

# Visualize number of peptides per protein (non-quantitative)
df2 %>%
  distinct(Peptide, Protein) %>%
  count(Protein, sort = TRUE) %>%
  ggplot(aes(x = n)) + 
  geom_bar() +
  xlim(c(0, 25))

# Visualize number of peptides per protein (quantitative only)
df2_only_quant %>%
  distinct(Peptide, Protein) %>%
  count(Protein, sort = TRUE) %>%
  ggplot(aes(x = n)) + 
  geom_bar() +
  xlim(c(0, 25))

# Boxplot of percent heavy per replicate (non-quantitative)
df2 %>%
  ggplot(aes(x = Replicate, y = Percent_heavy)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 12))

# Save high-resolution plot
ggsave("Ratio_boxplot.tiff", width = 8, height = 6, dpi = 600, compression = "lzw")

# Boxplot of percent heavy per replicate (quantitative only)
df2_only_quant %>%
  ggplot(aes(x = Replicate, y = Percent_heavy)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 12))

# Save high-resolution plot
ggsave("Ratio_boxplot_only_quantitative.tiff", width = 8, height = 6, dpi = 600, compression = "lzw")

# Export final quantitative percent heavy data
df2_only_quant %>%
  write_tsv("chromatogram_library_grouped_peptides_ratios_only_quant.tsv")

