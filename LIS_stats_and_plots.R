######## R code to perform statistical analyses and generate figures for the manuscript "Localized environmental heterogeneity enhances diversity in seasonally shaped microbial communities in an urban estuary".

### Load libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(indicspecies)

## Load Phyloseq objects
ps_18Sr <- readRDS("ps_rarefied18S")
ps_16Sr <- readRDS("ps_rarefied16S")
ps_18S <- readRDS("ps_18S")
ps_16S <- readRDS("ps_16S")

## Figure settings
# Libraries for plots
library(ggrepel)
library(tibble)
library(QsRutils)

# Set orders
orderSeason <- c("Fall", "Winter", "Spring", "Summer")
orderSites <- c("A4","B3", "D3", "AC", "CI", "RP", "CP", "SR")

site_shapes <- c(
  "AC" = 21,
  "CI"   = 22,
  "RP" = 23,
  "CP" = 24,
  "SR" = 25,
  "A4" = 3,
  "B3" = 4,
  "D3" = 8
)

# Define colors
season_colors <- c(
  "Fall"   = "#D55E00",
  "Winter" = "#0072B2",
  "Spring" = "#009E73",
  "Summer" = "E69F00"
)

barplots_colors <- c(
  '#40004b','#9970ab','#e7d4e8','#00441b', '#1b7837', '#a6dba0',
  '#662506','#bf812d', '#dfc27d','#01665e','#80cdc1','#c7eae5',
  '#c51b7d','#de77ae','#f1b6da','#fde0ef', '#993404', '#cc4c02')


diversity_colors <- c(
  "#E69F00",
  "#0072B2"
)

# THeme
theme_figures <- function() {
  
  theme_bw(base_size = 10) +
    
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.8),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}


### Fig 1. Environmental PCA

## Extract metadata
metadata <- data.frame(sample_data(ps_16S))
# Select variables for the PCA and standardize
env <- select_if(metadata, is.numeric)
env <- env %>%
  select(-c(Depth_m,Latitude,Longitude))
envStd <- decostand(env, method='standardize')

## PCA
pca_env <- prcomp(
  envStd,
  center = TRUE,
  scale. = FALSE
)

## Check contributions
# Rotation matrix
pca_env$rotation
# The first two components
pca_env$rotation[, 1:2]
# Sort by the absolute value of the PC1 loadings
sort(abs(pca_env$rotation[, 1]), decreasing = TRUE)
# Sort by the absolute value of the PC2 loadings
sort(abs(pca_env$rotation[, 2]), decreasing = TRUE)

## For plotting
# Extract scores
sitesPCA_env <- as.data.frame(pca_env$x)
sitesPCA_env$Season <- metadata$Season
sitesPCA_env$Site <- metadata$Site
# Sites order
sitesPCA_env$Site <- factor(
  sitesPCA_env$Site,
  levels = orderSites)
# Extract loadings
vars <- as.data.frame(pca_env$rotation)
vars$Variable <- rownames(vars)
# Percent variance explained
eig <- eigenvals(pca_env)
var_exp <- round(eig / sum(eig) * 100, 1)
pc1 <- var_exp[1]
pc2 <- var_exp[2]

# Plot PCA
PCAplot <- ggplot() +
  geom_point(
    data = sitesPCA_env,
    aes(
      PC1,
      PC2,
      colour = Season,
      shape = Site
    ),
    size = 4.5,
    stroke = 1.3
  ) +
  
  geom_segment(
    data = vars,
    aes(
      x = 0,
      y = 0,
      xend = PC1 * 2,
      yend = PC2 * 2
    ),
    arrow = arrow(length = unit(0.2,"cm")),
    colour = "grey30"
  ) +
  
  geom_text_repel(
    data = vars,
    aes(
      PC1 * 2,
      PC2 * 2,
      label = Variable
    ),
    size = 4
  ) +
  
  scale_colour_manual(values = season_colors) +
  
  scale_shape_manual(values = site_shapes) +
  
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)")
  ) +
  
  coord_equal() +
  
  theme_figures()

# Plot figure
PCAplot


### Fig 2. Microbial communities barplots

## Protistan community barplot data prep
# Convert to relative abundance and melt to long format
df_18S <- ps_18Sr %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt()
# Rename the Picozoa Division and Subdivision
df_18S$Subdivision <- ifelse(df_18S$Division == "Picozoa", "Picozoa", df_18S$Subdivision)
df_18S$Subdivision <- ifelse(df_18S$Subdivision == "Picozoa_subdivision_1", "Picozoa", df_18S$Subdivision)
# Calculate relative abundance of each Subdivision within each sample
df_18S <- df_18S %>%
  group_by(Sample, Subdivision) %>%
  mutate(
    taxon_abund_sample = sum(Abundance)
  ) %>%
  ungroup()

# Create plotting categories
# Define the diatom and dinoflagellate classes within Subdivision
diatom_classes <- c("Mediophyceae", "Coscinodiscophyceae", "Bacillariophyceae")
dino_class <- c("Dinophyceae")

df_18S <- df_18S %>%
  mutate(Plot_Subdivision = case_when(
    # Subdivisions representing <1% within each sample
    taxon_abund_sample < 0.01 ~ "Other (<1%)",
    # Gyrista: separate diatoms
    Subdivision == "Gyrista" & Class %in% diatom_classes ~ "Gyrista-Diatoms",
    Subdivision == "Gyrista" & !(Class %in% diatom_classes) ~ "Gyrista-Other",
    # Dinoflagellata: separate Dinophyceae
    Subdivision == "Dinoflagellata" & Class %in% dino_class ~ "Dinoflagellata-Dinophyceae",
    Subdivision == "Dinoflagellata" & !(Class %in% dino_class) ~ "Dinoflagellata-Syndiniales",
    # For all other use Subdivision
    TRUE ~ as.character(Subdivision)
  ))

## Bacterial-archaeal community barplot data prep
# Agglomerate to Order, convert to relative abundance and melt to long format
df_16S <- ps_16Sr %>%
  tax_glom(taxrank = "order") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt()
# Calculate relative abundance of each Phylum within each sample
df_16S <- df_16S %>%
  group_by(Sample, phylum) %>%
  mutate(
    phylum_abund_sample = sum(Abundance)
  ) %>%
  ungroup()

# Create plotting categories
df_16S <- df_16S %>%
  mutate(Plot_phylum = case_when(
    # Phylq representing <1% within each sample
    phylum_abund_sample < 0.01 ~ "Other (<1%)",
    # Separte Orders
    phylum == "Pseudomonadota" & order == "Pelagibacterales" ~ "Pseudomonadota-Pelagibacterales",
    phylum == "Pseudomonadota" & order == "Rhodobacterales" ~ "Pseudomonadota-Rhodobacterales",
    phylum == "Pseudomonadota" ~ "Pseudomonadota-Other",
    phylum == "Bacteroidota" & order == "Flavobacteriales" ~ "Bacteroidota-Flavobacteriales",
    phylum == "Bacteroidota" ~ "Bacteroidota-Other",
    phylum == "Actinomycetota" & order == "Actinomarinales" ~ "Actinomycetota-Actinomarinales",
    phylum == "Actinomycetota" ~ "Actinomycetota-Other",
    phylum == "Cyanobacteriota" & order == "Synechococcales" ~ "Cyanobacteriota-Synechococcales",
    phylum == "Cyanobacteriota" ~ "Cyanobacteriota-Other",
    # For all other use Phylum
    TRUE ~ as.character(phylum)
  ))

## Protistan community barplot
# Order the Subdivisions and sites
df_18S_plot <- df_18S
df_18S_plot$Plot_Subdivision <- df_18S_plot$Plot_Subdivision %>%
  fct_relevel(c("Gyrista-Diatoms", "Gyrista-Other"), after = 10) %>%
  fct_relevel("Other (<1%)", after = Inf)
ordered_sites_df_18S <- factor(df_18S_plot$Site, levels = orderSites)

barplot_18S <- ggplot(df_18S_plot, 
  aes(x = ordered_sites_df_18S, y = Abundance, fill = Plot_Subdivision)) +
  geom_col() +
  facet_grid(.~ factor(Season, levels = orderSeason), scales = "free_x", space = "free") +
  scale_fill_manual(values = barplots_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Site",
    y = "Relative Abundance",
    fill = NULL) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + 
  theme_figures()
  
## Bacterial-archaeal community barplot
# Order Phyla and sites
df_16S_plot <- df_16S
df_16S_plot$Plot_phylum <- df_16S_plot$Plot_phylum %>%
  fct_relevel(c("Cyanobacteriota-Synechococcales", "Cyanobacteriota-Other"), after = 9) %>%
  fct_relevel(c("Pseudomonadota-Pelagibacterales", "Pseudomonadota-Rhodobacterales", "Pseudomonadota-Other"), after = 12) %>%
  fct_relevel("Other (<1%)", after = Inf)
ordered_sites_df_16S <- factor(df_16S_plot$Site, levels = orderSites)

barplot_16S <- ggplot(df_16S_plot, 
                      aes(x = ordered_sites_df_16S, y = Abundance, fill = Plot_phylum)) +
  geom_col() +
  facet_grid(.~ factor(Season, levels = orderSeason), scales = "free_x", space = "free") +
  scale_fill_manual(values = barplots_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Site",
    y = "Relative Abundance",
    fill = NULL) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + 
  theme_figures()


### Fig 3. Apha-diversity (and Stats)

n_iterations <- 1000

## For 18S
min(sample_sums(ps_18S)) # 11995
rare_depth_18S <- 11995

alpha_iterations_18S <- lapply(seq_len(n_iterations), function(i) {
  ps_i <- rarefy_even_depth(
    ps_18S,
    sample.size = rare_depth_18S,
    rngseed = 20260922 + i,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
  alpha_i <- estimate_richness(
    ps_i,
    measures = c("Observed", "Shannon")
  )
  alpha_i %>%
    rownames_to_column("SampleID") %>%
    mutate(Iteration = i)
})

alpha_long_18S <- bind_rows(alpha_iterations_18S)

# Check that I have all the samples
setdiff(
  sample_names(ps_18S),
  unique(alpha_long_18S$SampleID)
)

sample_data_ps18S <- data.frame(sample_data(ps_18S)) %>%
  rownames_to_column("SampleID")

alpha_div_repeated_18S <- alpha_long_18S %>%
  group_by(SampleID) %>%
  summarise(
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_SD = sd(Observed, na.rm = TRUE),
    Observed_low = quantile(
      Observed,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Observed_high = quantile(
      Observed,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_SD = sd(Shannon, na.rm = TRUE),
    Shannon_low = quantile(
      Shannon,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Shannon_high = quantile(
      Shannon,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    
    Number_iterations = n(),
    .groups = "drop"
  ) %>%
  left_join(sample_data_ps18S, by = "SampleID")

# Check that the join worked:
nrow(alpha_div_repeated_18S)
nrow(sample_data_ps18S)
anyDuplicated(alpha_div_repeated_18S$SampleID)
any(is.na(sample_data_ps18S$Season))

# Retain the original names
alpha_div_repeated_18S <- alpha_div_repeated_18S %>%
  rename(
    Observed = Observed_mean,
    Shannon = Shannon_mean
  )


## For 16S
min(sample_sums(ps_16S)) # 16591
rare_depth_16S <- 16591

alpha_iterations_16S <- lapply(seq_len(n_iterations), function(i) {
  ps_i <- rarefy_even_depth(
    ps_16S,
    sample.size = rare_depth_16S,
    rngseed = 20260922 + i,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
  
  alpha_i <- estimate_richness(
    ps_i,
    measures = c("Observed", "Shannon")
  )
  
  alpha_i %>%
    rownames_to_column("SampleID") %>%
    mutate(Iteration = i)
})

alpha_long_16S <- bind_rows(alpha_iterations_16S)

# Check that I have all the samples
setdiff(
  sample_names(ps_16S),
  unique(alpha_long_16S$SampleID)
)

sample_data_ps16S <- data.frame(sample_data(ps_16S)) %>%
  rownames_to_column("SampleID")

alpha_div_repeated_16S <- alpha_long_16S %>%
  group_by(SampleID) %>%
  summarise(
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_SD = sd(Observed, na.rm = TRUE),
    Observed_low = quantile(
      Observed,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Observed_high = quantile(
      Observed,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_SD = sd(Shannon, na.rm = TRUE),
    Shannon_low = quantile(
      Shannon,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Shannon_high = quantile(
      Shannon,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    
    Number_iterations = n(),
    .groups = "drop"
  ) %>%
  left_join(sample_data_ps16S, by = "SampleID")

# Check that the join worked:
nrow(alpha_div_repeated_16S)
nrow(sample_data_ps16S)
anyDuplicated(alpha_div_repeated_16S$SampleID)
any(is.na(sample_data_ps16S$Season))

# Retain the original names
alpha_div_repeated_16S <- alpha_div_repeated_16S %>%
  rename(
    Observed = Observed_mean,
    Shannon = Shannon_mean
  )

#### Prepare for plots
alpha_div_repeated_18S <- alpha_div_repeated_18S %>%
  mutate(Group = "18S")

alpha_div_repeated_16S <- alpha_div_repeated_16S %>%
  mutate(Group = "16S")

all_df <- bind_rows(
  alpha_div_repeated_18S,
  alpha_div_repeated_16S
)

all_df %>%
  count(Group)
all_df %>%
  count(Group, SampleID) %>%
  filter(n != 1)


all_df$Site <- factor(all_df$Site,
                      levels = orderSites,
                      ordered = TRUE)

theme_figures <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.8),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Make a figure function
make_alpha_plot <- function(data, group_name, metric, y_label,
                            x_label = NULL) {
  
  ggplot(
    filter(data, Group == group_name),
    aes(
      x = Site,
      y = .data[[metric]],
      fill = Habitat
    )
  ) +
    geom_boxplot(
      staplewidth = 0.5,
      linewidth = 0.6,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitter(
        width = 0.10,
        height = 0,
        seed = 123
      ),
      shape = 21,
      size = 1.7,
      stroke = 0.3,
      colour = "black"
    ) +
    scale_fill_manual(
      values = diversity_colors,
      drop = FALSE
    ) +
    labs(
      x = x_label,
      y = y_label,
      fill = "Habitat"
    ) +
    theme_figures()
}

div_fig1 <- make_alpha_plot(
  all_df,
  group_name = "18S",
  metric = "Observed",
  y_label = "Observed ASV Richness"
)

div_fig2 <- make_alpha_plot(
  all_df,
  group_name = "18S",
  metric = "Shannon",
  y_label = "Shannon Diversity"
)

div_fig3 <- make_alpha_plot(
  all_df,
  group_name = "16S",
  metric = "Observed",
  y_label = "Observed ASV Richness",
  x_label = "Site"
)

div_fig4 <- make_alpha_plot(
  all_df,
  group_name = "16S",
  metric = "Shannon",
  y_label = "Shannon Diversity",
  x_label = "Site"
)

## Apha-diversity mixed effects model
library(lmerTest)

fit_alpha_model <- function(data, group_name, metric_name) {
  
  dat <- data %>%
    filter(Group == group_name) %>%
    filter(
      !is.na(.data[[metric_name]]),
      !is.na(Habitat),
      !is.na(Season),
      !is.na(Site)
    )
  
  if (metric_name == "Observed") {
    
    model_formula <- log(Observed) ~
      Habitat + Season + (1 | Site)
    
  } else {
    
    model_formula <- Shannon ~
      Habitat + Season + (1 | Site)
  }
  
  model <- lmer(
    model_formula,
    data = dat,
    REML = TRUE
  )
  
  model_anova <- anova(
    model,
    ddf = "Satterthwaite"
  )
  
  habitat_test <- data.frame(
    Group = group_name,
    Metric = metric_name,
    F_value = unname(
      model_anova["Habitat", "F value"]
    ),
    Numerator_df = unname(
      model_anova["Habitat", "NumDF"]
    ),
    Denominator_df = unname(
      model_anova["Habitat", "DenDF"]
    ),
    P_value = unname(
      model_anova["Habitat", "Pr(>F)"]
    )
  )
  
  list(
    model = model,
    anova = model_anova,
    habitat_test = habitat_test
  )
}

fit_18S_observed <- fit_alpha_model(
  all_df, "18S", "Observed"
)

fit_18S_shannon <- fit_alpha_model(
  all_df, "18S", "Shannon"
)

fit_16S_observed <- fit_alpha_model(
  all_df, "16S", "Observed"
)

fit_16S_shannon <- fit_alpha_model(
  all_df, "16S", "Shannon"
)

test_results <- bind_rows(
  fit_18S_observed$habitat_test,
  fit_18S_shannon$habitat_test,
  fit_16S_observed$habitat_test,
  fit_16S_shannon$habitat_test
) %>%
  mutate(
    P_Holm = p.adjust(
      P_value,
      method = "holm"
    )
  )

test_results


### Fig 4. NMDS and dbRDA (and Stats)

## Function to obtain averaged Bray-Curtis matrices
average_bray <- function(ps, iterations = 1000, seed = 2001) {
  
  taxa <- as(otu_table(ps), "matrix")
  
  # avgdist() requires samples in rows
  if (taxa_are_rows(ps)) {
    taxa <- t(taxa)
  }
  
  depth <- min(rowSums(taxa))
  
  message("Rarefaction depth: ", depth)
  
  set.seed(seed)
  
  avgdist(
    taxa,
    sample = depth,
    iterations = iterations,
    dmethod = "bray",
    diag = FALSE,
    upper = FALSE
  )
}

dist_Bray18S <- average_bray(
  ps_18S
)

dist_Bray16S <- average_bray(
  ps_16S
)

# Check the results
dist_Bray18S
dist_Bray16S

labels(dist_Bray18S)
labels(dist_Bray16S)

setequal(labels(dist_Bray18S), sample_names(ps_18S))
setequal(labels(dist_Bray16S), sample_names(ps_16S))

## Save distance matrices for subsequent analyses
bray_results <- list(
  dist_Bray18S = dist_Bray18S,
  dist_Bray16S = dist_Bray16S,
  depth_18S = min(sample_sums(ps_18S)),
  depth_16S = min(sample_sums(ps_16S)),
  iterations = 1000,
  seed_18S = 2001,
  seed_16S = 2002
)

saveRDS(
  bray_results,
  file = "average_Bray_distances_18S_16S.rds"
)

## NMDS 18S
set.seed(2001)
NMDS_18S <- ordinate(
  ps_18S,
  method = "NMDS",
  distance = dist_Bray18S,
  autotransform = FALSE,
  trymax = 500
)

NMDS_18S$stress

## PERMANOVA 18S
# Extract metadata
metadata_18S <- data.frame(
  sample_data(ps_18S),
  check.names = FALSE
)

# Make sure metadata is in the same order of the Bray–Curtis matrix
metadata_18S <- metadata_18S[
  labels(dist_Bray18S),
  ,
  drop = FALSE
]

# Confirm correspondence
stopifnot(
  identical(
    rownames(metadata_18S),
    labels(dist_Bray18S)
  )
)

# Adonis test
set.seed(123)
test_adonis_Season_18S <- adonis2(
  dist_Bray18S ~ Season, 
  data = metadata_18S, 
  permutations = 9999
)
test_adonis_Season_18S

set.seed(123)
test_adonis_Habitat_18S <- adonis2(
  dist_Bray18S ~ Habitat, 
  data = metadata_18S, 
  permutations = 9999
)
test_adonis_Habitat_18S

set.seed(123)
test_adonis_SeasonHabitat_18S <- adonis2(
  dist_Bray18S ~ Season + Habitat, 
  data = metadata_18S, 
  permutations = 9999
)
test_adonis_SeasonHabitat_18S

set.seed(123)
test.adonis_Site_18S <- adonis2(
  dist_Bray18S ~ Site, 
  data = metadata_18S, 
  permutations = 9999
)
test.adonis_Site_18S

# Homogeneity of dispersion test
beta_Season_18SS <- betadisper(
  dist_Bray18S, 
  metadata_18S$Season
)
set.seed(100)
permutest(beta_Season_18SS, 
          permutations = 9999
)

beta_Habitat_18SS <- betadisper(
  dist_Bray18S, 
  metadata_18S$Habitat
)
set.seed(100)
permutest(beta_Habitat_18SS, 
          permutations = 9999
)

beta_Site_18SS <- betadisper(
  dist_Bray18S, 
  metadata_18S$Site
)
set.seed(100)
permutest(beta_Site_18SS, 
          permutations = 9999
)

# Check dispersion distances
beta_Season_18SS$group.distances

beta_Habitat_18SS$group.distances

## NMDS 16S
set.seed(2001)
NMDS_16S <- ordinate(
  ps_16S,
  method = "NMDS",
  distance = dist_Bray16S,
  autotransform = FALSE
)

NMDS_16S$stress

## PERMANOVA 16S
# Extract metadata
metadata_16S <- data.frame(
  sample_data(ps_16S),
  check.names = FALSE
)

# Make sure metadata is in the same order of the Bray–Curtis matrix
metadata_16S <- metadata_16S[
  labels(dist_Bray16S),
  ,
  drop = FALSE
]

# Confirm correspondence
stopifnot(
  identical(
    rownames(metadata_16S),
    labels(dist_Bray16S)
  )
)

# Adonis tests
set.seed(123)
test_adonis_Season_16S <- adonis2(
  dist_Bray16S ~ Season, 
  data = metadata_16S, 
  permutations = 9999
)
test_adonis_Season_16S

set.seed(123)
test_adonis_Habitat_16S <- adonis2(
  dist_Bray16S ~ Habitat, 
  data = metadata_16S, 
  permutations = 9999
)
test_adonis_Habitat_16S

set.seed(123)
test_adonis_SeasonHabitat_16S <- adonis2(
  dist_Bray16S ~ Season + Habitat, 
  data = metadata_16S, 
  permutations = 9999
)
test_adonis_SeasonHabitat_16S

set.seed(123)
test.adonis_Site_16S <- adonis2(dist_Bray16S ~ Site, 
                                data = metadata_16S, 
                                permutations = 9999
)
test.adonis_Site_16S

# Homogeneity of dispersion test
beta_Season_16SS <- betadisper(
  dist_Bray16S, 
  metadata_16S$Season
)
set.seed(100)
permutest(beta_Season_16SS, 
          permutations = 9999
)

beta_Habitat_16SS <- betadisper(
  dist_Bray16S, 
  metadata_16S$Habitat
)
set.seed(100)
permutest(beta_Habitat_16SS, 
          permutations = 9999
)

beta_Site_16SS <- betadisper(
  dist_Bray16S, 
  metadata_16S$Site
)
set.seed(100)
permutest(beta_Site_16SS, 
          permutations = 9999
)

# Check dispersion distances
beta_Season_16SS$group.distances

beta_Habitat_16SS$group.distances

## dbRDA 18S
# Extract environmental data and standardize
env18S <- data.frame(sample_data(ps_18S))
env_con18S_std <- env18S %>%
  select(where(is.numeric)) %>%
  select(-c(Depth_m, Latitude, Longitude)) %>%
  decostand(method = "standardize")
env_con18S_std <- env_con18S_std[
  labels(dist_Bray18S),, drop = FALSE]

# dbRDA
dbRDA_18S <- dbrda(
  dist_Bray18S ~ .,
  data = env_con18S_std,
)

dbRDA_18S
summary(dbRDA_18S)

# Overall model
set.seed(2001)
anova(dbRDA_18S, 
      permutations = 9999
)

# Constrained axes
set.seed(2001)
anova(
  dbRDA_18S,
  by = "axis",
  permutations = 9999
)

RsquareAdj(dbRDA_18S)
# $r.squared
# [1] 0.6500681
# 
# $adj.r.squared
# [1] 0.4001167

# VIF
sort(vif.cca(dbRDA_18S))

set.seed(2001)
envfit_18S <- envfit(
  dbRDA_18S,
  env_con18S_std,
  permutations = 9999
)
envfit_18S

## dbRDA 16S
# Extract environmental data and standardize
env16S <- data.frame(sample_data(ps_16S))
env_con16S_std <- env16S %>%
  select(where(is.numeric)) %>%
  select(-c(Depth_m, Latitude, Longitude)) %>%
  decostand(method = "standardize")
env_con16S_std <- env_con16S_std[
  labels(dist_Bray16S),, drop = FALSE]

# dbRDA 
dbRDA_16S <- dbrda(
  dist_Bray16S ~ .,
  data = env_con16S_std
)

dbRDA_16S
summary(dbRDA_16S)

# Overall model
set.seed(2001)
anova(dbRDA_16S, 
      permutations = 9999
)

# Constrained axes
set.seed(2001)
anova(
  dbRDA_16S,
  by = "axis",
  permutations = 9999
)

RsquareAdj(dbRDA_16S)
# $r.squared
# [1] 0.6961018
# 
# $adj.r.squared
# [1] 0.5173381

# VIF
sort(vif.cca(dbRDA_16S))

set.seed(2001)
envfit_16S <- envfit(
  dbRDA_16S,
  env_con16S_std,
  permutations = 9999
)
envfit_16S

## Prepare for plots
# Set orders as factor levels
# 18S
sample_data(ps_18S)$Season <- factor(
  sample_data(ps_18S)$Season,
  levels = orderSeason)

sample_data(ps_18S)$Site <- factor(
  sample_data(ps_18S)$Site,
  levels = orderSites)

# 16S
sample_data(ps_16S)$Season <- factor(
  sample_data(ps_16S)$Season,
  levels = orderSeason)

sample_data(ps_16S)$Site <- factor(
  sample_data(ps_16S)$Site,
  levels = orderSites)

## Plot NMDS 18S
# Extract site scores
NMDS_scores_18S <- as.data.frame(
  scores(NMDS_18S, display = "sites")
) %>%
  rownames_to_column("Sample")

# Extract metadata for plotting
metadata_plot_18S <- metadata_18S %>%
  rownames_to_column("Sample")

# Combine scores and metadata
NMDS_df_18S <- NMDS_scores_18S %>%
  left_join(metadata_plot_18S, by = "Sample")

lim_NMDS_18S <- 1.05 * max(
  abs(c(
    NMDS_df_18S$NMDS1,
    NMDS_df_18S$NMDS2
  )),
  na.rm = TRUE
)

NMDS_plot_18S <- ggplot(
  NMDS_df_18S,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Season,
    shape = Site
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    stroke = 1.3
  ) +
  scale_color_manual(
    values = season_colors,
    breaks = orderSeason,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    breaks = orderSites,
    drop = FALSE
  ) +
  labs(
    color = "Season",
    shape = "Site"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-lim_NMDS_18S, lim_NMDS_18S),
    ylim = c(-lim_NMDS_18S, lim_NMDS_18S)
  ) +
  theme_figures()

## Plot NMDS 16S
# Extract site scores
NMDS_scores_16S <- as.data.frame(
  scores(NMDS_16S, display = "sites")
) %>%
  rownames_to_column("Sample")

# Extract metadata for plotting
metadata_plot_16S <- metadata_16S %>%
  rownames_to_column("Sample")

# Combine scores and metadata
NMDS_df_16S <- NMDS_scores_16S %>%
  left_join(metadata_plot_16S, by = "Sample")

lim_NMDS_16S <- 1.05 * max(
  abs(c(
    NMDS_df_16S$NMDS1,
    NMDS_df_16S$NMDS2
  )),
  na.rm = TRUE
)

NMDS_plot_16S <- ggplot(
  NMDS_df_16S,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Season,
    shape = Site
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    stroke = 1.3
  ) +
  scale_color_manual(
    values = season_colors,
    breaks = orderSeason,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    breaks = orderSites,
    drop = FALSE
  ) +
  labs(
    color = "Season",
    shape = "Site"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-lim_NMDS_16S, lim_NMDS_16S),
    ylim = c(-lim_NMDS_16S, lim_NMDS_16S)
  ) +
  theme_figures()

## Plot dbRDA 18S
# Extract environmental vector scores
vectors_18S <- as.data.frame(
  scores(
    envfit_18S,
    display = "vectors",
    choices = 1:2
  )
) %>%
  rownames_to_column("Variable")

# Add envfit p values
vectors_18S$p_value <- envfit_18S$vectors$pvals

# Retain significant vectors
vectors_18S <- vectors_18S %>%
  filter(p_value < 0.05)

sites_18S <- as.data.frame(
  scores(
    dbRDA_18S,
    display = "sites",
    choices = 1:2
  )
) %>%
  rownames_to_column("Sample") %>%
  left_join(
    env18S %>%
      rownames_to_column("Sample"),
    by = "Sample"
  )

# Scale arrows for plotting
if (nrow(vectors_18S) > 0) {
  arrow_mult_18S <- 0.8 * min(
    max(abs(sites_18S$dbRDA1), na.rm = TRUE) /
      max(abs(vectors_18S$dbRDA1), na.rm = TRUE),
    max(abs(sites_18S$dbRDA2), na.rm = TRUE) /
      max(abs(vectors_18S$dbRDA2), na.rm = TRUE)
  )
  vectors_18S <- vectors_18S %>%
    mutate(
      dbRDA1 = dbRDA1 * arrow_mult_18S,
      dbRDA2 = dbRDA2 * arrow_mult_18S,
      label_x = dbRDA1 * 1.10,
      label_y = dbRDA2 * 1.10
    )
}

lim_dbRDA_18S <- 1.10 * max(
  abs(c(
    sites_18S$dbRDA1,
    sites_18S$dbRDA2,
    vectors_18S$dbRDA1,
    vectors_18S$dbRDA2
  )),
  na.rm = TRUE
)

dbRDA_plot_18S <- ggplot(
  sites_18S,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Season,
    shape = Site
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    stroke = 1.3
  ) +
  geom_segment(
    data = vectors_18S,
    aes(
      x = 0,
      y = 0,
      xend = dbRDA1,
      yend = dbRDA2
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "darkgray",
    arrow = arrow(
      length = unit(0.18, "cm")
    )
  ) +
  geom_text_repel(
    data = vectors_18S,
    aes(
      x = dbRDA1,
      y = dbRDA2,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 3,
    colour = "black"
  ) +
  scale_color_manual(
    values = season_colors,
    breaks = orderSeason,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    breaks = orderSites,
    drop = FALSE
  ) +
  labs(
    x = ord_labels(dbRDA_18S)[1],
    y = ord_labels(dbRDA_18S)[2],
    color = "Season",
    shape = "Site"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-lim_dbRDA_18S, lim_dbRDA_18S),
    ylim = c(-lim_dbRDA_18S, lim_dbRDA_18S),
    clip = "off"
  ) +
  theme_figures()

## Plot dbRDA 16S
# Extract environmental vector scores
vectors_16S <- as.data.frame(
  scores(
    envfit_16S,
    display = "vectors",
    choices = 1:2
  )
) %>%
  rownames_to_column("Variable")

# Add envfit p values
vectors_16S$p_value <- envfit_16S$vectors$pvals

# Retain significant vectors
vectors_16S <- vectors_16S %>%
  filter(p_value < 0.05)

sites_16S <- as.data.frame(
  scores(
    dbRDA_16S,
    display = "sites",
    choices = 1:2
  )
) %>%
  rownames_to_column("Sample") %>%
  left_join(
    env16S %>%
      rownames_to_column("Sample"),
    by = "Sample"
  )

# Scale arrows for plotting
if (nrow(vectors_16S) > 0) {
  arrow_mult_16S <- 0.8 * min(
    max(abs(sites_16S$dbRDA1), na.rm = TRUE) /
      max(abs(vectors_16S$dbRDA1), na.rm = TRUE),
    max(abs(sites_16S$dbRDA2), na.rm = TRUE) /
      max(abs(vectors_16S$dbRDA2), na.rm = TRUE)
  )
  vectors_16S <- vectors_16S %>%
    mutate(
      dbRDA1 = dbRDA1 * arrow_mult_16S,
      dbRDA2 = dbRDA2 * arrow_mult_16S,
      label_x = dbRDA1 * 1.10,
      label_y = dbRDA2 * 1.10
    )
}

lim_dbRDA_16S <- 1.10 * max(
  abs(c(
    sites_16S$dbRDA1,
    sites_16S$dbRDA2,
    vectors_16S$dbRDA1,
    vectors_16S$dbRDA2
  )),
  na.rm = TRUE
)

dbRDA_plot_16S <- ggplot(
  sites_16S,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Season,
    shape = Site
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    stroke = 1.3
  ) +
  geom_segment(
    data = vectors_16S,
    aes(
      x = 0,
      y = 0,
      xend = dbRDA1,
      yend = dbRDA2
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "darkgray",
    arrow = arrow(
      length = unit(0.18, "cm")
    )
  ) +
  geom_text_repel(
    data = vectors_16S,
    aes(
      x = dbRDA1,
      y = dbRDA2,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 3,
    colour = "black"
  ) +
  scale_color_manual(
    values = season_colors,
    breaks = orderSeason,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    breaks = orderSites,
    drop = FALSE
  ) +
  labs(
    x = ord_labels(dbRDA_16S)[1],
    y = ord_labels(dbRDA_16S)[2],
    color = "Season",
    shape = "Site"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-lim_dbRDA_16S, lim_dbRDA_16S),
    ylim = c(-lim_dbRDA_16S, lim_dbRDA_16S),
    clip = "off"
  ) +
  theme_figures()


### Fig 5. Variance Partitioning (and RDAs)

### Variance Partitioning 18S
# Extract OTU table
taxa_18S <- as(otu_table(ps_18S), "matrix")
if (taxa_are_rows(ps_18S)) {
  taxa_18S <- t(taxa_18S)
}
# Hellinger transformation
taxa_18S_hel <- decostand(
  taxa_18S,
  method = "hellinger"
)

# Extract metadata
metadata_18S <- data.frame(sample_data(ps_18S))

## Predictor matrices
# Temporal component Season
season_18S <- model.matrix(
  ~ Season,
  data = metadata_18S
)[, -1, drop = FALSE]


# Spatial component Habitat
habitat_18S <- model.matrix(
  ~ Habitat,
  data = metadata_18S
)[, -1, drop = FALSE]

# Environmental component (only parameters retained after VIF)
env_18S <- metadata_18S[,c(
  "Temp",
  "Salinity",
  "NO3_NO2",
  "PO4",
  "NH4",
  "ChlTot",
  "DO",
  "SiO2",
  "pH",
  "Turbidity"
)]

# Standardize environmental variables
env_18S_std <- as.data.frame(
  decostand(
    env_18S,
    method = "standardize"
  )
)

## Variance partitioning
varpart_18S <- varpart(
  taxa_18S_hel,
  season_18S,
  habitat_18S,
  env_18S_std
)

varpart_18S

plot(
  varpart_18S,
  Xnames = c(
    "Season",
    "Habitat",
    "Environment"
  )
)

## Partial RDA (test pure fractions)
# Pure seasonal fraction
rda_18S_season <- rda(
  taxa_18S_hel,
  season_18S,
  cbind(
    habitat_18S,
    env_18S_std
  )
)

set.seed(321)
test_18S_season <- anova(
  rda_18S_season, 
  permutations = 9999
)
test_18S_season

# Pure habitat fraction
rda_18S_habitat <- rda(
  taxa_18S_hel, 
  habitat_18S, 
  cbind(
    season_18S, 
    env_18S_std
  )
)

set.seed(321)
test_18S_habitat <- anova(
  rda_18S_habitat,
  permutations = 9999
)
test_18S_habitat

# Pure environmental fraction
rda_18S_env <- rda(
  taxa_18S_hel,
  env_18S_std,
  cbind(
    season_18S,
    habitat_18S
  )
)

set.seed(321)
test_18S_env <- anova(
  rda_18S_env,
  permutations = 9999
)
test_18S_env

### Variance Partitioning 16S
# Extract OTU table
taxa_16S <- as(otu_table(ps_16S), "matrix")
if (taxa_are_rows(ps_16S)) {
  taxa_16S <- t(taxa_16S)
}
# Hellinger transformation
taxa_16S_hel <- decostand(
  taxa_16S,
  method = "hellinger"
)

# Extract metadata
metadata_16S <- data.frame(sample_data(ps_16S))

## Predictor matrices
# Temporal component Season
season_16S <- model.matrix(
  ~ Season,
  data = metadata_16S
)[, -1, drop = FALSE]


# Spatial component Habitat
habitat_16S <- model.matrix(
  ~ Habitat,
  data = metadata_16S
)[, -1, drop = FALSE]

# Environmental component (only parameters retained after VIF)
env_16S <- metadata_16S[,c(
  "Temp",
  "Salinity",
  "NO3_NO2",
  "PO4",
  "NH4",
  "ChlTot",
  "DO",
  "SiO2",
  "pH",
  "Turbidity"
)]

# Standardize environmental variables
env_16S_std <- as.data.frame(
  decostand(
    env_16S,
    method = "standardize"
  )
)

## Variance partitioning
varpart_16S <- varpart(
  taxa_16S_hel,
  season_16S,
  habitat_16S,
  env_16S_std
)

varpart_16S

plot(
  varpart_16S,
  Xnames = c(
    "Season",
    "Habitat",
    "Environment"
  )
)

## Partial RDA (test pure fractions)
# Pure seasonal fraction
rda_16S_season <- rda(
  taxa_16S_hel,
  season_16S,
  cbind(
    habitat_16S,
    env_16S_std
  )
)

set.seed(321)
test_16S_season <- anova(
  rda_16S_season, 
  permutations = 9999
)
test_16S_season

# Pure habitat fraction
rda_16S_habitat <- rda(
  taxa_16S_hel, 
  habitat_16S, 
  cbind(
    season_16S, 
    env_16S_std
  )
)

set.seed(321)
test_16S_habitat <- anova(
  rda_16S_habitat,
  permutations = 9999
)
test_16S_habitat

# Pure environmental fraction
rda_16S_env <- rda(
  taxa_16S_hel,
  env_16S_std,
  cbind(
    season_16S,
    habitat_16S
  )
)

set.seed(321)
test_16S_env <- anova(
  rda_16S_env,
  permutations = 9999
)
test_16S_env


### Procrustes analysis and Mantel tests

## Load averaged Bray-Curtis distance matrices
bray_results <- readRDS(
  "average_Bray_distances_18S_16S.rds"
)
dist_Bray18S <- bray_results$dist_Bray18S
dist_Bray16S <- bray_results$dist_Bray16S

# Get samples present in both datasets
common_samples <- labels(dist_Bray18S)[
  labels(dist_Bray18S) %in% labels(dist_Bray16S)
]

dist_Bray18S_proc <- as.dist(
  as.matrix(dist_Bray18S)[
    common_samples,
    common_samples,
    drop = FALSE
  ]
)

dist_Bray16S_proc <- as.dist(
  as.matrix(dist_Bray16S)[
    common_samples,
    common_samples,
    drop = FALSE
  ]
)

stopifnot(
  identical(
    labels(dist_Bray18S_proc),
    labels(dist_Bray16S_proc)
  )
)

ps_18S_proc <- prune_samples(
  sample_names(ps_18S) %in% common_samples,
  ps_18S
)

ps_16S_proc <- prune_samples(
  sample_names(ps_16S) %in% common_samples,
  ps_16S
)

# Generate ordinations
set.seed(2001)
NMDS_proc_18S <- ordinate(
  ps_18S_proc,
  method = "NMDS",
  distance = dist_Bray18S_proc,
  autotransform = FALSE,
  trymax = 500
)

set.seed(2001)
NMDS_proc_16S <- ordinate(
  ps_16S_proc,
  method = "NMDS",
  distance = dist_Bray16S_proc,
  autotransform = FALSE,
  trymax = 500
)

NMDS_proc_18S$stress
NMDS_proc_16S$stress

scores_proc_18S <- scores(
  NMDS_proc_18S,
  display = "sites",
  choices = 1:2
)

scores_proc_16S <- scores(
  NMDS_proc_16S,
  display = "sites",
  choices = 1:2
)

scores_proc_18S <- scores_proc_18S[
  common_samples,
  ,
  drop = FALSE
]

scores_proc_16S <- scores_proc_16S[
  common_samples,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(scores_proc_18S),
    rownames(scores_proc_16S)
  )
)

## Procrustes
proc_18S_16S <- procrustes(
  scores_proc_18S,
  scores_proc_16S,
  symmetric = TRUE
)

set.seed(2001)
protest_18S_16S <- protest(
  scores_proc_18S,
  scores_proc_16S,
  permutations = 9999
)

proc_18S_16S
protest_18S_16S

## Mantel test
set.seed(2001)
mantel_18S_16S <- mantel(
  dist_Bray18S_proc,
  dist_Bray16S_proc,
  method = "spearman",
  permutations = 9999
)

mantel_r <- unname(mantel_18S_16S$statistic)
mantel_p <- mantel_18S_16S$signif

mantel_results <- data.frame(
  Samples = length(labels(dist_Bray18S_proc)),
  Mantel_r = mantel_r,
  P_value = mantel_p
)

mantel_results


### ISA

## For 18S
counts_18S <- data.frame(otu_table(ps_18Sr))
counts_18S <- t(counts_18S)
env_16S <- data.frame(sample_data(ps_18Sr))
sites_16S <- env$Site

set.seed(123)
ISA_18S <- multipatt(counts_18S, sites_18S, func = "IndVal.g", control = how(nperm=9999))
summary(ISA_18S)

## For 16S
counts_16S <- data.frame(otu_table(ps_16Sr))
counts_16S <- t(counts_16S)
env_16S <- data.frame(sample_data(ps_16Sr))
sites_16S <- env$Site

set.seed(123)
ISA_16S <- multipatt(counts_16S, sites_16S, func = "IndVal.g", control = how(nperm=9999))
summary(ISA_16S)
