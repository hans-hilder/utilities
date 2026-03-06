library(tidyverse)
library(ggplot2)
library(myplotfunction)

setwd("C:/Users/aarfer/OneDrive - NOC/Documents/BIOCARBON+/COMPARISON")

#--- load raw data---
detailed_data <- readr::read_tsv(
  "C:/Users/aarfer/OneDrive - NOC/Documents/BIOCARBON+/COMPARISON/UVP6_DETAILED_TSVS/export_detailed_20260109_11_43_PAR_cal2rcf19.tsv",
  locale = locale(encoding = "windows-1252")
)
# -------------------------------
# 0. Choose binning structure
# -------------------------------
binning_scheme <- "uvp6"   # options: "global_log", "uvp6"


# -------------------------------
# 1. Define bin edges for both schemes
# -------------------------------

# --- Global log bins ---
log10_bin_width <- 0.065
global_min_um <- 1
global_max_um <- 26000
global_bin_edges_um <- 10^seq(
  floor(log10(global_min_um) / log10_bin_width) * log10_bin_width,
  ceiling(log10(global_max_um) / log10_bin_width) * log10_bin_width,
  by = log10_bin_width
)

# --- UVP6 bins ---
uvp_bins <- detailed_data %>%
  # select only LPM columns and exclude biovolume
  select(starts_with("LPM"), -contains("biovolume")) %>%
  pivot_longer(
    everything(),
    names_to = "bin_label",
    values_to = "dummy"
  ) %>%
  distinct(bin_label) %>%
  mutate(
    # extract min and max from label
    bin_min = as.numeric(str_extract(bin_label, "(?<=\\()[0-9.]+")),
    bin_max = as.numeric(str_extract(bin_label, "(?<=-)[0-9.]+(?=\\s)")),
    is_mm   = str_detect(bin_label, "mm\\)")
  ) %>%
  mutate(
    bin_min = ifelse(is_mm, bin_min * 1000, bin_min),
    bin_max = ifelse(is_mm, bin_max * 1000, bin_max)
  ) %>%
  arrange(bin_min) %>%
  filter(!is.na(bin_min))   # remove any bins that failed to parse

# --- UVP6 bin edges, adding open-ended last bin ---
uvp6_bin_edges_um <- c(uvp_bins$bin_min, 26000)

uvp6_bin_edges_um

# -------------------------------# -------------------------------# -------------------------------
# 2. Select which binning scheme to use
# -------------------------------
if (binning_scheme == "global_log") {
  bin_edges_um <- global_bin_edges_um
} else if (binning_scheme == "uvp6") {
  bin_edges_um <- uvp6_bin_edges_um
} else {
  stop("Unknown binning_scheme")
}

bin_mids_um   <- sqrt(bin_edges_um[-length(bin_edges_um)] * bin_edges_um[-1])
bin_widths_um <- diff(bin_edges_um)





#######################################################
#constants
pxl <- 73 # this is the pixel size in the export summary


# -------------------------------
# 1. Load and prepare PAR particles only
# -------------------------------
par <- readr::read_tsv(
  "C:/Users/aarfer/OneDrive - NOC/Documents/BIOCARBON+/COMPARISON/UVP6_RAW_TSVS/20240621230846_CAL2rcf19_PAR_raw_20260109_10_47.tsv",
  locale = locale(encoding = "windows-1252")
)

Aa  <- 0.0023    # check your instrument's UVPdb/calibration record
Exp <- 1.136   # check your instrument's UVPdb/calibration record

particles_par <- par %>%
  mutate(
    # Step 1: calibrated surface area in mm²
    Sm_mm2 = Aa * (area ^ Exp),
    # Step 2: ESD in mm from circular area assumption
    esd_mm = 2 * sqrt(Sm_mm2 / pi),
    # Step 3: convert to µm
    esd_um = esd_mm * 1000,
    nbr = nbr,
    source = "PAR"
  ) %>%
  select(depth, esd_um, nbr, source)

cat("=== Particle counts ===\n")
cat(sprintf("PAR particles: %d\n", nrow(particles_par)))

# Check size range
cat("\n=== Size distribution ===\n")
particles_par %>%
  summarise(
    count = n(),
    min_esd = min(esd_um),
    median_esd = median(esd_um),
    max_esd = max(esd_um)
  ) %>%
  print()


sum(particles_par$nbr)

# -------------------------------
# 2. Depth bins (from uvp6 detailed export)
# -------------------------------
sample_vol <- detailed_data

depth_bins <- sample_vol %>%
  select(Profile, depth_mid = `Depth [m]`, sampled_volume_L = `Sampled volume [L]`) %>%
  arrange(depth_mid)



particles_binned <- particles_par %>%
  rowwise() %>%
  mutate(
    # find the index of the closest depth_mid
    depth_bin_index = which.min(abs(depth - depth_bins$depth_mid))
  ) %>%
  ungroup() %>%
  # join to depth_bins via row number
  mutate(row_num = depth_bin_index) %>%
  left_join(
    depth_bins %>% mutate(row_num = row_number()),
    by = "row_num"
  ) %>%
  select(-row_num, -depth_bin_index)




stopifnot(all(c("Profile","depth_mid","sampled_volume_L") %in% colnames(particles_binned)))



# -------------------------------
# 4. Depth-resolved spectra (Lisst/Holo style)
# -------------------------------
uvp_spectra_depth_par <- particles_binned %>%
  filter(!is.na(depth_mid), source == "PAR") %>%
  group_by(Profile, depth_mid) %>%
  group_modify(~{
    
    vol_L <- unique(.x$sampled_volume_L)
    if (length(vol_L) != 1 || vol_L <= 0) vol_L <- NA_real_
    
    # Expand by nbr
    esd_vector <- rep(.x$esd_um, .x$nbr)
    
    h <- hist(
      esd_vector,
      breaks = bin_edges_um,
      plot = FALSE,
      right = TRUE,
      include.lowest = TRUE
    )
    
    bin_min <- bin_edges_um[-length(bin_edges_um)]
    bin_max <- bin_edges_um[-1]
    bin_mid <- sqrt(bin_min * bin_max)
    bin_width <- bin_max - bin_min
    
    tibble(
      bin_min    = bin_min,
      bin_max    = bin_max,
      bin_mid    = bin_mid,
      bin_width  = bin_width,
      count      = h$counts,
      count_norm = if (!is.na(vol_L)) h$counts / bin_width / vol_L else NA_real_,
      sampled_volume_L = vol_L
    )
    
  }) %>%
  ungroup()




# -------------------------------
# 5. Vertically integrated spectra (Lisst/Holo style)
# -------------------------------

# Expand all particles by nbr
esd_all <- rep(particles_binned$esd_um, particles_binned$nbr)

# Compute histogram
h_integrated <- hist(
  esd_all,
  breaks = bin_edges_um,
  plot = FALSE,
  right = TRUE,
  include.lowest = TRUE
)

# Total sampled volume
total_vol_L <- depth_bins %>%
  distinct(Profile, depth_mid, sampled_volume_L) %>%
  pull(sampled_volume_L) %>%
  sum()

# Define bin edges, mids, widths
bin_min   <- bin_edges_um[-length(bin_edges_um)]
bin_max   <- bin_edges_um[-1]
bin_mid   <- sqrt(bin_min * bin_max)
bin_width <- bin_max - bin_min

# Build tibble
uvp_spectra_integrated_par <- tibble(
  bin_min        = bin_min,
  bin_max        = bin_max,
  bin_mid        = bin_mid,
  bin_width      = bin_width,
  total_particles = h_integrated$counts,
  total_volume_L  = total_vol_L
) %>%
  mutate(
    count_norm = total_particles / bin_width / total_volume_L
  )



cat("Total particles (depth-resolved):", sum(uvp_spectra_depth_par$count), "\n")
cat("Total particles (integrated):", sum(uvp_spectra_integrated_par$total_particles), "\n")

#sanity check
stopifnot(
  sum(uvp_spectra_integrated_par$total_particles) ==
    sum(particles_par$nbr)
)

# -------------------------------
# 6. Load uvp6 official detailed PAR export
# -------------------------------
uvp_official_par <- detailed_data

# Select relevant columns: Profile, depth, sampled volume, LPM bins
uvp_counts <- uvp_official_par %>%
  select(Profile, depth_mid = `Depth [m]`, sampled_volume_L = `Sampled volume [L]`,
         starts_with("LPM")) %>%
  select(-matches("mm3"))  # remove any volume-normalized columns

# Convert to long format
uvp_counts_long <- uvp_counts %>%
  pivot_longer(
    cols = starts_with("LPM"),
    names_to = "bin_label",
    values_to = "lpm_value"
  ) %>%
  filter(!is.na(lpm_value)) %>%
  mutate(
    total_particles = lpm_value * sampled_volume_L   # particles per bin
  )

# Sum total particles over all bins and depths
total_particles <- sum(uvp_counts_long$total_particles, na.rm = TRUE)
cat("Total particles in uvp6 official PAR:", total_particles, "\n")


uvp_official_long <- uvp_official_par %>%
  select(Profile, depth_mid = `Depth [m]`, sampled_volume_L = `Sampled volume [L]`,
         starts_with("LPM")) %>%
  select(-matches("mm3")) %>%
  pivot_longer(cols = starts_with("LPM"), names_to = "bin_label", values_to = "lpm_value") %>%
  mutate(
    bin_min   = as.numeric(str_extract(bin_label, "(?<=\\()[0-9.]+")),
    bin_max   = as.numeric(str_extract(bin_label, "(?<=-)[0-9.]+(?=\\s)")),
    is_mm     = str_detect(bin_label, "mm\\)"),
    bin_min   = ifelse(is_mm, bin_min * 1000, bin_min),
    bin_max   = ifelse(is_mm, bin_max * 1000, bin_max),
    bin_mid   = sqrt(bin_min * bin_max),
    bin_width = bin_max - bin_min
  ) %>%
  filter(!is.na(bin_min), !is.na(lpm_value)) %>%   # removed lpm_value > 0
  mutate(count_norm = lpm_value / bin_width)

uvp_official_integrated <- uvp_official_long %>%
  group_by( bin_min, bin_max) %>%           # removed bin_mid and bin_width from group_by
  summarise(
    total_particles = sum(lpm_value * sampled_volume_L, na.rm = TRUE),
    total_volume_L  = sum(sampled_volume_L,             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bin_mid    = sqrt(bin_min * bin_max),
    bin_width  = bin_max - bin_min,
    count_norm = total_particles / bin_width / total_volume_L
  )


# -------------------------------
# 8. Plot: Vertically integrated comparison
# -------------------------------
p_integrated <- ggplot() +
  geom_point(data = uvp_spectra_integrated_par, aes(bin_mid, count_norm, color = "PAR raw"), size = 2.5) +
  geom_line(data = uvp_spectra_integrated_par, aes(bin_mid, count_norm, color = "PAR raw"), alpha = 0.4) +
  geom_point(data = uvp_official_integrated, aes(bin_mid, count_norm, color = "uvp6 detailed"), size = 2.5, shape = 17) +
  geom_line(data = uvp_official_integrated, aes(bin_mid, count_norm, color = "uvp6 detailed"), alpha = 0.5, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10()+
  scale_color_manual(values = c("PAR raw" = "blue", "uvp6 detailed" = "red")) +
  theme_custom() +
  labs(
    x = "ESD (µm)",
    y = expression(N / (Delta*D * V[total])),
    title = "Vertically Integrated uvp6 Size Spectra Comparison",
    color = "Source") +
  theme(legend.position = "top")

# Print plots
print(p_integrated)



