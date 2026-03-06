
################################################################################
# 1. LIBRARIES AND WORKING DIRECTORY ------------------------------------------
library(tidyverse)
library(ggplot2)
library(myplotfunction)
library(broom)
library(latex2exp)
library(purrr)
library(emmeans)
library(patchwork)
library(gridExtra)
library(grid)
library(cowplot)
library(grid)
library(readxl)

setwd("C:/Users/aarfer/OneDrive - NOC/Documents/BIOCARBON+/COMPARISON")


################################################################################
# 2. Global size bins ------------------------------------------
log10_bin_width <- 0.065

global_min_um <- 1
global_max_um <- 48000

global_bin_edges_um <- 10^seq(
  floor(log10(global_min_um) / log10_bin_width) * log10_bin_width,
  ceiling(log10(global_max_um) / log10_bin_width) * log10_bin_width,
  by = log10_bin_width
)

global_bin_mids_um   <- sqrt(global_bin_edges_um[-length(global_bin_edges_um)] *
                               global_bin_edges_um[-1])
global_bin_widths_um <- diff(global_bin_edges_um)


cat("GLOBAL BINS\n")
cat("Range:", min(global_bin_edges_um), "-", max(global_bin_edges_um), "µm\n")
cat("Number of bins:", length(global_bin_widths_um), "\n\n")


################################################################################
# 3. Constants ------------------------------------------
#uvp5 constants
uvp5_pxl <- 94  # µm per pixel (uvp5 export summary value)

#uvp6 constants
uvp6_pxl <- 73

# LISST-Holo constants
lisst_pxl <- 4.4              # µm per pixel
lisst_pxl_cm <- 4.4e-4        # cm per pixel
lisst_v_cm3 <- 1.5            # sample volume in cm³
lisst_v_l <- lisst_v_cm3 / 1000  # sample volume in L


# WeeHolo constants
wh_pxl <- 3.45           # µm per pixel
wh_pcl_cm <- 3.45e-4     # cm per pixel
wh_v_cm3 <- 12           # sample volume in cm³
wh_v_l <- wh_v_cm3 / 1000  # sample volume in L




################################################################################
# 4. Function to get uvp5_raw_PSD ------------------------------------------
# and also pull sample_volumes from the detailed data
process_uvp_raw <- function(raw_file, detailed_file, Aa, Exp, instrument_name) {
  
  # --- Extract raw file "middle section" for profile_id ---
  raw_basename <- tools::file_path_sans_ext(basename(raw_file))
  profile_id_clean <- stringr::str_extract(raw_basename, "(?<=\\d{14}_)[^_]+")
  profile_id <- paste0(profile_id_clean, "_downcast")
  if (is.na(profile_id_clean)) stop("Could not extract profile_id from: ", raw_file)
  
  # --- Load raw particles ---
  par <- read_tsv(raw_file, show_col_types = FALSE)
  particles_par <- par %>%
    mutate(
      Sm_mm2 = Aa * (area ^ Exp),       # calibrated surface area in mm²
      esd_mm = 2 * sqrt(Sm_mm2 / pi),   # ESD in mm
      esd_um = esd_mm * 1000,            # ESD in µm
      nbr    = nbr,
      source = "PAR"
    ) %>%
    select(depth, esd_um, nbr, source)
  
  # --- Load detailed export for volume info ---
  sample_vol <- read_tsv(detailed_file, locale = locale(encoding = "Windows-1252"))
  depth_bins <- sample_vol %>%
    select(Profile, depth_mid = `Depth [m]`, sampled_volume_L = `Sampled volume [L]`) %>%
    arrange(depth_mid)
  
  # --- Assign particles to depth bins ---
  particles_binned <- particles_par %>%
    rowwise() %>%
    mutate(depth_bin_index = which.min(abs(depth - depth_bins$depth_mid))) %>%
    ungroup() %>%
    mutate(row_num = depth_bin_index) %>%
    left_join(depth_bins %>% mutate(row_num = row_number()), by = "row_num") %>%
    select(-row_num, -depth_bin_index)
  
  # --- Depth-resolved spectra ---
  uvp5_spectra_depth <- particles_binned %>%
    filter(!is.na(depth_mid)) %>%
    group_by(Profile, depth_mid) %>%
    group_modify(~{
      vol_L <- unique(.x$sampled_volume_L)
      if(length(vol_L) != 1 || vol_L <= 0) vol_L <- NA_real_
      esd_vector <- rep(.x$esd_um, .x$nbr)
      h <- hist(esd_vector, breaks = global_bin_edges_um, plot = FALSE, right = TRUE, include.lowest = TRUE)
      tibble(
        bin_min    = global_bin_edges_um[-length(global_bin_edges_um)],
        bin_max    = global_bin_edges_um[-1],
        bin_mid    = global_bin_mids_um,
        bin_width  = global_bin_widths_um,
        count      = h$counts,
        count_norm = if(!is.na(vol_L)) h$counts / global_bin_widths_um / vol_L else NA_real_,
        sampled_volume_L = vol_L,
        profile_id = profile_id
      )
    }) %>% ungroup()
  
  # --- Vertically integrated spectra ---
  esd_all <- rep(particles_binned$esd_um, particles_binned$nbr)
  h_integrated <- hist(esd_all, breaks = global_bin_edges_um, plot = FALSE, right = TRUE, include.lowest = TRUE)
  total_vol_L <- depth_bins %>%
    distinct(Profile, depth_mid, sampled_volume_L) %>%
    pull(sampled_volume_L) %>%
    sum()
  
  uvp5_spectra_integrated <- tibble(
    instrument      = instrument_name,
    profile_id      = profile_id,
    bin_min         = global_bin_edges_um[-length(global_bin_edges_um)],
    bin_max         = global_bin_edges_um[-1],
    bin_mid         = global_bin_mids_um,
    bin_width       = global_bin_widths_um,
    total_particles = h_integrated$counts,
    total_volume_L  = total_vol_L
  ) %>%
    mutate(count_norm = total_particles / bin_width / total_volume_L)
  
  return(list(
    depth_resolved = uvp5_spectra_depth,
    integrated     = uvp5_spectra_integrated
  ))
}



################################################################################
# 5. Run function on all profiles ------------------------------------------
# UVP5
uvp5_raw_files    <- sort(list.files("UVP5_RAW_TSVS",    pattern = "\\.tsv$", full.names = TRUE))
uvp5_detailed_files <- list.files("UVP5_DETAILED_TSVS", pattern = "\\.tsv$", full.names = TRUE)

  all_uvp5_spectra <- map2(uvp5_raw_files, uvp5_detailed_files,
                         ~process_uvp_raw(.x, .y, Aa = 0.0051, Exp = 1.0645,
                                          instrument_name = "uvp5"))

# UVP6
uvp6_raw_files    <- sort(list.files("UVP6_RAW_TSVS",    pattern = "\\.tsv$", full.names = TRUE))
uvp6_detailed_files <- list.files("UVP6_DETAILED_TSVS", pattern = "\\.tsv$", full.names = TRUE)
  
  all_uvp6_spectra <- map2(uvp6_raw_files, uvp6_detailed_files,
                         ~process_uvp_raw(.x, .y, Aa = 0.0023, Exp = 1.136,
                                          instrument_name = "uvp6"))


# Combine integrated spectra for plotting
uvp5_integrated_all <- map_dfr(all_uvp5_spectra, "integrated")

#adjust capitilisation to amtch holo isntruments
uvp5_integrated_all$profile_id <-
  gsub("CAL2rcf19", "cal2rcf19", uvp5_integrated_all$profile_id)

# Combine integrated spectra for plotting
uvp6_integrated_all <- map_dfr(all_uvp6_spectra, "integrated")

#adjust capitilisation to amtch holo isntruments
uvp6_integrated_all$profile_id <-
  gsub("CAL2rcf19", "cal2rcf19", uvp6_integrated_all$profile_id)


################################################################################
# 6. Check particle counts ------------------------------------------
# Get the total nbr per profile from the raw TSVs
total_nbr <- map_dfr(uvp5_raw_files, ~{
  df <- read_tsv(.x, show_col_types = FALSE)
  profile_id_clean <- str_extract(tools::file_path_sans_ext(basename(.x)), "(?<=\\d{14}_)[^_]+")
  tibble(
    profile_id = tolower(paste0(profile_id_clean, "_downcast")),  # add tolower()
    expected_particles = sum(df$nbr, na.rm = TRUE)
  )
})

# Compare to integrated spectra counts
check <- uvp5_integrated_all %>%
  group_by(profile_id) %>%
  summarise(sum_particles = sum(total_particles, na.rm = TRUE), .groups = "drop") %>%
  left_join(total_nbr, by = "profile_id") %>%
  mutate(diff = sum_particles - expected_particles)

print(check)




################################################################################
# 7. LOAD LISST-HOLO DATA ------------------------------------------------------
load_lisst_data <- function(folder_path) {
  cat("LOADING LISST DATA\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
  tsv_files <- tsv_files[!grepl("_black|_summary", tsv_files)]
  
  cat("Found", length(tsv_files), "TSV file(s)\n\n")
  
  all_data <- map_df(tsv_files, function(file_path) {
    df <- read_tsv(file_path, show_col_types = FALSE, col_types = cols(.default = "c"))
    first_val <- tolower(trimws(as.character(df[[1]][1])))
    if (first_val %in% c("[t]", "[f]")) df <- df[-1, ]
    df %>% mutate(tsv_file = basename(file_path))
  })
  
  # Convert pixels to µm
  all_data <- all_data %>%
    type_convert(col_types = cols()) %>%
    mutate(
      object_major_um = as.numeric(object_major) * lisst_pxl,
      object_ESD_um = as.numeric(object_equivalent_diameter) * lisst_pxl,
      depth = as.numeric(object_depth_min),
      acq_num_images_profile = as.numeric(acq_num_images_profile)
    ) %>%
    filter(!is.na(depth))
  
  # Extract sample_id and cast_id
  all_data <- all_data %>%
    mutate(
      sample_id = as.character(sample_id),
      cast_id = str_split(sample_id, "_", simplify = TRUE)[, 1],
      direction = str_extract(acq_id, "(upcast|downcast)"),
      sample_id = paste(sample_id, direction, sep = "_")
    )
  
  # -----------------------------
  # RAW PARTICLE COUNTS
  # -----------------------------
  cat("\nRAW PARTICLE COUNTS (before filtering)\n")
  cat("--------------------------------------------------\n")
  raw_total <- nrow(all_data)
  raw_by_cast <- all_data %>%
    group_by(cast_id) %>%
    summarise(n_raw = n(), .groups = "drop")
  cat("Total raw particles:", raw_total, "\n\n")
  print(raw_by_cast)
  
  # -----------------------------
  # REMOVE DIPPING PARTICLES
  # -----------------------------
  if ("object_flag_dipping" %in% names(all_data)) {
    dipping_removed <- all_data %>%
      filter(object_flag_dipping == "dipping") %>%
      group_by(cast_id) %>%
      summarise(n_dipping_removed = n(), .groups = "drop")
    dipping_total <- sum(dipping_removed$n_dipping_removed)
    all_data <- all_data %>% filter(object_flag_dipping != "dipping")
    
    cat("\nDIPPING REMOVAL\n")
    cat("--------------------------------------------------\n")
    cat("Total removed (dipping):", dipping_total, "\n")
    cat("Percent of raw:", round(100 * dipping_total / raw_total, 2), "%\n\n")
    print(dipping_removed)
  }
  
  # -----------------------------
  # REMOVE BAD CATEGORIES
  # -----------------------------
  bad_categories <- c("duplicate", "duplicate 1", "duplicate 2", "duplicate 3",
                      "duplicate 4", "duplicate 5", "blank", "noise", "bubble")
  
  if ("object_annotation_status" %in% names(all_data) &&
      "object_annotation_category" %in% names(all_data)) {
    bad_removed <- all_data %>%
      filter(
        grepl("validated", object_annotation_status, ignore.case = TRUE) &
          object_annotation_category %in% bad_categories
      ) %>%
      group_by(cast_id) %>%
      summarise(n_bad_removed = n(), .groups = "drop")
    bad_total <- sum(bad_removed$n_bad_removed)
    all_data <- all_data %>%
      filter(
        !(grepl("validated", object_annotation_status, ignore.case = TRUE) &
            object_annotation_category %in% bad_categories))
    cat("\nBAD CATEGORY REMOVAL (validated only)\n")
    cat("--------------------------------------------------\n")
    cat("Total removed:", bad_total, "\n")
    cat("Percent of raw:", round(100 * bad_total / raw_total, 2), "%\n\n")
    print(bad_removed)
  }

  # -----------------------------
  # FULLCAST ENTRIES
  # -----------------------------
  cat("\nCreating fullcast entries for LISST...\n")
  fullcast_list <- list()
  
  for (cast in unique(all_data$cast_id)) {
    group_data <- all_data %>% filter(cast_id == cast)
    directions <- unique(group_data$direction)
    if ("upcast" %in% directions && "downcast" %in% directions) {
      cat("  Found complete cast pair for", cast, "\n")
      upcast_data <- group_data %>% filter(direction == "upcast")
      downcast_data <- group_data %>% filter(direction == "downcast")
      upcast_images <- first(upcast_data$acq_num_images_profile)
      downcast_images <- first(downcast_data$acq_num_images_profile)
      total_images <- upcast_images + downcast_images
      fullcast_df <- bind_rows(upcast_data, downcast_data) %>%
        mutate(
          sample_id = paste0(cast, "_fullcast"),
          direction = "fullcast",
          acq_num_images_profile = total_images)
      fullcast_list[[length(fullcast_list) + 1]] <- fullcast_df
      cat("    Upcast images:", upcast_images, ", Downcast images:", downcast_images,
          ", Total:", total_images, "\n")
    }
  }
  if (length(fullcast_list) > 0) {
    all_data <- bind_rows(all_data, fullcast_list)
    cat("  Created", length(fullcast_list), "fullcast entries\n")
  }
  
  # -----------------------------
  # FINAL RETAINED PARTICLES
  # -----------------------------
  cat("\nFINAL RETAINED PARTICLES\n")
  cat("--------------------------------------------------\n")
  final_total <- nrow(all_data)
  final_by_cast <- all_data %>%
    group_by(cast_id) %>%
    summarise(n_retained = n(), .groups = "drop")
  cat("Total retained:", final_total, "\n")
  cat("Percent retained:", round(100 * final_total / raw_total, 2), "%\n\n")
  print(final_by_cast)
  
  return(all_data)
}
# Reload the data
data_raw_lisst <- load_lisst_data("LISSTHOLO_TSVS")



################################################################################
# 8. LOAD WEEHOLO DATA ---------------------------------------------------------
load_weeholo_data <- function(folder_path) {
  cat("LOADING WEEHOLO DATA\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
  tsv_files <- tsv_files[!grepl("_black|_summary", tsv_files)]
  cat("Found", length(tsv_files), "TSV file(s)\n\n")
  all_data <- map_df(tsv_files, function(file_path) {
    df <- read_tsv(file_path, show_col_types = FALSE, col_types = cols(.default = "c"))
    first_val <- tolower(trimws(as.character(df[[1]][1])))
    if (first_val %in% c("[t]", "[f]")) df <- df[-1, ]
    df %>% mutate(tsv_file = basename(file_path))
  })
  # Convert pixels to µm
  all_data <- all_data %>%
    type_convert(col_types = cols()) %>%
    mutate(
      object_major_um = as.numeric(object_major) * wh_pxl,
      object_ESD_um = as.numeric(object_equivalent_diameter) * wh_pxl,
      sample_id = paste0(as.character(sample_id), "_fullcast"),
      cast_id = sample_id
    )
  
  # -----------------------------
  # RAW PARTICLE COUNTS
  # -----------------------------
  cat("\nRAW PARTICLE COUNTS (before filtering)\n")
  cat("--------------------------------------------------\n")
  raw_total <- nrow(all_data)
  raw_by_cast <- all_data %>%
    group_by(cast_id) %>%
    summarise(n_raw = n(), .groups = "drop")
  cat("Total raw particles:", raw_total, "\n\n")
  print(raw_by_cast)
  
  # -----------------------------
  # REMOVE BAD CATEGORIES
  # -----------------------------
  bad_categories <- c("duplicate", "duplicate 1", "duplicate 2", "duplicate 3",
                      "duplicate 4", "duplicate 5", "blank", "noise", "bubble")
  if ("object_annotation_status" %in% names(all_data) &&
      "object_annotation_category" %in% names(all_data)) {
    bad_removed <- all_data %>%
      filter(
        grepl("validated", object_annotation_status, ignore.case = TRUE) &
          object_annotation_category %in% bad_categories
      ) %>%
      group_by(cast_id) %>%
      summarise(n_bad_removed = n(), .groups = "drop")
    bad_total <- sum(bad_removed$n_bad_removed)
    all_data <- all_data %>%
      filter(
        !(grepl("validated", object_annotation_status, ignore.case = TRUE) &
            object_annotation_category %in% bad_categories) )
    cat("\nBAD CATEGORY REMOVAL (validated only)\n")
    cat("--------------------------------------------------\n")
    cat("Total removed:", bad_total, "\n")
    cat("Percent of raw:", round(100 * bad_total / raw_total, 2), "%\n\n")
    print(bad_removed)
  }
  
  # -----------------------------
  # FINAL RETAINED PARTICLES
  # -----------------------------
  cat("\nFINAL RETAINED PARTICLES\n")
  cat("--------------------------------------------------\n")
  final_total <- nrow(all_data)
  final_by_cast <- all_data %>%
    group_by(cast_id) %>%
    summarise(n_retained = n(), .groups = "drop")
  cat("Total retained:", final_total, "\n")
  cat("Percent retained:", round(100 * final_total / raw_total, 2), "%\n\n")
  print(final_by_cast)
  
  return(all_data)
}

data_raw_weeholo <- load_weeholo_data("WEEHOLO_TSVS")




################################################################################
# 9. Function to processes Holo images -----------------------------------------
process_holo_integrated <- function(df,
                                    id_col,
                                    size_col,
                                    volume_L_per_image,
                                    image_count_col = NULL,
                                    instrument = "LISST") {
  
  # ------------------------------------------------------------
  # 9.1 Profile ID + sample_volume
  # ------------------------------------------------------------
  profile_id <- unique(df[[id_col]])
  if (length(profile_id) != 1)
    stop("Multiple profile IDs found")
  
  n_images <- unique(df[[image_count_col]])
  n_images <- n_images[!is.na(n_images)]
  
  if (length(n_images) != 1)
    stop("Multiple image counts found for profile ", profile_id)
  
  total_volume_L <- n_images * volume_L_per_image
  
  
  # ------------------------------------------------------------
  # 9.2 Histogram particle sizes onto GLOBAL bins
  # ------------------------------------------------------------
  esd <- df[[size_col]]
  esd <- esd[is.finite(esd)]
  
  h <- hist(
    esd,
    breaks = global_bin_edges_um,
    plot  = FALSE,
    right = TRUE,
    include.lowest = TRUE
  )
  
  # ------------------------------------------------------------
  # 9.3 uvp5-compatible output
  # ------------------------------------------------------------
  tibble(
    instrument        = instrument,
    profile_id        = profile_id,
    bin_min           = global_bin_edges_um[-length(global_bin_edges_um)],
    bin_max           = global_bin_edges_um[-1],
    bin_mid           = global_bin_mids_um,
    bin_width         = global_bin_widths_um,
    total_particles   = h$counts,
    total_volume_L    = total_volume_L,
    count_norm        = h$counts / bin_width / total_volume_L
  )
}




################################################################################
# 10. Apply function------------------------------------------
# ------------------------------------------------------------
# 10.1 LISSTHOLO
# ------------------------------------------------------------
lisst_integrated <- data_raw_lisst %>%
  group_by(sample_id) %>%
  group_modify(
    ~ process_holo_integrated(
      df                 = .x,
      id_col             = "sample_id",
      size_col           = "object_ESD_um",
      volume_L_per_image = lisst_v_l,
      image_count_col    = "acq_num_images_profile",
      instrument         = "LISST-Holo"
    ),
    .keep = TRUE
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 10.2 WEEHOLO
# ------------------------------------------------------------

#manually adjusting image number = time in water minus 2 minutes which are discarded by thanga
## note dipping takes 1 minute so this excludes this
weeholo_image_counts <- tribble(
  ~cast_id,          ~acq_num_images_profile,
  "r2rcf06_fullcast",  23980,
  "s3rcf07_fullcast",  20830,
  "s3rcf09_fullcast",  20260 #this is still low look into it... we are missing 6 mins of sampling!
)


weeholo_integrated <- data_raw_weeholo %>%
  select(-acq_num_images_profile) %>%           # drop old column
  left_join(weeholo_image_counts, by = "cast_id") %>%  # join per-cast counts
  group_by(cast_id) %>%
  group_modify(
    ~ process_holo_integrated(
      df                 = .x,
      id_col             = "cast_id",
      size_col           = "object_ESD_um",
      volume_L_per_image = wh_v_l,
      image_count_col    = "acq_num_images_profile",
      instrument         = "WeeHolo"
    ),
    .keep = TRUE
  ) %>%
  ungroup()





################################################################################
# 11. Combine all dsp together ------------------------------------------
all_integrated_psd <- bind_rows(
  uvp5_integrated_all %>% mutate(instrument = "uvp5"),
  uvp6_integrated_all,
  lisst_integrated,
  weeholo_integrated
)


################################################################################
# 12. Simple plot ------------------------------------------
simple= ggplot(all_integrated_psd, aes(bin_mid, count_norm, color = instrument)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_custom() +
  facet_wrap(~profile_id)+
  labs(
    x = "ESD (µm)",
    y = expression(N / (Delta*D * V[total])),
    title = "Size Spectra (Raw, Global Bins)"
  )



################################################################################
# 13. Plots ------------------------------------------
# ------------------------------------------------------------------
# 13.1 Remove zero counts
# ------------------------------------------------------------------
psd_clean <- all_integrated_psd %>%
  filter(count_norm > 0)


# ------------------------------------------------------------------
# 13.2 Extract station + cast info
# ------------------------------------------------------------------
psd_clean <- psd_clean %>%
  extract(
    col = profile_id,
    into = c("station_id", "rcf#", "direction"),
    regex = "^([a-zA-Z]+[0-9]+)(r?cf?[0-9]+)_(downcast|upcast|fullcast)$",
    remove = FALSE
  )


# ------------------------------------------------------------------
# 13.3 Remove groups with < 10 total particles
# ------------------------------------------------------------------
psd_clean <- psd_clean %>%
  group_by(instrument, station_id, `rcf#`, direction) %>%
  filter(sum(total_particles) >= 10) %>%
  ungroup()


# ------------------------------------------------------------------
# 13.4 Keep preferred cast type (fullcast preferred over downcast)
# ------------------------------------------------------------------
preferred <- psd_clean %>%
  filter(direction %in% c("fullcast", "downcast")) %>%
  
  
  # Force cal2 rcf19 LISST-Holo to use downcast
  filter(!(instrument == "LISST-Holo" &
             station_id == "cal2" &
             `rcf#` == "rcf19" &
             direction == "fullcast")) %>%
  
  group_by(instrument, station_id, `rcf#`) %>%
  mutate(direction_rank = ifelse(direction == "fullcast", 1, 2)) %>%
  filter(direction_rank == min(direction_rank)) %>%
  ungroup() %>%
  select(-direction_rank)


# ------------------------------------------------------------------
# 13.5 Determine modal bin to start fits
# ------------------------------------------------------------------
# Compute fit limits from ALL clean data (not just preferred)
fit_limits_all <- psd_clean %>%
  group_by(profile_id, instrument) %>%
  slice_max(total_particles, n = 1, with_ties = FALSE) %>%
  select(profile_id, instrument, fit_start_bin = bin_mid)

# Apply to both preferred and psd_clean
preferred <- preferred %>%
  left_join(fit_limits_all, by = c("profile_id", "instrument")) %>%
  mutate(fit_flag = bin_mid >= fit_start_bin)

psd_clean <- psd_clean %>%
  left_join(fit_limits_all, by = c("profile_id", "instrument")) %>%
  mutate(fit_flag = bin_mid >= fit_start_bin)





psd_clean <- psd_clean %>%
  group_by(instrument, station_id, `rcf#`, direction, bin_min) %>%
  filter(sum(total_particles) >= 10) %>%
  ungroup()

preferred <- preferred %>%
  group_by(instrument, station_id, `rcf#`, direction, bin_min) %>%
  filter(sum(total_particles) >= 10) %>%
  ungroup()
# ------------------------------------------------------------------
# 13.6 Standard plot
# ------------------------------------------------------------------
# --- set colours
instrument_cols <- c(
  "LISST-Holo" = "#0072B2",
  "WeeHolo"    = "#009E73",
  "uvp5"       = "#E69F00",
  "uvp6"       = "#D81B60"
)

#compute limits of x axis scale
global_xaxis <- range(
  c(preferred$bin_mid ,
          psd_clean$bin_mid),
  na.rm = TRUE
)
global_xaxis

global_yaxis <- range(
  c(preferred$count_norm ,
    psd_clean$count_norm),
  na.rm = TRUE
)
global_yaxis
# --- HELPER FUNCTION
# df:          use `preferred` for auto cast selection, `psd_clean` to specify direction manually
# instruments: character vector of instruments to include
# directions:  NULL (use whatever is in df) or e.g. "downcast" / c("downcast", "fullcast")
plot_psd <- function(df, station, rcf, instruments, directions = NULL, title = NULL) {
  
  d <- df %>%
    filter(station_id == station, `rcf#` == rcf, instrument %in% instruments)
  
  if (!is.null(directions)) d <- d %>% filter(direction %in% directions)
  
  ggplot(d, aes(bin_mid, count_norm, color = instrument, shape = direction)) +
    geom_point(size = 2) +
    geom_smooth(data = d %>% filter(fit_flag), method = "lm", formula = y ~ x,
                se = FALSE, linewidth = 0.7, linetype = "dashed") +
    scale_x_log10( limits = c(25, 10000), breaks = c(30, 100, 300, 1000, 3000, 10000)) +
    scale_y_log10( limits = global_yaxis) +
    scale_color_manual(name = "Instrument", values = instrument_cols) +
    scale_shape_manual(name = "Cast Direction", values = c("fullcast" = 17, "downcast" = 16)) +
    labs(x = "ESD (µm)", y = expression(n~(L^{-1}~mu*m^{-1}))) +
    theme_custom() +
    theme(legend.position = "right")
}

# --- inspect available combinations
preferred %>%
  distinct(station_id, `rcf#`, instrument, direction) %>%
  arrange(station_id, `rcf#`, instrument, direction)

# --- make plots
# Use `preferred` when happy with auto cast selection
# Use `psd_clean` + directions = "downcast" when overriding
p_cal2rcf19    <- plot_psd(preferred, "cal2", "rcf19", c("LISST-Holo", "uvp5", "uvp6"))
p_r2rcf06      <- plot_psd(preferred, "r2",   "rcf06", c("LISST-Holo", "WeeHolo", "uvp5"))
p_r2rcf06_down <- plot_psd(psd_clean, "r2",   "rcf06", c("LISST-Holo", "uvp5"), directions = "downcast")
p_r4rcf14      <- plot_psd(preferred, "r4",   "rcf14", c("LISST-Holo", "uvp5"))
p_s3rcf07_full <- plot_psd(preferred, "s3",   "rcf07", c("LISST-Holo", "WeeHolo", "uvp5"))
p_s3rcf07_down <- plot_psd(psd_clean, "s3",   "rcf07", c("LISST-Holo", "uvp5"), directions = "downcast")
p_s3rcf08      <- plot_psd(psd_clean, "s3",  "rcf08", c("LISST-Holo", "uvp5"))
p_s3rcf09_full <- plot_psd(preferred, "s3",  "rcf09", c("LISST-Holo", "WeeHolo", "uvp5"))
p_s3rcf09_down <- plot_psd(psd_clean, "s3",  "rcf09", c("LISST-Holo", "uvp5"), directions = "downcast")
p_s4rcf11      <- plot_psd(preferred, "s4",  "rcf11", c("LISST-Holo", "uvp5"))
p_s4rcf12      <- plot_psd(preferred, "s4",  "rcf12", c("LISST-Holo", "uvp5"))


p_r2rcf06
p_s3rcf07_full
p_s3rcf09_full


p_cal2rcf19
p_r2rcf06_down
p_r4rcf14
p_s3rcf07_down
p_s3rcf08
p_s3rcf09_down
p_s4rcf11
p_s4rcf12

plots <- c(
  "p_r2rcf06",
  "p_s3rcf07_full",
  "p_s3rcf09_full",
  "p_cal2rcf19",
  "p_r2rcf06_down",
  "p_r4rcf14",
  "p_s3rcf07_down",
  "p_s3rcf08",
  "p_s3rcf09_down",
  "p_s4rcf11",
  "p_s4rcf12"
)

for (p in plots) {
  ggsave(
    filename = paste0("figures/corrected/", p, ".png"),
    plot = get(p),
    width = 7,
    height = 9,
    dpi = 300
  )
}

###############################################################################
# --- weeholo comparison ---#

# --- clean plots (remove legends from main panels)
p1 <- p_r2rcf06 + theme(legend.position = "none")
p2 <- p_s3rcf07_full + theme(legend.position = "none", axis.title.y = element_blank())
p3 <- p_s3rcf09_full + theme(legend.position = "none")

# --- extract all legends at once from one of the plots
get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legend_index]]
}

# Extract combined legend
legend_grob <- get_legend(p_r2rcf06)

# --- arrange top row (2 plots)
top_row <- plot_grid(
  p1, p2,
  nrow = 1,
  labels = c("a","b"),
  label_fontface = "bold",
  label_size = 14,
  rel_widths = c(1.07,1)
)

# --- arrange bottom row (plot + legend)
bottom_row <- plot_grid(
  p3, legend_grob,
  nrow = 1,
  labels = c("c",""),
  label_fontface = "bold",
  label_size = 14,
  rel_widths = c(1.07,1)
)

# --- combine top and bottom row into final 2x2 grid
weeholo_2x2 <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(3,3)  # adjust heights if needed
)

# --- display
weeholo_2x2

# --- save
#ggsave("figures/weeholo_2x2_corrected.png", plot = weeholo_2x2, width = 7, height = 6, dpi = 600)



##############################################################################
# --- LisstHolo comparison ---#
#NOTE: run original plots before this
# --- remove y-axis titles for consistency
p_s3rcf07_down <- p_s3rcf07_down      + theme(axis.title.y = element_blank())
p_s3rcf09_down <- p_s3rcf09_down + theme(axis.title.y = element_blank())
p_s4rcf12      <- p_s4rcf12      + theme(axis.title.y = element_blank())

# --- remove 'direction' from legend by overriding scale
# This works if direction is mapped in shape or color
remove_direction <- function(p) {
  p + guides(shape = "none")  # removes direction legend
}

p_r2rcf06_down <- remove_direction(p_r2rcf06_down)
p_s3rcf07_down <- remove_direction(p_s3rcf07_down)
p_s3rcf08      <- remove_direction(p_s3rcf08)
p_s3rcf09_down <- remove_direction(p_s3rcf09_down)
p_s4rcf11      <- remove_direction(p_s4rcf11)
p_s4rcf12      <- remove_direction(p_s4rcf12)

# --- combine into a 2x3 grid (2 rows, 3 columns)
LisstHolo_comparison <- 
  (p_r2rcf06_down | p_s3rcf07_down)/
  (p_s3rcf08 | p_s3rcf09_down) /
  (p_s4rcf11 | p_s4rcf12) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') &  # use capital letters if you like
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 1)
  )

LisstHolo_comparison

# --- save
#ggsave("figures/lisstholo_2x3_corrected.png", plot = LisstHolo_comparison, width = 8, height = 10, dpi = 600)



##############################################################################
# --- UVP6 comparison ---#
p_cal2rcf19 <- remove_direction(p_cal2rcf19)

p_cal2rcf19
# --- save
#ggsave("figures/uvp6_corrected.png", plot = p_cal2rcf19, width = 7, height = 5, dpi = 600)












##############################################################################
#Dataset generation
###############################################################################
# --- ensure folder exists
if (!dir.exists("outputs")) dir.create("outputs")

# --- read in log details
rcf_log <- read_excel(
  "C:/Users/aarfer/OneDrive - NOC/Documents/BIOCARBON+/COMPARISON/RCF_log_summary.xlsx"
)

# --- clean all_integrated_psd for export
all_integrated_export <- all_integrated_psd %>%
  filter(count_norm > 0) %>%
  select(-any_of(c("sample_id", "cast_id"))) %>%  # remove if present
  rename(
    bin_min_um   = bin_min,
    bin_max_um   = bin_max,
    bin_mid_um   = bin_mid,
    bin_width_um = bin_width
  ) %>%
  extract(
    col = profile_id,
    into = c("station", "cast_number", "cast_direction"),
    regex = "^([a-zA-Z0-9]+)(rcf[0-9]{1,2})_(.*)$",
    remove = FALSE
  )%>%
  mutate(
    station = tolower(station),
    cast_number = as.numeric(gsub("^rcf", "", cast_number)),
    flag = ifelse(total_particles < 10, 1, 0)
  ) %>%
  mutate(
    flag = ifelse(total_particles < 10, 1, "")
  )%>%
  left_join(
    rcf_log %>%
      select(
        Station,
        RCF,
        DateTime,
        Event,
        `Profile depth`,
        Latitude,
        Longitude
      ),
    by = c("station" = "Station",
           "cast_number" = "RCF")
  )


#data used for fits
all_integrated_export_10 <- all_integrated_export %>%
  filter(total_particles >9)

#number of size bin containing particles
n_distinct(all_integrated_export_10$bin_mid_um)

# range of sizes
all_integrated_export_10 %>%
  summarise(
    min_size_um = min(bin_min_um),
    max_size_um = max(bin_max_um)
  )
# --- save CSV
#write_csv(all_integrated_export, "outputs/all_integrated_psd_clean_corrected.csv")

