# ============================================================
# 06_figure2_heatmap.R
#
# Purpose:
# Generate a publication-quality two-panel heatmap showing
# pooled meropenem and imipenem resistance by WHO region and
# year, 2020–2023.
#
# Outputs:
# output/figures/Figure2_Heatmap.pdf
# output/figures/Figure2_Heatmap.png
# output/figures/Figure2_Heatmap.tiff
# ============================================================


# ---- Packages ------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "tidyr",
  "ggplot2",
  "viridis"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)


# ---- File paths ----------------------------------------------

input_file <-
  "data/processed/carbapenem_pooled_summary.csv"

pdf_output <-
  "output/figures/Figure2_Heatmap.pdf"

png_output <-
  "output/figures/Figure2_Heatmap.png"

tiff_output <-
  "output/figures/Figure2_Heatmap.tiff"


# ---- Check input ---------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file,
    "\nRun scripts/02_pooled_resistance_wilson_ci.R first."
  )
}

dir.create(
  "output/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Load pooled data ----------------------------------------

pooled <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---- Validate columns ----------------------------------------

required_columns <- c(
  "WHORegionName",
  "Year",
  "AntibioticName",
  "Resistance"
)

missing_columns <- setdiff(
  required_columns,
  names(pooled)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---- Ordering -------------------------------------------------

region_order <- c(
  "African Region",
  "Eastern Mediterranean Region",
  "European Region",
  "Region of the Americas",
  "South-East Asia Region",
  "Western Pacific Region"
)

antibiotic_order <- c(
  "Meropenem",
  "Imipenem"
)


# ---- Prepare heatmap data ------------------------------------

heatmap_data <- pooled %>%
  transmute(
    WHORegionName = as.character(WHORegionName),
    Year = as.integer(Year),
    AntibioticName = as.character(AntibioticName),
    Resistance = as.numeric(Resistance)
  ) %>%
  filter(
    WHORegionName %in% region_order,
    AntibioticName %in% antibiotic_order,
    Year %in% 2020:2023
  ) %>%
  mutate(
    WHORegionName = factor(
      WHORegionName,
      levels = rev(region_order)
    ),
    AntibioticName = factor(
      AntibioticName,
      levels = antibiotic_order
    ),
    Year = factor(
      Year,
      levels = 2020:2023
    )
  ) %>%
  complete(
    WHORegionName,
    AntibioticName,
    Year
  ) %>%
  mutate(
    CellLabel = if_else(
      is.na(Resistance),
      "—",
      sprintf("%.1f", Resistance)
    ),
    LabelColour = case_when(
      is.na(Resistance) ~ "black",
      Resistance >= 58 ~ "white",
      TRUE ~ "black"
    )
  )


# ---- Validate values -----------------------------------------

if (
  any(
    heatmap_data$Resistance < 0 |
    heatmap_data$Resistance > 100,
    na.rm = TRUE
  )
) {
  stop(
    "Resistance values must lie between 0 and 100."
  )
}


# ---- Create Figure 2 -----------------------------------------

figure2 <- ggplot(
  heatmap_data,
  aes(
    x = Year,
    y = WHORegionName,
    fill = Resistance
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 1.1,
    width = 0.98,
    height = 0.98
  ) +
  geom_text(
    aes(
      label = CellLabel,
      colour = LabelColour
    ),
    size = 4.1,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ AntibioticName,
    nrow = 1
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = 1,
    limits = c(25, 85),
    breaks = seq(
      30,
      80,
      by = 10
    ),
    na.value = "grey90"
  ) +
  scale_colour_identity() +
  labs(
    x = "Year",
    y = NULL,
    fill = "Resistance (%)"
  ) +
  coord_fixed(
    ratio = 0.82
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 13,
      face = "bold",
      margin = margin(
        b = 9
      )
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(
        t = 10
      )
    ),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.height = grid::unit(
      1.3,
      "cm"
    ),
    
    panel.spacing = grid::unit(
      1.4,
      "lines"
    ),
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  )


# ---- Export ---------------------------------------------------

ggsave(
  filename = pdf_output,
  plot = figure2,
  width = 11.5,
  height = 6.2,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = png_output,
  plot = figure2,
  width = 11.5,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = tiff_output,
  plot = figure2,
  width = 11.5,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ---- Confirm outputs -----------------------------------------

stopifnot(
  file.exists(pdf_output),
  file.exists(png_output),
  file.exists(tiff_output)
)

message("\nFigure 2 generated successfully.")
message("PDF: ", pdf_output)
message("PNG: ", png_output)
message("TIFF: ", tiff_output)

print(figure2)