# Script to generate results graph for power and bias + se for the DEFINITIVE simulations

rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

# read excel results file ------------------------------------------------------

read_scenario_workbook <- function(path) {
  raw <- read_excel(path, sheet = 1, col_names = FALSE)
  raw <- raw[rowSums(!is.na(raw)) > 0, ]
  
  # scenario title rows: first column filled, all other columns empty
  is_scenario_row <- !is.na(raw[[1]]) &
    rowSums(!is.na(raw[, -1, drop = FALSE])) == 0
  
  # fill scenario names downward
  raw$scenario <- NA_character_
  raw$scenario[is_scenario_row] <- as.character(raw[[1]][is_scenario_row])
  raw$scenario <- tidyr::fill(raw, scenario, .direction = "down")$scenario
  
  # get col indexes of missingness percentages
  miss_idx <- which(grepl("missing|dropouts", raw[raw[ ,1] == "method", ]) & !(grepl("percent", raw[raw[ ,1] == "method", ])))
  miss_nams <- unique(raw[raw[1,1] == "method", ])[miss_idx]
  
  # data rows: exclude title rows and header rows
  dat <- raw[!is_scenario_row & as.character(raw[[1]]) != "method", ]
  
  # rename columns
  names(dat)[1:6] <- c("method", "n_sim", "bias", "se_est", "power", "ci_coverage")
  
  dat <- dat %>%
    dplyr::select(scenario, method, n_sim, bias, se_est, power, ci_coverage) %>%
    mutate(
      scenario = as.character(scenario),
      method   = as.character(method),
      n_sim    = as.numeric(n_sim),
      bias    = as.numeric(bias),
      se_est  = as.numeric(se_est),
      power   = as.numeric(power),
      ci_coverage = as.numeric(ci_coverage)
    )
  
  dat
}


file_path <- "Z:/EU_Projekt_DEFINITIVE/R_Imputation/sim_results_June25nd_pmm_n50000.xlsx"
df <- read_scenario_workbook(file_path)

# order scenarios in the order they appear in the file
scenario_levels <- unique(df$scenario)
df <- df %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    # scenarios like referencecase and referencecase_H0 get the same base color
    scenario_group = str_remove(as.character(scenario), "_H0$")
  )


df <- df %>%
  mutate(
    scenario = factor(scenario, levels = unique(scenario)),
    
    scenario_block =
      scenario %>%
      str_remove("_H0$") %>%
      str_remove("_[0-9]+percent$")
  )

df$method[grepl("spline", df$method)] <- gsub("_", "", df$method[grepl("spline", df$method)])
df$method <- gsub("_", " ", df$method)

method_cols <- c(
  "raw data"   = "red",
  "mean timeframe"  = "orange",
  "nn"              = "olivedrab2",
  "nn timeframe"    = "olivedrab",
  "complete case"        = "orchid1",
  "spline"          = "skyblue1",
  "spline2"        = "skyblue3",
  "spline3"        = "skyblue4"
)


combined_plot <- function(data, powerline, powerlims, powernam, power_se = sqrt((0.5*0.5)/n_sim_runs)){
  # order scenarios exactly as they appear
  scenario_order <- levels(data$scenario)
  
  # identify block boundaries
  scenario_info <- data %>%
    distinct(scenario, scenario_block) %>%
    mutate(x = seq_len(n()))
  
  scenario_info <- scenario_info %>%
    mutate(
      scenario_block = factor(scenario_block, levels = unique(scenario_block))
    )
  scenario_info <- scenario_info %>%
    mutate(
      tick_label = c("10%", "25%")[((row_number()-1) %% 2) + 1]
    )
  
  # block_info <- scenario_info %>%
  #   group_by(scenario_block) %>%
  #   summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, center = mean(x), .groups = "drop")
  
  block_info <- scenario_info %>%
    group_by(scenario_block) %>%
    summarise(
      center = mean(x),
      group_label = sub("_(10|25)percent.*", "", first(scenario)),
      .groups = "drop"
    )
  
  block_info$group_label <- gsub("_", " \n", block_info$group_label)
  

  
  # alternating background shading
  block_info$shade <- rep(c(TRUE, FALSE),
                          length.out = nrow(block_info))
  
  data <- merge(data, scenario_info[, c("scenario", "x")])
  
  power_low_up <- c(min(powerlims), c(max(powerlims)))
  
  pointsize <- 1.3

  # power plot
  
  p_power <- ggplot(data, aes(x = x, y = power, colour = method, group = method)) +
    
    geom_hline(aes(yintercept = powerline)) +
    
    geom_rect(data = block_info %>% filter(shade),
              aes(xmin = center-1, xmax = center+1, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.08) +
    
    geom_point(position = position_dodge(width = 0.6), size = pointsize) +
    
    scale_colour_manual(values = method_cols) +
    
    labs(y = powernam, x = NULL, colour = "Method") +
    
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm"),
          legend.spacing.x = unit(0.6, "cm"),
          panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.6),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 1.2
          )) +
    guides(
      colour = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 2.5))
    )+
    
    scale_y_continuous(breaks = powerlims, limits = power_low_up, minor_breaks = F) +
    scale_x_continuous(
      breaks = c(scenario_info$x-0.5, max(scenario_info$x)+0.5),
      minor_breaks = F,
      expand = expansion(mult = 0, add = 0),
      limits = c(0.5,12.5)
    )
    
  #p_power
  
  if (!is.na(power_se)){
    p_power <- p_power + 
      geom_hline(
        yintercept = c(powerline - power_se, powerline + power_se),
        linetype = "dashed",
        linewidth = 0.2
      )
  }
  
  
  # bias plot
  
  step <- 0.2
  
  biaslims <- seq(floor(min(data$bias)/step)*step, ceiling(max(data$bias)/step)*step , step)
  
  p_bias <- ggplot(data, aes(x = x, y = bias, colour = method, group = method)) +
    
    geom_hline(yintercept = 0) +
    
    geom_rect(data = block_info %>% filter(shade),
              aes(xmin = center-1, xmax = center+1, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.08) +
    
    geom_errorbar(aes(ymin = bias - se_est, ymax = bias + se_est),
                  position = position_dodge(width = 0.6), width = 0.2) +
    
    geom_point(position = position_dodge(width = 0.6), size = pointsize) +
    
    scale_colour_manual(values = method_cols) +
    
    labs(y = "Bias ± SE", x = NULL, colour = "Method") +
    
    theme_bw() +
    theme(legend.position = "none", 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_line(colour = "grey80", linewidth = 0.6),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.2
      )
    ) +
    
    scale_y_continuous(breaks = biaslims, minor_breaks = F) +
    scale_x_continuous(
      breaks = scenario_info$x,
      labels = scenario_info$tick_label,
      minor_breaks = c(scenario_info$x-0.5, max(scenario_info$x)+0.5),
      expand = expansion(mult = 0, add = 0),
      limits = c(0.5,12.5)
    ) #+
    
    # annotate(
    #   "text",
    #   x = block_info$center,
    #   y = min(biaslims)-0.35,
    #   label = block_info$group_label,
    #   fontface = "bold"
    # ) +
    # 
    # coord_cartesian(
    #   ylim = c(min(biaslims)-0.1, max(biaslims)),
    #   clip = "off"
    # ) +
    # 
    # theme(
    #   plot.margin = margin(5.5, 5.5, 30, 5.5)
    # )
  
    #p_bias
  
  coverage_se <- sqrt((0.95*0.05)/n_sim_runs)
  
  p_coverage <- ggplot(data, aes(x = x, y = ci_coverage, colour = method, group = method)) +
    
    geom_hline(yintercept = 0.95) +
    geom_hline(
      yintercept = c(0.95 - coverage_se, 0.95 + coverage_se),
      linetype = "dashed",
      linewidth = 0.2
    ) +
    
    geom_rect(data = block_info %>% filter(shade),
              aes(xmin = center-1, xmax = center+1, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.08) +
    
    geom_point(position = position_dodge(width = 0.6), size = pointsize) +
    
    scale_colour_manual(values = method_cols) +
    
    labs(y = "CI Coverage", x = NULL, colour = "Method") +
    
    theme_bw() +
    theme(legend.position = "none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(colour = "grey80", linewidth = 0.6)
    ) +
    
    scale_y_continuous(breaks = seq(0.9,1,0.025), minor_breaks = F) +
    scale_x_continuous(
      breaks = scenario_info$x,
      labels = scenario_info$tick_label,
      minor_breaks = c(scenario_info$x-0.5, max(scenario_info$x)+0.5),
      expand = expansion(mult = 0, add = 0),
      limits = c(0.5,12.5)
    ) +
    
    # annotate(
    #   "text",
    #   x = block_info$center,
    #   y = 0.9,
    #   label = block_info$group_label,
    #   fontface = "bold"
    # ) +
    # 
    # coord_cartesian(
    #   ylim = c(0.9, 1),
    #   clip = "off"
    # ) +
    
    geom_text(
      data = block_info,
      aes(x = center, y = -Inf, label = group_label),
      inherit.aes = FALSE,
      vjust = 2,
      fontface = "bold"
    ) +
    
    
    coord_cartesian(
      ylim = c(0.9, 1),
      clip = "off"
    ) +

    
    theme(
      plot.margin = margin(5.5, 5.5, 35, 5.5),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.01
      )
    )
  

  #p_coverage
  
  # combine plots
  
  combined_plot <- (p_power / p_bias / p_coverage) +  plot_layout(heights = c(1, 1, 1))#, guides = "collect") #& theme(legend.position = "bottom")
  
  combined_plot
  
}

# method order
method_order <- c("raw data", "nn", "nn timeframe", "mean timeframe", "spline", "spline2", "spline3", "complete case")

n_sim_runs <- 50000

## H1 graphs

df_H1 <- df[!(grepl("_H0", df$scenario)),]
df_H1$method <- factor(df_H1$method, levels = method_order)
plot_H1 <- combined_plot(df_H1, powerlims = seq(0.6, 1, 0.05), powerline = 0.8, powernam = "Power", power_se = NA)
ggsave("Z:/EU_Projekt_DEFINITIVE/R_Imputation/H1_July6th_pmm_n50000_wCI.pdf", plot_H1, height = 6.5, width = 7)


## H0 graphs
# Rejection rate statt Power
df_H0 <- df[(grepl("_H0", df$scenario)),]
df_H0$method <- factor(df_H0$method, levels = method_order)
plot_H0 <- combined_plot(df_H0, powerlims = seq(0.01, 0.04, 0.005), powerline = 0.025, powernam = "Rejection rate", power_se = sqrt((0.025*(1-0.025))/n_sim_runs))
ggsave("Z:/EU_Projekt_DEFINITIVE/R_Imputation/H0_July6th_pmm_n50000_wCI.pdf", plot_H0, height = 6.5, width = 7)









## plot demonstrating the different visit regimes
regime1 <- data.frame(time = seq(0, 5*3*7, 3*7), y = 2, label = "Regime 1")
regime2 <- data.frame(time = seq(0, 9*2*7, 2*7), y = 1, label = "Regime 2")

regime_df <- rbind(regime1, regime2)
regime_df$shape <- "Neoadjuvant visit"
regime_df$shape[regime_df$time %in% c(105,126)] <- "Pre-surgery visit"
regime_df$shape[regime_df$time %in% c(0)] <- "BL visit"

p <- ggplot(data = regime_df, aes(x = time, y = y, color = label, group = label, shape = shape)) +
  geom_vline(xintercept = seq(3*7, 12*7, 3*7), colour = "grey85", linewidth = 2) +# this is the grid
  
  geom_point(aes(size = shape)) +
  geom_line() +
  
  theme_bw() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("lightpink2", "darkolivegreen3")) +
  scale_shape_manual(breaks = c("BL visit", "Neoadjuvant visit", "Pre-surgery visit"), values = c(15, 16, 18)) +
  scale_size_manual(breaks = c("BL visit", "Neoadjuvant visit", "Pre-surgery visit"), values = c(5, 3, 7)) +
  scale_y_continuous(breaks = c(1,2), labels = c("Regime 2", "Regime 1"), minor_breaks = F, limits = c(0.8, 2.2)) +
  scale_x_continuous(breaks = seq(0,20*7, 7), labels = seq(0,20), minor_breaks = F) + 
  labs(color = "", shape = "", y = "", x = "weeks") +
  guides(
    color = "none",
    size = "none",
    shape = guide_legend(
      override.aes = list(size = c(5,3,7)  # make legend symbols bigger
    )))
p

ggsave("Z:/EU_Projekt_DEFINITIVE/R_Imputation/schematic_visit_regimes.pdf", p, height = 2.5, width = 6.5)

















