# required packages
library(ggplot2)
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)
library(RColorBrewer)
library(scales)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(detectCores() - 2) # update this if you want to use more cores

# load cost change cases
source("code/results/vary_movement_params_ebola.R")
# load plotting utilities
source("code/results/plotting_utils.R")

#### FIGURE 8: Ebola controls with changing movement params -----------------
chol_titles <- c("(b) increase movement from patch 1 to patch 2", 
                 "(c) increase movement from patch 2 to  patch 1", 
                 "(d) increase movement in both directions",
                 "(a) base case")
names(chol_titles) <- c("m1",
                        "m2", 
                        "m1_m2",
                        "base")

# reorder columns to plot base, m1, m2, m1_m1
plot_order = c("base", "m1", "m2", "m1_m2")


# plot controls over time in each patch
p<- list()
for(i in 1:(length(plot_order))){
  p[[i]] <- states %>%
    filter(id == plot_order[i], 
           #control_type == "uniform",
           variable %in% c("v", "u")) %>%
    mutate(patch_control = paste(patch, control_type, sep = "-")) %>%
    ggplot(aes(x = time, y = value, color = patch_control, linetype = patch_control)) +
    geom_line() +
    facet_grid(rows = vars(variable),
               labeller = labeller(variable = ebola_control_labs),
               switch = "y",
               scales = "free") +
    labs(color = "Patch:", 
         linetype = "Control type:",
         subtitle = chol_titles[names(chol_titles) == plot_order[i]]
    ) +
    scale_color_manual(values = ebola_cols, labels = ebola_gov_labs, name = "Control:") +
    scale_linetype_manual(values = ebola_ltys, labels = ebola_gov_labs, name = "Control:") +
    theme_minimal_grid(12) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(2,"lines"),
      legend.title = element_blank(),
      panel.spacing = unit(c(0), "lines"),
      plot.margin = unit(c(0.1,0.1,0.1,2), "lines"),
      strip.background = element_blank(),
      strip.placement = "outside" 
    ) +
    panel_border()
}

# fix legend
# make dummy plot with the correct legend
new_lty <-  c("solid", "23", "solid", "23")
names(new_lty) <- c("1-unique", "1-uniform",  "2-unique","2-uniform")
leg <- states  %>% 
  filter(variable == "I" & control_type != "none") %>%
  mutate(patch_control = paste(patch, control_type, sep = "-")) %>% 
  ggplot(aes(x = time, y = value, color = patch_control, linetype = patch_control)) +
  geom_line() +
  scale_color_manual(values = ebola_cols, labels = ebola_gov_labs, name = "Control:") +
  scale_linetype_manual(values = new_lty, labels = ebola_gov_labs, name = "Control:") +
  theme_bw() +
  theme(legend.position = "bottom")
leg <- get_legend(leg)


p <- lapply(p, function(i){i + theme(legend.position = "none")})
p[[2]] <- p[[2]] + theme(strip.text.y = element_blank())
p[[3]] <- p[[3]] + theme(strip.text.y = element_blank())
p[[4]] <- p[[4]] + theme(strip.text.y = element_blank())
#p <- lapply(p[[2:3]], function(i){i + theme(axis.title.y = element_blank())})
bl <- ggplot() + theme_nothing()

ebola_vary_mvmt_plots <- plot_grid(plot_grid(plotlist = p, ncol = 4),
                                   leg,
                                   ncol = 1, 
                                   rel_heights = c(0.95,0.05))
ebola_vary_mvmt_plots

ggsave("figures/figure8-Ebola_vary_movement.eps",
       ebola_vary_mvmt_plots,
       width = 8, height = 3, scale = 2)

