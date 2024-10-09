Sys.setlocale("LC_ALL", "English")

library(tidyverse)
library(lubridate)
library(patchwork)
library(ggtext)


read.csv("LRV_data.csv") %>%
  filter(Group == "Phytoplankton") %>%
  ggplot(aes(x = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)[DOY], y = value)) +
  scale_x_date("Date", breaks = "1 month", date_labels = "%b %d", limits = as.Date(c("2019-02-28", "2019-11-18")), expand = c(0,0)) +
  scale_y_continuous("Concentration\n(µg/L)") +
  scale_color_manual(values = c("#4f8c1d", "#442e9e")) +
  scale_fill_manual(values = c("#4f8c1d", "#442e9e")) +
  coord_cartesian(clip = "off") +
  ggtitle("(a) Phytoplankton") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  LRV_data %>%
  filter(Group == "Daphnia spp.") %>%
  ggplot(aes(x = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)[DOY], y = value)) +
  geom_path(aes(color = name)) +
  geom_point(aes(fill = name), color = "white", size = 1.5, shape = 21, stroke = 0.5) +
  scale_x_date("Date", breaks = "1 month", date_labels = "%b %d", limits = as.Date(c("2019-02-28", "2019-11-18")), expand = c(0,0)) +
  scale_y_log10("Density\n(ind./L)") +
  scale_color_manual(values = c("#3969AC","#E73F74","#F2B701","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  scale_fill_manual(values = c("#3969AC","#E73F74","#F2B701","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  coord_cartesian(clip = "off") +
  ggtitle("<b>(b) Zooplankton:</b> <i>Daphnia</i> spp.") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        plot.title = element_markdown(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  LRV_data %>%
  filter(Group == "Cladocera spp. (other)") %>%
  ggplot(aes(x = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)[DOY], y = value)) +
  geom_path(aes(color = name)) +
  geom_point(aes(fill = name), color = "white", size = 1.5, shape = 21, stroke = 0.5) +
  scale_x_date("Date", breaks = "1 month", date_labels = "%b %d", limits = as.Date(c("2019-02-28", "2019-11-18")), expand = c(0,0)) +
  scale_y_log10("Density\n(ind./L)") +
  scale_color_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99"), guide = guide_legend(ncol = 4)) +
  scale_fill_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99"), guide = guide_legend(ncol = 4)) +
  coord_cartesian(clip = "off") +
  #facet_grid(.~ Year, scales = "free_y") +
  ggtitle("<b>\\(c\\) Zooplankton:</b> Cladocera spp. (other)") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        plot.title = element_markdown(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  LRV_data %>%
  filter(Group == "Copepoda spp.") %>%
  ggplot(aes(x = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)[DOY], y = value)) +
  geom_path(aes(color = name)) +
  geom_point(aes(fill = name), color = "white", size = 1.5, shape = 21, stroke = 0.5) +
  scale_x_date("Date", breaks = "1 month", date_labels = "%b %d", limits = as.Date(c("2019-02-28", "2019-11-18")), expand = c(0,0)) +
  scale_y_log10("Density\n(ind./L)") +
  scale_color_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  scale_fill_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  coord_cartesian(clip = "off") +
  ggtitle("<b>(d) Zooplankton:</b> Copepoda spp.") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        plot.title = element_markdown(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  LRV_data %>%
  filter(Group == "Ostracoda sp.") %>%
  ggplot(aes(x = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)[DOY], y = value)) +
  geom_path(aes(color = name)) +
  geom_point(aes(fill = name), color = "white", size = 1.5, shape = 21, stroke = 0.5) +
  scale_x_date("Date", breaks = "1 month", date_labels = "%b %d", limits = as.Date(c("2019-02-28", "2019-11-18")), expand = c(0,0)) +
  scale_y_log10("Density\n(ind./L)") +
  scale_color_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  scale_fill_manual(values = c("#3969AC","#F2B701","#E73F74","#7F3C8D","#11A579","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")) +
  coord_cartesian(clip = "off") +
  ggtitle("<b>(e) Zooplankton:</b> Ostracoda sp.") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        plot.title = element_markdown(size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  plot_layout(ncol = 1)
ggsave("Seasonal_dynamics.png", width = 16, height = 22, units = "cm", dpi = 600)
