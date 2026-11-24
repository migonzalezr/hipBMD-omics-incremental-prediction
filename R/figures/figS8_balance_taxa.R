# ======================================================================
# figS8_balance_taxa.R — Figure S8: taxonomic composition of principal balances 1 and 136
#
# Inputs : results/balance1_balance136_species_loadings.csv (script 35)
# Outputs: figures/png/figS8_balance_composition.png, figures/pdf/figS8_balance_composition.pdf
# Paper  : Figure S8
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)
library(patchwork)

source("config/paths.R")
input <- file.path(BASE, "results", "balance1_balance136_species_loadings.csv")
png_dir <- file.path(BASE, "figures", "png")
pdf_dir <- file.path(BASE, "figures", "pdf")
members <- fread(input)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
sides <- c("Denominator", "Numerator")
colors <- c(Denominator = "#685B9D", Numerator = "#118372")
ink <- "#223044"
theme_base <- theme_classic(base_size = 9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        legend.position = "none",
        panel.grid.major.x = element_line(colour = "#E6E9EF", linewidth = .25),
        plot.margin = margin(3, 3, 3, 2))

balance_panel <- function(id, letter) {
  selected <- members[balance == id]
  stopifnot(nrow(selected) == if (id == 1L) 117L else 18L,
            uniqueN(selected$n_total) == 1L, selected$n_total[1] == 1430L)
  phyla <- selected[, .N, by = .(phylum, side)]
  phylum_order <- phyla[, .(total = sum(N)), by = phylum][order(total, phylum), phylum]
  phyla[, phylum := factor(gsub("_", "–", phylum), levels = gsub("_", "–", phylum_order))]
  phyla[, signed_count := ifelse(side == "Denominator", -N, N)]
  count_ticks <- if (id == 1L) seq(-40, 40, by = 20) else seq(-4, 4, by = 2)
  count_limit <- if (id == 1L) 59 else 6.1
  counts <- ggplot(phyla, aes(signed_count, phylum, fill = side)) +
    geom_col(width = .63) +
    geom_vline(xintercept = 0, colour = "#778091", linewidth = .34) +
    geom_text(aes(label = N, hjust = ifelse(signed_count < 0, 1.27, -.27)),
              colour = ink, size = 2.55) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = count_ticks, labels = abs, limits = c(-count_limit, count_limit),
                       expand = expansion(mult = 0)) +
    labs(title = sprintf("%s   Principal balance %d (%d taxa)", letter, id, nrow(selected)),
         x = "Taxa per phylum", y = NULL) +
    theme_base +
    theme(plot.title.position = "plot",
          plot.title = element_text(size = 11, face = "bold", colour = ink, margin = margin(b = 8)))

  prevalence_strip <- function(group, show_axis) {
    all_taxa <- selected[side == group][order(-prevalence_percent, species)]
    top <- head(all_taxa, 10L)
    top[, display := gsub("_", " ", species, fixed = TRUE)]
    top[, row := .N:1L]
    title <- sprintf("%s · %d taxa · loading %s%.5f", group, nrow(all_taxa),
                     if (group == "Denominator") "−" else "+", abs(all_taxa$loading[1]))
    stopifnot(uniqueN(all_taxa$loading) == 1L, all(top$prevalence_percent <= 100))
    ggplot(top, aes(prevalence_percent, row)) +
      geom_segment(aes(x = 0, xend = prevalence_percent, yend = row),
                   colour = colors[[group]], linewidth = .46, alpha = .49) +
      geom_point(colour = colors[[group]], size = 1.75) +
      geom_text(aes(label = sprintf("%.1f%%", prevalence_percent)), hjust = -.20,
                size = 2.35, colour = ink) +
      scale_y_continuous(limits = c(.5, 10.5), breaks = top$row, labels = top$display,
                         expand = expansion(mult = 0)) +
      scale_x_continuous(limits = c(0, 125), breaks = c(0, 50, 100),
                         labels = function(x) paste0(x, "%"), expand = expansion(mult = 0)) +
      labs(title = title, x = if (show_axis) "Detection prevalence" else NULL, y = NULL) +
      theme_base +
      theme(axis.text.y = element_text(face = "italic", size = 8),
            plot.title.position = "plot",
            plot.title = element_text(face = "bold", colour = colors[[group]], size = 8.6,
                                      margin = margin(b = 5))) +
      if (!show_axis) theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.x = element_blank()) else theme()
  }

  (counts / prevalence_strip("Denominator", FALSE) / prevalence_strip("Numerator", TRUE)) +
    plot_layout(heights = c(1.52, 2.06, 2.06))
}

figure <- balance_panel(1L, "a") | plot_spacer() | balance_panel(136L, "b")
figure <- figure + plot_layout(widths = c(1, .055, 1))
ggsave(file.path(png_dir, "figS8_balance_composition.png"), figure,
       width = 7.65, height = 6.8, dpi = 450, bg = "white")
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin" && capabilities("aqua")) {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(pdf_dir, "figS8_balance_composition.pdf"), figure,
       width = 7.65, height = 6.8, device = pdf_device, bg = "white")
