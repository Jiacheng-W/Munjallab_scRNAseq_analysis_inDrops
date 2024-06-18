library(ggrepel)
DEG$BH <- p.adjust(DEG$pvalue, method = "BH")
out <- volcano_plot(DEG,
                    title = "Group 1 vs Group 2",
                    # file = paste0(output_dir, cell_type_label, "_advanced_vs_middle_volcano.pdf"),
                    padj.thresh = 0.01, lfc.thresh = 2)
# Export volcano plot.
out[[1]]
# Export DEGs
degs <- out[[2]]

# Volcano function --------------------------------------------------------

volcano_plot <- function(data, file = NULL, title = NULL,
         padj.lim = NULL, lfc.lim = NULL, lfc.asymmetric = NULL,
         padj.thresh = 0.01, lfc.thresh = 0.25,
         x_title = "log2FC",y_title = "-log2(P_adj)",
         width=unit(4, "inch"), height=unit(4, "inch")) {
  # Create gene name column
  data$gene <- rownames(data)
  
  # Generate PFC scores
  data$PFC <- -log10(data$BH) * abs(data$avg_log2FC)
  PFC <- unique(data$PFC)
  data$PFC[data$PFC == Inf] <- sort(PFC, partial=length(PFC)-1)[length(PFC)-1]
  
  #Generate log padj column
  data$log_padj <- -log10(data$BH)
  
  # Define limits if not provided
  if (is.null(padj.lim)) {
    log.padj.lim <- unique(data$log_padj)[order(-unique(data$log_padj))][2]
  } else {
    log.padj.lim <- -log10(padj.lim)
  }
  if (is.null(lfc.lim)) {
    lfc.lim <- abs(data[order(-abs(data$avg_log2FC)),"avg_log2FC"][1])
  }
  
  # Generate color column
  data$color <- rep("black", nrow(data))
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC > lfc.thresh), 'color'] <- "red"
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC < -lfc.thresh), 'color'] <- "blue"
  
  # Scale down genes outside of bounds
  data[which(data$log_padj > log.padj.lim), 'log_padj'] <- log.padj.lim
  data[which(data$avg_log2FC > lfc.lim), 'avg_log2FC'] <- lfc.lim
  data[which(data$avg_log2FC < -lfc.lim), 'avg_log2FC'] <- -lfc.lim
  
  # Generate asymmetric lfc limits if necessary
  if (is.null(lfc.asymmetric)) {
    lfc_lims <- c(-lfc.lim, lfc.lim)
  } else {
    lfc_lims <- lfc.asymmetric
  }
  
  # Plot data
  p <-
    ggplot(data,
           aes(
             x = avg_log2FC,
             y = log_padj,
             color = color,
             label = gene,
             size = PFC
           )) +
    theme(panel.background = element_blank(),
          legend.key=element_blank(),
          plot.tag=element_blank(),
          plot.caption=element_blank())+
    geom_hline(
      yintercept = -log10(padj.thresh),
      linetype = 2,
      color = "gray"
    ) +
    geom_vline(xintercept = lfc.thresh,
               linetype = 2,
               color = "gray") +
    geom_vline(
      xintercept = -lfc.thresh,
      linetype = 2,
      color = "gray"
    ) +
    geom_point(aes(size = PFC), alpha = 0.5) +
    scale_color_manual(values = c("red" = "red",
                                  "black" = "black",
                                  "blue" = "blue")) +
    geom_text_repel(
      data = data[which(data$color != 'black'), ],
      inherit.aes = T,
      color = 'black',
      size = 4,
      force = 3
    ) +
    theme(legend.position = "none") +
    labs(title = title,
         x = x_title,
         y = y_title) +
    scale_x_continuous(limits = lfc_lims) +
    scale_y_continuous(limits = c(0, log.padj.lim)) 
  theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12))
  
  # Export plot
  if (is.null(file)) {
    out <- list(one=p, two=data)
    return(out)
  } else {
    set_panel_size(
      p,
      file = file,
      width = width,
      height = height
    )
    return(data)
  }
}