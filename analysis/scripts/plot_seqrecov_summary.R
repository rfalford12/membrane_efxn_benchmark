#@file: plot_seqrecov_summary.R
#@brief: Compare sequence recovery given different energy functions
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Volumes/ralford/membrane_efxn_benchmark/analysis/batch_2.0_four_way_compare"
overview.df <- read.table( paste( dir, "seqrecov_overview.txt", sep = "/" ), header = T )

p <- ggplot( data = overview.df, aes( x = subset, y = recovery, fill = efxn ) ) + 
  geom_bar( width = 0.7, stat = "identity", position = position_dodge( width = 0.7), color = "gray60" ) + 
  scale_x_discrete( "", expand = c(0,0) ) + 
  scale_y_continuous( "% Recovered", expand = c(0,0), limits = c(0, 0.5) ) + 
  background_grid() + 
  scale_fill_brewer( palette = "Pastel1", "Method")
print(p)

save_plot( "~/Desktop/sequence_recovery.tiff", p, units = "in", base_width = 5, base_height = 3 )