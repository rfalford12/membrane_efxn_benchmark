#@file: plot_seqrecov_summary.R
#@brief: Compare sequence recovery given different energy functions
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
overview.df <- read.table( paste( dir, "seqrecov_overview.txt", sep = "/" ), header = T )

p <- ggplot( data = overview.df, aes( x = subset, y = recovery, fill = efxn) ) + 
  geom_bar( stat = "identity", position = "dodge" ) + 
  scale_x_discrete( "", expand = c(0,0) ) + 
  scale_y_continuous( "% Recovered", expand = c(0,0), limits = c(0, 0.5) ) + 
  theme( legend.position = "none" ) + 
  scale_fill_manual( values = c( "#F8766D", "#619CFF", "#00BA38" ) ) + 
  background_grid()
print(p)

save_plot( "~/Desktop/sequence_recovery.tiff", p, units = "in", base_width = 3.5, base_height = 3 )