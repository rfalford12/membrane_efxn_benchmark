#@file: plot_docking_disc.R
#@brief: Compare decoy discrimination from docking runs
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Volumes/ralford/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
large.set <- read.table( paste( dir, "docking_large_set.disc", sep = "/" ), header = T )
large.set$category <- "large proteins"
small.set <- read.table( paste( dir, "docking_small_set.disc", sep = "/" ), header = T )
small.set$category <- "small proteins" 

combined.df <- data.frame()
combined.df <- rbind( combined.df, large.set )
combined.df <- rbind( combined.df, small.set )

docking.plot <- ggplot( data = combined.df, aes( x = efxn, y = SampledRMS, fill = efxn, alpha = category ) ) +
  geom_boxplot() + 
  scale_x_discrete( "" ) + 
  scale_y_continuous( "Max RMS for Lowest Scoring 2%", limits = c(0, 2), expand = c(0,0) ) + 
  theme( legend.position = "none" ) + 
  background_grid() + 
  scale_alpha_manual( values = c( 1, 0.4 ) ) + 
  scale_fill_manual( values = c( "#1b9e77", "#7570b3", "#e7298a" ) )
print(docking.plot)

