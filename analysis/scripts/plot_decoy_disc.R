#@file: plot_decoy_disc.R
#@brief: Compare decoy discrimination from refinement runs
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
yarov.yaravoy.set <- read.table( paste( dir, "decoy_disc_large_set.disc", sep = "/" ), header = T )
yarov.yaravoy.set$category <- "Ayarov-yaravoy"
dutagaci.set <- read.table( paste( dir, "decoy_disc_small_set.disc", sep = "/" ), header = T )
dutagaci.set$category <- "Bdutagaci" 

combined.df <- data.frame()
combined.df <- rbind( combined.df, yarov.yaravoy.set )
combined.df <- rbind( combined.df, dutagaci.set )

decoy.disc.plot <- ggplot( data = combined.df, aes( x = efxn, y = SampledRMS, fill = efxn, alpha = category ) ) +
  geom_boxplot() + 
  scale_x_discrete( "" ) + 
  scale_y_continuous( "Max RMS for Lowest Scoring 2%", limits = c(0,20), expand = c(0,0) ) + 
  #theme( legend.position = "none" ) + 
  background_grid() + 
  scale_alpha_manual( values = c( 1, 0.4 ) ) + 
  scale_fill_manual( values = c( "#1b9e77", "#d95f02", "#7570b3" ) )
print(decoy.disc.plot)

save_plot( "~/Desktop/decoydisc.png", decoy.disc.plot, units = "in", base_width = 7.25, base_height = 4 )

