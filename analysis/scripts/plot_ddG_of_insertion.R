#@file: plot_ddG_of_insertion.R
#@brief: Plot ddG of insertion as part of the benchmark
#@author: Rebecca Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
df.neutral <- read.table( paste( dir, "ddG_of_insertion.dat", sep = "/"), header = T )
df.pH <- read.table( paste(dir, "ddG_of_pH_insertion.dat", sep = "/"), header = T )

ddG.of.insertion.plot <- ggplot( data = df.neutral, aes( x = exp, y = ddG_from_inserted, color = efxn) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_point( size = 2 ) + 
  geom_abline( slope = 2.7412, intercept = -2.49, color = "#F8766D", linetype = "dashed" ) + 
  geom_abline( slope = 0.522, intercept = -3.099, color = "#619CFF", linetype = "dashed" ) + 
  theme( legend.position = "none" ) +
  scale_color_brewer( palette = "Set1" ) +  
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(-2,3), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)",  limits = c(-8, 4), expand = c(0,0) ) + 
  background_grid()
print(ddG.of.insertion.plot)

ddG.of.pH.insertion.plot <- ggplot( data = df.pH, aes( x = exp, y = ddG_from_inserted, color = efxn) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_point( size = 2 ) + 
  theme( legend.position = "none" ) +
  scale_color_brewer( palette = "Set1" ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(0,3), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", limits = c(-10, 15), expand = c(0,0) ) + 
  background_grid()
print(ddG.of.pH.insertion.plot)

composite <- plot_grid( ddG.of.insertion.plot, ddG.of.pH.insertion.plot, nrow = 1, ncol = 2 )
save_plot( "~/Desktop/ddG_of_insertion_B3.tiff", composite, units = "in", base_width = 6, base_height = 2.5)

