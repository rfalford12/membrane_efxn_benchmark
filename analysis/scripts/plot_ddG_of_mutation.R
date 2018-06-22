#@file: plot_ddG_of_mutation.R
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Plot ddG of mutation from benchmarks

library(cowplot)

dir <- "/Volumes/ralford/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
df.OmpLA <- read.table( paste(dir, "ddG_prediction_OmpLA.dat", sep = "/"), header = T )
df.OmpLA_aro <- read.table( paste(dir, "ddG_prediction_OmpLA_aro.dat", sep = "/"), header = T)
df.PagP <- read.table( paste(dir, "ddG_prediction_PagP.dat", sep = "/"), header = T)

OmpLA.correl.plot <- ggplot( data = df.OmpLA, aes( x = experimental_ddG, y = predicted_ddG, color = efxn ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) +
  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(-4,4), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", limits = c(-4,6), expand = c(0,0) ) + 
  background_grid()
print(OmpLA.correl.plot)

OmpLA_aro.correl.plot <- ggplot( data = df.OmpLA_aro, aes( x = experimental_ddG, y = predicted_ddG, color = efxn ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) + 
  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", expand = c(0,0), limits = c(-4, 1) ) + 
  scale_y_continuous( "Predicted (REU)", expand = c(0,0), limits = c(-4, 6) ) + 
  background_grid()
print(OmpLA_aro.correl.plot)

PagP.correl.plot <- ggplot( data = df.PagP, aes( x = experimental_ddG, y = predicted_ddG, color = efxn ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) + 
  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(-4, 4), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", limits = c(-4, 6), expand = c(0,0) ) + 
  background_grid()
print(PagP.correl.plot)

composite <- plot_grid( OmpLA.correl.plot, PagP.correl.plot, OmpLA_aro.correl.plot, ncol = 3, nrow = 1 )
save_plot( "~/Desktop/ddG_of_mutation_B3.tiff", composite, units = "in", base_width = 8, base_height = 3 )

### New data
df.OmpLA.new <- read.table( paste(dir, "meb2018_v1_OmpLA.dat", sep = "/"), header = T )
df.OmpLA_aro.new <- read.table( paste(dir, "meb2018_v1_OmpLA_aro.dat", sep = "/"), header = T)
df.PagP.new <- read.table( paste(dir, "meb2018_v1_PagP.dat", sep = "/"), header = T)

OmpLA.correl.plot.new <- ggplot( data = df.OmpLA.new, aes( x = experimental_ddG, y = predicted_ddG*0.001, color = class ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) +
#  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(-4,4), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", expand = c(0,0), limits = c(-3, 2) ) + 
  background_grid()
print(OmpLA.correl.plot.new)

OmpLA_aro.correl.plot.new <- ggplot( data = df.OmpLA_aro.new, aes( x = experimental_ddG, y = predicted_ddG*0.001, color = Mut ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) + 
#  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", expand = c(0,0), limits = c(-4, 1) ) + 
  scale_y_continuous( "Predicted (REU)", expand = c(0,0), limits = c(-2,2) ) + 
  background_grid()
print(OmpLA_aro.correl.plot.new)

PagP.correl.plot.new <- ggplot( data = df.PagP.new, aes( x = experimental_ddG, y = predicted_ddG*0.001, color = class ) ) + 
  geom_hline( yintercept = 0, color = "gray80" ) + 
  geom_vline( xintercept = 0, color = "gray80" ) + 
  geom_abline( color = "gray60" ) + 
  geom_point() + 
  theme( legend.position = "none" ) + 
#  scale_color_manual( values = c( "#F8766D", "#619CFF" ) ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(-4, 4), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", expand = c(0,0), limits = c(-1,1) ) + 
  background_grid()
print(PagP.correl.plot.new)

composite.new <- plot_grid( OmpLA.correl.plot.new, PagP.correl.plot.new, OmpLA_aro.correl.plot.new, ncol = 3, nrow = 1 )
save_plot( "~/Desktop/ddG_of_mutation_new.tiff", composite.new, units = "in", base_width = 8, base_height = 3 )


