#@file: plot_docking_disc.R
#@brief: Compare decoy discrimination from docking runs
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_1.0_three_way_compare"
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
  scale_fill_manual( values = c( "#F8766D", "#619CFF", "#00BA38" ) )
print(docking.plot)

# Pull and compare examples from the large and small set
df.2nq2.m07 <- read.table( paste(dir, "2nq2_large_set_docking_m07.sc", sep = "/" ), header = T)
df.2nq2.m12 <- read.table( paste(dir, "2nq2_large_set_docking_m12.sc", sep = "/" ), header = T)
df.2nq2.r15 <- read.table( paste(dir, "2nq2_large_set_docking_r15.sc", sep = "/" ), header = T)

funnel.2nq2.m07.plot <- ggplot() +   
  geom_point( data = df.2nq2.m07, aes( x = rms, y = total_score ), color = "#F8766D", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(-200,50) ) + 
  background_grid()

funnel.2nq2.m12.plot <- ggplot() +   
  geom_point( data = df.2nq2.m12, aes( x = rms, y = total_score ), color = "#619CFF", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(-800,-600) ) + 
  background_grid()

funnel.2nq2.r15.plot <- ggplot() +   
  geom_point( data = df.2nq2.r15, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(-100,100) ) + 
  background_grid()

subplot1 <- plot_grid( funnel.2nq2.m07.plot, funnel.2nq2.m12.plot, funnel.2nq2.r15.plot, ncol = 1, nrow = 3 )
print(subplot1)

df.1afo.m07 <- read.table( paste(dir, "1afo_small_set_docking_m07.sc", sep = "/"), header = T )
df.1afo.m12 <- read.table( paste(dir, "1afo_small_set_docking_m12.sc", sep = "/"), header = T )
df.1afo.r15 <- read.table( paste(dir, "1afo_small_set_docking_r15.sc", sep = "/"), header = T )

funnel.1afo.m07.plot <- ggplot() +   
  geom_point( data = df.1afo.m07, aes( x = rms, y = total_score ), color = "#F8766D", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(75, 150)) + 
  background_grid()

funnel.1afo.m12.plot <- ggplot() +   
  geom_point( data = df.1afo.m12, aes( x = rms, y = total_score ), color = "#619CFF", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(0, 100) ) + 
  background_grid()

funnel.1afo.r15.plot <- ggplot() +   
  geom_point( data = df.1afo.r15, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) + 
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 25) ) + 
  scale_y_continuous( "Total Score", expand = c(0,0), limits = c(250, 325) ) + 
  background_grid()

subplot2 <- plot_grid( funnel.1afo.m07.plot, funnel.1afo.m12.plot, funnel.1afo.r15.plot, ncol = 1, nrow = 3 )
print(subplot2)

subplot <- plot_grid( subplot1, subplot2, ncol = 2, nrow = 1 )
composite.plot <- plot_grid( docking.plot, subplot, ncol = 2, nrow = 1 )
print(composite.plot)
save_plot( "~/Desktop/docking.tiff", composite.plot, units = "in", base_width = 9, base_height = 4.25 )

