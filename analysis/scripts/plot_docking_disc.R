#@file: plot_docking_disc.R
#@brief: Compare decoy discrimination from docking runs
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_2.0_five_way_compare/docking"
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
  background_grid() + 
  scale_fill_brewer( palette = "Set1" ) + 
  scale_alpha_manual( values = c( 1, 0.4 ) ) 
print(docking.plot)
save_plot( "~/Desktop/docking_summary.tiff", docking.plot, units = "in", base_width = 9, base_height = 4 )


# Pull and compare examples from the large and small set
df.2nq2.m07 <- read.table( paste(dir, "2nq2_large_set_docking_m07.sc", sep = "/" ), header = T)
df.2nq2.m12 <- read.table( paste(dir, "2nq2_large_set_docking_m12.sc", sep = "/" ), header = T)
df.2nq2.m18.v1 <- read.table( paste(dir, "2nq2_large_set_docking_memb2018_v1.wts", sep = "/" ), header = T)
df.2nq2.m18.v2 <- read.table( paste(dir, "2nq2_large_set_docking_memb2018_v2.wts", sep = "/" ), header = T)
df.2nq2.r15 <- read.table( paste(dir, "2nq2_large_set_docking_r15.sc", sep = "/" ), header = T)

funnel.2nq2.m07.plot <- ggplot() +
  geom_point( data = df.2nq2.m07, aes( x = Irms, y = I_sc, color = fa_mpenv+fa_mpsolv ), size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Reds" ) + 
  theme( legend.position = "none" )

funnel.2nq2.m12.plot <- ggplot() +
  geom_point( data = df.2nq2.m12, aes( x = Irms, y = I_sc, color = fa_mpenv+fa_mpenv_smooth+fa_mpsolv ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Blues" ) + 
  theme( legend.position = "none" )

funnel.2nq2.m18.v1.plot <- ggplot() +
  geom_point( data = df.2nq2.m18.v1, aes( x = Irms, y = I_sc, color = fa_water_to_bilayer ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Greens" ) + 
  theme( legend.position = "none" )

funnel.2nq2.m18.v2.plot <- ggplot() +
  geom_point( data = df.2nq2.m18.v2, aes( x = Irms, y = I_sc, color = fa_water_to_bilayer), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Purples" ) + 
  theme( legend.position = "none" )

funnel.2nq2.r15.plot <- ggplot() +
  geom_point( data = df.2nq2.r15, aes( x = Irms, y = I_sc ), color = "#ff7f00", size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid()

subplot1 <- plot_grid( funnel.2nq2.m07.plot, funnel.2nq2.m12.plot, funnel.2nq2.m18.v1.plot, funnel.2nq2.m18.v2.plot, funnel.2nq2.r15.plot, ncol = 5, nrow = 1 )

df.1afo.m07 <- read.table( paste(dir, "1afo_small_set_docking_m07.sc", sep = "/"), header = T )
df.1afo.m12 <- read.table( paste(dir, "1afo_small_set_docking_m12.sc", sep = "/"), header = T )
df.1afo.m18.v1 <- read.table( paste(dir, "1afo_small_set_docking_memb2018_v1.wts", sep = "/" ), header = T)
df.1afo.m18.v2 <- read.table( paste(dir, "1afo_small_set_docking_memb2018_v2.wts", sep = "/" ), header = T)
df.1afo.r15 <- read.table( paste(dir, "1afo_small_set_docking_r15.sc", sep = "/"), header = T )

funnel.1afo.m07.plot <- ggplot() +
  geom_point( data = df.1afo.m07, aes( x = Irms, y = I_sc, color = fa_mpenv+fa_mpsolv ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Reds" ) + 
  theme( legend.position = "none" )

funnel.1afo.m12.plot <- ggplot() +
  geom_point( data = df.1afo.m12, aes( x = Irms, y = I_sc, color = fa_mpenv+fa_mpenv_smooth+fa_mpsolv ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Blues" ) + 
  theme( legend.position = "none" )

funnel.1afo.m18.v1.plot <- ggplot() +
  geom_point( data = df.1afo.m18.v1, aes( x = Irms, y = I_sc, color = fa_water_to_bilayer ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Greens" ) + 
  theme( legend.position = "none" )

funnel.1afo.m18.v2.plot <- ggplot() +
  geom_point( data = df.1afo.m18.v2, aes( x = Irms, y = I_sc, color = fa_water_to_bilayer ), size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid() + 
  scale_color_distiller( palette = "Purples" ) + 
  theme( legend.position = "none" )

funnel.1afo.r15.plot <- ggplot() +
  geom_point( data = df.1afo.r15, aes( x = Irms, y = I_sc ), color = "#ff7f00", size = 0.6 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score" ) +
  background_grid()

subplot2 <- plot_grid( funnel.1afo.m07.plot, funnel.1afo.m12.plot, funnel.1afo.m18.v1.plot, funnel.1afo.m18.v2.plot, funnel.1afo.r15.plot, ncol = 5, nrow = 1 )

subplot <- plot_grid( subplot2, subplot1, ncol = 1, nrow = 2 )
save_plot( "~/Desktop/docking_examples.tiff", subplot, units = "in", base_width = 11, base_height = 4 )
