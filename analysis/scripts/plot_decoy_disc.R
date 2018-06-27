#@file: plot_decoy_disc.R
#@brief: Compare decoy discrimination from refinement runs
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)

dir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/batch_2.0_five_way_compare/decoy-disc"

yarov.yaravoy.set <- read.table( paste( dir, "decoy_disc_large_set.disc", sep = "/" ), header = T )
yarov.yaravoy.set$category <- "Ayarov-yaravoy"
dutagaci.set <- read.table( paste( dir, "decoy_disc_small_set.disc", sep = "/" ), header = T )
dutagaci.set$category <- "Bdutagaci" 

combined.df <- data.frame()
combined.df <- rbind( combined.df, yarov.yaravoy.set )
combined.df <- rbind( combined.df, dutagaci.set )

max.rms.plot <- ggplot( data = combined.df, aes( x = efxn, y = SampledRMS, fill = efxn, alpha = category ) ) +
  geom_boxplot() + 
  scale_x_discrete( "" ) + 
  scale_y_continuous( "Max RMS for Lowest Scoring 2%", limits = c(0,20), expand = c(0,0) ) + 
  theme( legend.position = "none" ) + 
  background_grid() + 
  scale_alpha_manual( values = c( 1, 0.4 ) )
print(max.rms.plot)


# Pull and compare examples from the large and small set
df.yy.m07 <- read.table( paste(dir, "vatp_yaravoy_refined_models_m07.sc", sep = "/" ), header = T)
df.yy.m12 <- read.table( paste(dir, "vatp_yaravoy_refined_models_m12.sc", sep = "/" ), header = T)
df.yy.m18.v1 <- read.table( paste(dir, "vatp_yaravoy_refined_models_memb2018_v1.sc", sep = "/" ), header = T)
df.yy.m18.v2 <- read.table( paste(dir, "vatp_yaravoy_refined_models_memb2018_v2.sc", sep = "/" ), header = T)
df.yy.r15 <- read.table( paste(dir, "vatp_yaravoy_refined_models_r15.sc", sep = "/" ), header = T)

funnel.yy.m07.plot <- ggplot() +
  geom_point( data = df.yy.m07, aes( x = rms, y = total_score ), color = "#F8766D", size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score", expand = c(0,0) ) +
  background_grid()

funnel.yy.m12.plot <- ggplot() +
  geom_point( data = df.yy.m12, aes( x = rms, y = total_score ), color = "#619CFF", size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score", expand = c(0,0) ) +
  background_grid()

funnel.yy.m18.v1.plot <- ggplot() +
  geom_point( data = df.yy.m18.v1, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score", expand = c(0,0) ) +
  background_grid()

funnel.yy.m18.v2.plot <- ggplot() +
  geom_point( data = df.yy.m18.v2, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score", expand = c(0,0) ) +
  background_grid()

funnel.yy.r15.plot <- ggplot() +
  geom_point( data = df.yy.r15, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) +
  scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) +
  scale_y_continuous( "Total Score", expand = c(0,0) ) +
  background_grid()

subplot1 <- plot_grid( funnel.yy.m07.plot, funnel.yy.m12.plot, funnel.yy.m18.v1.plot, funnel.yy.m18.v2.plot, funnel.yy.r15.plot, ncol = 1, nrow = 5 )
print(subplot1)

# df.du.m07 <- read.table( paste(dir, "vatp_dutagaci_refined_models_m07.sc", sep = "/"), header = T )
# df.du.m12 <- read.table( paste(dir, "vatp_dutagaci_refined_models_m12.sc", sep = "/"), header = T )
# df.du.r15 <- read.table( paste(dir, "vatp_dutagaci_refined_models_r15.sc", sep = "/"), header = T )
# 
# funnel.du.m07.plot <- ggplot() +   
#   geom_point( data = df.du.m07, aes( x = rms, y = total_score ), color = "#F8766D", size = 0.4 ) + 
#   scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) + 
#   scale_y_continuous( "Total Score", expand = c(0,0)) + 
#   background_grid()
# 
# funnel.du.m12.plot <- ggplot() +   
#   geom_point( data = df.du.m12, aes( x = rms, y = total_score ), color = "#619CFF", size = 0.4 ) + 
#   scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) + 
#   scale_y_continuous( "Total Score", expand = c(0,0) ) + 
#   background_grid()
# 
# funnel.du.r15.plot <- ggplot() +   
#   geom_point( data = df.du.r15, aes( x = rms, y = total_score ), color = "#00BA38", size = 0.4 ) + 
#   scale_x_continuous( "RMS (Å)", expand = c(0,0), limits = c(0, 10) ) + 
#   scale_y_continuous( "Total Score", expand = c(0,0) ) + 
#   background_grid()
# 
# subplot2 <- plot_grid( funnel.du.m07.plot, funnel.du.m12.plot, funnel.du.r15.plot, ncol = 1, nrow = 3 )
# print(subplot2)
# 
# subplot <- plot_grid( subplot2, subplot1, ncol = 2, nrow = 1 )
# composite.plot <- plot_grid( decoy.disc.plot, subplot, ncol = 2, nrow = 1 )
# print(composite.plot)
# save_plot( "~/Desktop/decoy_discrimination.tiff", composite.plot, units = "in", base_width = 9, base_height = 4.25 )
# 

