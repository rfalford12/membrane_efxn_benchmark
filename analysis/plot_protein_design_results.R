# @file: plot_protein_design_results.R
# @brief: Plot protein desgin results
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)
library(ggrepel)
library(scales)
library(dplyr)
library(viridisLite)
library(viridis)

# Current working directory
workdir <- "/Volumes/ralford/membrane_efxn_benchmark/analysis/protein_design"

# Files containing all of the data
total.seq.recov <- "total_sequence_recovery.dat"
total.kl.divergence <- "total_sequence_divergence.dat"
per.aa.seq.recov <- "per_amino_acid_recovery.dat"
per.aa.kl.divergence <- "per_amino_acid_divergence.dat"
per.class.seq.recov <- "per_amino_acid_class_recovery.dat"
per.class.kl.divergence <- "per_amino_acid_class_divergence.dat"

# Some predefined target subsets
all.targets <-  c( "m07_all", "m12_all", "r15_all", "m18_all_DLPC_37", "m18_all_DOPC_37")
bacteria.targets <- c( "m07_bacteria", "m12_all", "r15_bacteria", "m18_bacteria_DLPE", "m18_bacteria_DOPC" )
ecoli.targets <- c( "m07_ecoli", "m12_ecoli", "r15_ecoli", "m18_ecoli_DLPE", "m18_ecoli_DLPE" )
human.targets <- c( "m07_humans", "m12_humans", "r15_humans", "m18_humans_DLPC", "m18_humans_DOPC" )
eukaryote.targets <- c( "m07_eukaryotes", "m12_eukaryotes", "r15_eukaryotes", "m18_eukaryotes_DLPC", "m18_eukaryotes_DOPC")

load.design.performance.data <- function( recov.fn, div.fn, run.list ) {
  # Read the recovery and divergence data into melted dataframes
  df.total.recov <- read.table( paste( workdir, recov.fn, sep = "/"), header = T )
  trcov <- melt(df.total.recov, id.vars = c("efxn"), value.name = "recovery", variable.name = "solvation" )
  df.total.div <- read.table( paste( workdir, div.fn, sep = "/"), header = T )
  tdiv <- melt(df.total.div, id.vars = c("efxn"), value.name = "divergence", variable.name = "solvation" )
  
  # Merge into a single dataframe
  total.stats <- trcov
  total.stats$divergence <- tdiv$divergence

  # Use the "tag" to log the energy function version, species/taxonomy, and lipid types
  efxn.versions <- c()
  subset <- c()
  lipid.types <- c()
  for (efxn.type in total.stats$efxn) {
    efxn.array <- strsplit(efxn.type, "_")
    efxn.versions <- c( efxn.versions, efxn.array[[1]][1] )
    subset <- c( subset, efxn.array[[1]][2] )
    if ( length(efxn.array[[1]]) > 2 ) {
      lipid.types <- c( lipid.types, efxn.array[[1]][3])
    } else {
      lipid.types <- c( lipid.types, "none")
    }
  } 
  total.stats$efxn.version <- efxn.versions
  total.stats$subset <- subset
  total.stats$lipid.types <- lipid.types
  
  # Pull specific energy function versions for this data frame
  total.stats.subset <- total.stats[ which( is.element( total.stats$efxn, run.list ) ), ]
  return(total.stats.subset)
}

plot.overall.design.performance <- function( df ) {
  design.performance <- ggplot( data = df, aes( x = divergence, y = recovery, fill = lipid.types ) ) + 
    geom_point( aes( color = lipid.types ), size = 2 ) + 
    geom_label_repel( aes( label = efxn.version ), force = 1 ) + 
    geom_hline( yintercept = 0.05, color = "gray60", linetype = "dashed" ) + 
    geom_vline( xintercept = 0, color = "gray60", linetype = "dashed" ) +
    scale_y_continuous( "Sequence Recovery (%)", limits = c(0, NA), expand = c(0.05,0.05) ) + 
    scale_x_continuous( "KL Divergence", expand = c(0.1, 0.1) ) + 
    scale_fill_manual( values = c("#fbb4ae", "#b3cde3", "#bdbdbd" ) ) + 
    scale_color_manual( values = c("#fbb4ae", "#b3cde3", "#bdbdbd" ) ) + 
    facet_wrap(  ~ solvation, scales = "free", ncol = 3, nrow = 2 ) + 
    theme( legend.position = "none" ) 
  return(design.performance)
}

# Plot the overall design performance for all five subsets
all.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, all.targets )
all.targets.overall.design <- plot.overall.design.performance( all.targets.df )
print(all.targets.overall.design)
save_plot( paste( workdir, "overall_design_performance.png", sep= "/"), all.targets.overall.design, units = "in", base_width = 8, base_height = 5 )

bacteria.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, bacteria.targets )
bacteria.targets.overall.design <- plot.overall.design.performance( bacteria.targets.df)

ecoli.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, ecoli.targets )
ecoli.targets.overall.design <- plot.overall.design.performance( ecoli.targets.df )

eukaryote.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, eukaryote.targets )
eukaryote.targets.overall.design <- plot.overall.design.performance( eukaryote.targets.df )

human.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, human.targets )
human.targets.overall.design <- plot.overall.design.performance( human.targets.df )

# Per AA design performance
df.aa.recov <- read.table( paste( workdir, "per_amino_acid_recovery.dat", sep = "/"), header = T )
aa.rcov <- melt(df.aa.recov, id.vars = c("efxn", "residue", "category"), value.name = "recovery", variable.name = "solvation" )
df.aa.div <- read.table( paste( workdir, "per_amino_acid_divergence.dat", sep = "/"), header = T )
aa.div <- melt(df.aa.div, id.vars = c("efxn", "residue", "category"), value.name = "kl.divergence", variable.name = "solvation" )
per.aa.stats <- aa.rcov
per.aa.stats$kl.divergence <- aa.div$kl.divergence

design.all.list <- c("r15_all", "m07_all", "m12_all", "m18_all_DOPC_37")
per.aa.stats.all <- per.aa.stats[ which( is.element( per.aa.stats$efxn, design.all.list ) ), ]

per.aa.design.perf <- ggplot( data = per.aa.stats.all[ which( per.aa.stats.all$solvation != "surface"),], aes( x = recovery, y = kl.divergence, fill = category ) ) +
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) +
  geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) +
  geom_point( size = 0.75 ) +
  geom_label( aes( label = residue ), size = 2.5, label.padding = unit(0.1, "lines"), force = 0.2, segment.size = 0.35 ) +
  scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5), expand = c(0.05,0.05) ) +
  scale_y_continuous( "KL Divergence" ) +
  scale_fill_brewer( palette = "Pastel1") +
  facet_grid( efxn ~ solvation, scales = "free_y" ) +
  theme_bw() +
  theme(  legend.position = "right",
          axis.text = element_text( size = 10 ),
          text = element_text( size = 10 ),
          axis.line.x = element_line( size = 0.35 ),
          axis.line.y = element_line( size = 0.35 ),
          axis.ticks.x = element_line( size = 0.35 ),
          axis.ticks.y = element_line( size = 0.35) )

save_plot( paste( workdir, "per_aa_design_perf_all.png", sep= "/"), per.aa.design.perf, units = "in", base_width = 10, base_height = 6)

per.aa.design.perf.density <- ggplot( data = per.aa.stats.all[ which( per.aa.stats.all$solvation != "surface"),], aes( x = recovery, fill = efxn, color = efxn ) ) +
  background_grid() + 
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) +
  geom_density( alpha = 0.1 ) + 
  scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5), expand = c(0.05,0.05) ) +
  scale_y_continuous( "Density", expand = c(0,0), limits = c(0, 6) ) +
  scale_color_viridis( option = "D", end = 0.8, direction = -1, discrete = TRUE ) + 
  scale_fill_viridis(  option = "D", end = 0.8, direction = -1, discrete = TRUE ) + 
  facet_grid( ~ solvation, scales = "free_x" ) +
  theme(  legend.position = "bottom",
          axis.text = element_text( size = 10 ),
          text = element_text( size = 10 ),
          axis.line.x = element_line( size = 0.35 ),
          axis.line.y = element_line( size = 0.35 ),
          axis.ticks.x = element_line( size = 0.35 ),
          axis.ticks.y = element_line( size = 0.35) )
print(per.aa.design.perf.density)
save_plot( paste( workdir, "per_aa_design_distribution.png", sep = "/"), per.aa.design.perf.density, units = "in", base_width = 8, base_height = 4 )

# Plot per-class sequence recover vs. divergence
df.class.recov <- read.table( paste( workdir, "per_amino_acid_class_recovery.dat", sep = "/"), header = T )
class.rcov <- melt(df.class.recov, id.vars = c("efxn", "class"), value.name = "recovery", variable.name = "solvation" )
df.class.div <- read.table( paste( workdir, "per_amino_acid_class_divergence.dat", sep = "/"), header = T )
class.div <- melt(df.class.div, id.vars = c("efxn", "class"), value.name = "kl.divergence", variable.name = "solvation" )
per.class.stats <- class.rcov
per.class.stats$kl.divergence <- class.div$kl.divergence

design.all.list <- c("r15_all", "m07_all", "m12_all", "m18_all_DOPC_37")
per.class.stats.all <- per.class.stats[ which( is.element( per.class.stats$efxn, design.all.list ) ), ]

# (not sure this is the right plot for this example)
per.class.design.perf <- ggplot( data = per.class.stats.all[ which( per.class.stats.all$solvation != "surface"),], aes( x = recovery, y = kl.divergence, fill = class ) ) +
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) +
  geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) +
  scale_fill_brewer( palette = "Pastel1") +
  geom_point( size = 0.75 ) +
  geom_label_repel( aes( label = class ), size = 2.5, label.padding = unit(0.1, "lines"), force = 0.2, segment.size = 0.35 ) +
  scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.6), expand = c(0.05,0.05) ) +
  scale_y_continuous( "KL Divergence" ) +
  facet_grid( efxn ~ solvation, scales = "free_x" ) +
  theme_bw() +
  theme(  legend.position = "right",
          axis.text = element_text( size = 10 ),
          text = element_text( size = 10 ),
          axis.line.x = element_line( size = 0.35 ),
          axis.line.y = element_line( size = 0.35 ),
          axis.ticks.x = element_line( size = 0.35 ),
          axis.ticks.y = element_line( size = 0.35) )

save_plot( paste( workdir, "per_class_design_performance.png", sep = "/"), per.class.design.perf, units = "in", base_width = 10, base_height = 6 )

## Analyze the dependence of recovery of human sequences on lipid composition
# human.overall.stats <- total.stats[ which( total.stats$subset == "humans" & total.stats$efxn.version == "m18" & total.stats$lipid.types != "pm" ),]
# human.lipid.comp.design.perf <- ggplot( data = human.overall.stats, aes( x = lipid.types, y = recovery ) ) +
#   background_grid() +
#   geom_bar( position = "dodge", stat = "identity", fill = "gray70" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "overall"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "overall"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "buried"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "buried"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "surface"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "surface"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "lipid_facing"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "lipid_facing"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "interfacial"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "interfacial"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_bar( data = subset( human.overall.stats[ which( human.overall.stats$solvation == "aqueous"), ],
#                            recovery==max( human.overall.stats[ which( human.overall.stats$solvation == "aqueous"), ]$recovery )),
#             aes(lipid.types, recovery), fill = "#fb6a4a", stat = "identity" ) +
#   geom_hline( yintercept = 0.05, color = "gray30", linetype = "dashed" ) +
#   scale_x_discrete( "" ) +
#   scale_y_continuous( "Sequence Recovery (%)", limits = c(0,0.4), expand = c(0,0) ) +
#   facet_wrap( ~ solvation, scales = "free_x" ) +
#   theme( axis.text.x = element_text( angle = 90, hjust = 1 ) )
# save_plot( paste( workdir, "human_lipid_comp_design_perf.png", sep = "/"), human.lipid.comp.design.perf, units = "in", base_width = 5, base_height = 4 )
# 
# # Going to do one alternate version of this plot which compares the thickness
# dlpc.thk <- 15.351
# dmpc.thk <- 17.974
# dppc.thk <- 20.029
# dopc.thk <- 18.639
# popc.thk <- 19.1
# 
# thickness <- c()
# for (lipid.type in human.overall.stats$lipid.types) {
#   if ( lipid.type == "DLPC" ) {
#     thickness <- c( thickness, dlpc.thk )
#   } else if ( lipid.type == "DMPC" ) {
#     thickness <- c( thickness, dmpc.thk )
#   } else if ( lipid.type == "DPPC" ) {
#     thickness <- c( thickness, dppc.thk )
#   } else if ( lipid.type == "DOPC" ) {
#     thickness <- c( thickness, dopc.thk )
#   } else if ( lipid.type == "POPC" ) {
#     thickness <- c( thickness, popc.thk )
#   }
# }
# human.overall.stats$thickness <- thickness
# max.recov.by.solv <- human.overall.stats %>% group_by(solvation) %>% summarise(Max = max(recovery))
# human.lipid.comp.thk.design.perf <- ggplot( data = human.overall.stats, aes( x = thickness, y = recovery ) ) +
#   background_grid() +
#   geom_point() +
#   geom_line() +
#   geom_label_repel( aes( label = lipid.types ), force = 3, fill = "gray80", size = 2.5 ) +
#   geom_hline( yintercept = 0.05, color = "gray30", linetype = "dashed" ) +
#   scale_x_continuous( "Hydrophobic Thickness (Å)", limits = c(15, 21), expand = c(0,0) ) +
#   scale_y_continuous( "Sequence Recovery (%)", limits = c(0,0.4), expand = c(0,0) ) +
#   facet_wrap( ~ solvation )
# save_plot( paste( workdir, "human_lipid_comp_thk_design_perf.png", sep = "/"), human.lipid.comp.thk.design.perf, units = "in", base_width = 5, base_height = 4 )

