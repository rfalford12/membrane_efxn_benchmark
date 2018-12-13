# @file: plot_protein_design_results.R
# @brief: Plot protein desgin results
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)
library(ggrepel)
library(scales)

workdir <- "/Volumes/ralford/membrane_efxn_benchmark/analysis/protein_design"

df.total.recov <- read.table( paste( workdir, "total_sequence_recovery.dat", sep = "/"), header = T )
trcov <- melt(df.total.recov, id.vars = c("efxn"), value.name = "recovery", variable.name = "solvation" )
df.total.div <- read.table( paste( workdir, "total_sequence_divergence.dat", sep = "/"), header = T )
tdiv <- melt(df.total.div, id.vars = c("efxn"), value.name = "divergence", variable.name = "solvation" )
total.stats <- trcov
total.stats$divergence <- tdiv$divergence
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

design.all.list <- c("r15_all", "m07_all", "m12_all", "m18_all_DOPC_37", "m18_all_DLPC_37" )
total.stats.all <- total.stats[ which( is.element( total.stats$efxn, design.all.list ) ), ]

r <- ggplot( data = total.stats.all, aes( x = recovery, y = divergence, fill = lipid.types ) ) + 
  geom_point( aes( color = lipid.types ), size = 2 ) + 
  geom_label_repel( aes( label = efxn.version ), force = 1 ) + 
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) + 
  geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) +
  scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5), expand = c(0,0) ) + 
  scale_y_continuous( "KL Divergence", expand = c(0.1, 0.1) ) + 
  scale_fill_manual( values = c("#fbb4ae", "#b3cde3", "#bdbdbd" ) ) + 
  scale_color_manual( values = c("#fbb4ae", "#b3cde3", "#bdbdbd" ) ) + 
  facet_wrap(  ~ solvation, scales = "free", ncol = 3, nrow = 2 ) + 
  theme( legend.position = "none" ) 
print(r)

# df.aa.recov <- read.table( paste( workdir, "per_amino_acid_recovery.dat", sep = "/"), header = T )
# aa.rcov <- melt(df.aa.recov, id.vars = c("efxn", "residue", "category"), value.name = "recovery", variable.name = "solvation" )
# df.aa.div <- read.table( paste( workdir, "per_amino_acid_divergence.dat", sep = "/"), header = T )
# aa.div <- melt(df.aa.div, id.vars = c("efxn", "residue", "category"), value.name = "kl.divergence", variable.name = "solvation" )
# per.aa.stats <- aa.rcov
# per.aa.stats$kl.divergence <- aa.div$kl.divergence
# 
# design.all.list <- c("r15_all", "m07_all", "m12_all", "m18_all_DOPC_37")
# per.aa.stats.all <- per.aa.stats[ which( is.element( per.aa.stats$efxn, design.all.list ) ), ]
# 
# # (not sure this is the right plot for this example)
# p <- ggplot( data = per.aa.stats.all[ which( per.aa.stats.all$solvation != "surface"),], aes( x = recovery, y = kl.divergence, fill = category ) ) +
#   geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) + 
#   geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) + 
#   geom_point( size = 0.75 ) + 
#   geom_label( aes( label = residue ), size = 2.5, label.padding = unit(0.1, "lines"), force = 0.2, segment.size = 0.35 ) +
#   scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5), expand = c(0.05,0.05) ) +
#   scale_y_continuous( "KL Divergence", limits = c(-4, 3) ) +
#   scale_fill_brewer( palette = "Pastel1") + 
#   facet_grid( efxn ~ solvation, scales = "free_x" ) + 
#   theme_bw() + 
#   theme(  legend.position = "none",  
#           axis.text = element_text( size = 10 ), 
#           text = element_text( size = 10 ),
#           axis.line.x = element_line( size = 0.35 ),
#           axis.line.y = element_line( size = 0.35 ), 
#           axis.ticks.x = element_line( size = 0.35 ), 
#           axis.ticks.y = element_line( size = 0.35) )  
# print(p)
# 
# # Plot overall sequence recover vs. divergence
# df.class.recov <- read.table( paste( workdir, "per_amino_acid_class_recovery.dat", sep = "/"), header = T )
# class.rcov <- melt(df.class.recov, id.vars = c("efxn", "class"), value.name = "recovery", variable.name = "solvation" )
# df.class.div <- read.table( paste( workdir, "per_amino_acid_class_divergence.dat", sep = "/"), header = T )
# class.div <- melt(df.class.div, id.vars = c("efxn", "class"), value.name = "kl.divergence", variable.name = "solvation" )
# per.class.stats <- class.rcov
# per.class.stats$kl.divergence <- class.div$kl.divergence
# 
# design.all.list <- c("r15_all", "m07_all", "m12_all", "m18_all_DOPC_37")
# per.class.stats.all <- per.class.stats[ which( is.element( per.class.stats$efxn, design.all.list ) ), ]
# 
# # (not sure this is the right plot for this example)
# q <- ggplot( data = per.class.stats.all[ which( per.class.stats.all$solvation != "surface"),], aes( x = recovery, y = kl.divergence, fill = class ) ) +
#   geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) + 
#   geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) + 
#   scale_fill_brewer( palette = "Pastel1") + 
#   geom_point( size = 0.75 ) + 
#   geom_label_repel( aes( label = class ), size = 2.5, label.padding = unit(0.1, "lines"), force = 0.2, segment.size = 0.35 ) +
#   scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.75), expand = c(0.05,0.05) ) +
#   scale_y_continuous( "KL Divergence", limits = c(-20, 10) ) +
#   facet_grid( efxn ~ solvation, scales = "free_x" ) + 
#   theme_bw() + 
#   theme(  legend.position = "none",  
#           axis.text = element_text( size = 10 ), 
#           text = element_text( size = 10 ),
#           axis.line.x = element_line( size = 0.35 ),
#           axis.line.y = element_line( size = 0.35 ), 
#           axis.ticks.x = element_line( size = 0.35 ), 
#           axis.ticks.y = element_line( size = 0.35) )  
# print(q)
# 
# save_plot( "~/Desktop/per_class_seqrecov.png", q, units = "in", base_width = 8, base_height = 6.5)

## TODO: need to do a bunch of specific subset comparisons

## Analyze the dependence of recovery of human sequences on lipid composition
## TODO - maybe plot these on the same plot? 
human.overall.stats <- total.stats[ which( total.stats$subset == "humans" & total.stats$efxn.version == "m18" ),] 
t <- ggplot( data = human.overall.stats, aes( x = lipid.types, y = recovery, fill = efxn.version ) ) + 
  geom_point( aes( color = lipid.types ), size = 2 ) + 
  #geom_label_repel( aes( label = lipid.types ), force = 1 ) + 
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) + 
  geom_hline( yintercept = 0, color = "gray60", linetype = "dashed" ) +
 #scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5), expand = c(0,0) ) + 
#  scale_y_continuous( "KL Divergence", expand = c(0.1, 0.1) ) + 
  facet_wrap(  ~ solvation, scales = "free", ncol = 3, nrow = 2 ) + 
  theme( legend.position = "none" ) 
print(t)
