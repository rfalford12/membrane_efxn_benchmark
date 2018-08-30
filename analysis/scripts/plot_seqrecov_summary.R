#@file: plot_seqrecov_summary.R
#@brief: Compare sequence recovery given different energy functions
#@author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)

workdir <- "/Users/ralford/research/membrane_efxn_benchmark/analysis/seqrecov-compare"

summary.df <- read.table( paste( workdir, "seqrecov_summary_8-29-18.txt", sep = "/"), header = T )
m07.detailed.df <- read.table( paste( workdir, "mpframework_fa_2007_seqrecov.txt", sep = "/"), header = T )
m07.detailed.df$efxn <- "m07"
m12.detailed.df <- read.table( paste( workdir, "mpframework_smooth_fa_2012_seqrecov.txt", sep = "/"), header = T )
m12.detailed.df$efxn <- "m12" 
r15.detailed.df <- read.table( paste( workdir, "ref2015_seqrecov.txt", sep = "/"), header = T )
r15.detailed.df$efxn <- "r15"
new.detailed.df <- read.table( paste( workdir, "menv-franklin2018_test_wtbe_wt_v9_seqrecov.txt", sep = "/"), header = T )
new.detailed.df$efxn <- "new"

detailed.df <- data.frame()
detailed.df <- rbind( m07.detailed.df, m12.detailed.df, r15.detailed.df, new.detailed.df )
detailed.df <- melt(detailed.df, id=c("Residue", "efxn"))

# Plot Overall Sequence Recovery
p <- ggplot( data = summary.df, aes( x = efxn, y = recovery, fill = efxn ) ) +
  geom_bar( stat = "identity" ) +
  scale_x_discrete( "", expand = c(0,0) ) +
  scale_y_continuous( "% Recovered" ) +
  background_grid() +
  facet_wrap( ~region, scales = "free_y")
print(p)

# Plot Detailed Sequence Recovery by residue
detailed.df <- detailed.df[ which( detailed.df$variable == "Correct.vs.native.core"), ]
core.plot <- ggplot( data = detailed.df, aes( x = efxn, y = value ) ) + 
  geom_bar( stat = "identity" ) +
  facet_wrap( ~ Residue )
print(core.plot)

