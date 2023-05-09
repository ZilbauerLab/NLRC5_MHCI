# PLOTTING COLOUR SCHEMES AND THEMES (NLRC5 PAPER SPECIFIC)
# 12Oct22

library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(reshape)


myTheme_meanmeth_label <- theme(legend.position="bottom", axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_blank(), panel.border=element_blank(), legend.text=element_text(size=rel(2)))
myTheme_meanmeth_meth <- theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

myColors <- c("#01665e","#74AE97","#e6f5d0")
color_possibilities <- c("High","Intermediate","Low")
names(myColors) <- color_possibilities
fillscale_3_cluster <- scale_fill_manual(name="Cluster", values=myColors, drop=T, limits=c("High","Intermediate","Low"), na.value="white")
colscale_3_cluster <- scale_color_manual(name="Cluster", values=myColors, drop=T)

myColors <- c("#01665e","#e6f5d0")
color_possibilities <- c("Other","Low")
names(myColors) <- color_possibilities
fillscale_2_cluster <- scale_fill_manual(name="Cluster", values=myColors, drop=T, limits=c("Other","Low"), na.value="white")
colscale_2_cluster <- scale_color_manual(name="Cluster", values=myColors, drop=T)

myColours <- c("darkred","tomato","grey75")
colour_possibilities <- c("Macroscopic","Microscopic","Normal")
names(myColours) <- colour_possibilities
fillscale_inflammation <- scale_fill_manual(name="Inflammation", values=myColours, drop=T, limits=force, na.value="white")
colscale_inflammation <- scale_color_manual(name="Inflammation", values=myColours, drop=T)

myColours <- c("darkgreen","grey75")
colour_possibilities <- c("Y","N")
names(myColours) <- colour_possibilities
fillscale_biologics <- scale_fill_manual(name="Biologics", values=myColours, drop=T, limits=force, na.value="white")
colscale_biologics <- scale_color_manual(name="Biologics", values=myColors, drop=T)

myColours <- c("navyblue","grey75")
colour_possibilities <- c("Y","N")
names(myColours) <- colour_possibilities
fillscale_surgery <- scale_fill_manual(name="Surgery", values=myColours, drop=T, limits=force, na.value="white")
colscale_surgery <- scale_color_manual(name="Surgery", values=myColours, drop=T)

myColours <- c("orchid4","grey75")
colour_possibilities <- c("Y","N")
names(myColours) <- colour_possibilities
fillscale_immunosuppressor <- scale_fill_manual(name="AZA", values=myColours, drop=T, limits=force, na.value="white")
colscale_immunosuppressor <- scale_color_manual(name="AZA", values=myColours, drop=T)

myColours <- c("slateblue4","grey75")
colour_possibilities <- c("Y","N")
names(myColours) <- colour_possibilities
fillscale_perianal <- scale_fill_manual(name="Perianal Disease", values=myColours, drop=T, limits=force, na.value="white")
colscale_perianal <- scale_color_manual(name="Perianal Disease", values=myColours, drop=T)

myColours <- c("black","grey30","grey60","grey96")
colour_possibilities <- c("severe","moderate","mild","control")
names(myColours) <- colour_possibilities
fillscale_severity <- scale_fill_manual(name="Disease Severity", values=myColours, drop=T, limits=force, na.value="white")
colscale_severity <- scale_color_manual(name="Disease Severity", values=myColours, drop=T)


