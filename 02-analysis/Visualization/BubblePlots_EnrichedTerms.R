# aim, create a bubble plot for the clustering results from David
# load dependencies

library(grid)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(ggplot2)  
library(shades)
library(gridExtra)
library(here)


# =======================================
# =======================================
# =======================================
sink("sessionInfo_BubbplePlots2.txt")
sessionInfo()
sink()

# define color range, red and blue
colors2 <- brewer.pal(9,"Set1") 
myColors <- colors2[2:1]

#  ALL AGES
# load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
file <- here("Data","DavidOutput","DAVID_Results_Aug2022_ForPlotting_AllAges.txt")
df1 <- read.table(file, sep="\t", header = TRUE)

# separate the data frame in up and down (needed for the pvalue shading)
df1_UP <- df1[df1$Direction==1,] 
df1_DOWN <- df1[df1$Direction==-1,]

data <- rep(c(" "), times = c(1))
groups <- matrix(data, ncol = 1, byrow = TRUE)
df1_UP$groups <- groups

data <- rep(c("Cluster1"), times = c(6))
groups <- matrix(data, ncol = 1, byrow = TRUE)

df1_DOWN$groups <- groups

plot1 <- ggplot2::ggplot(data = df1_UP, mapping = aes(x=FoldEnrichment, #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = Count,
                                                        color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 10) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(1,10),breaks=c(3,6,9),range = c(5,15)) + # set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Reds"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("AllAges")+
  xlim(1,11.75)

plot2 <- ggplot2::ggplot(data = df1_DOWN, mapping = aes(x=FoldEnrichment, #x axis
                                             y=reorder(Term, -Order), # y axis
                                             size = Count,
                                             color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_continuous(limits=c(1,10),breaks=c(3,6,9),range = c(5,15)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Blues"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  xlim(1,11.75)+
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free")

# grid.arrange(plot1, plot2, nrow = 2) #stack the two plots together, some additional cosmetics can be done in illustrator

grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size="last"))

##################################
# Shared between P24, P30 and P90
#################################
file <- here("Data","DavidOutput","DAVID_Results_Aug2022_ForPlotting_P24-P30-P90.txt")
df2 <- read.table(file, sep="\t", header = TRUE)

# separate the data frame in up and down (needed for the pvalue shading)
df2_UP <- df2[df2$Direction==1,] 
df2_DOWN <- df2[df2$Direction==-1,]

data <- rep(c("Cluster 1","Cluster 2","unclustererd"), times = c(4,13,7))
groups <- matrix(data, ncol = 1, byrow = TRUE)
df2_UP$groups <- groups

data <- rep(c("unclustered","Cluster1"), times = c(3,4))
groups <- matrix(data, ncol = 1, byrow = TRUE)

df2_DOWN$groups <- groups

plot3 <- ggplot2::ggplot(data = df2_UP, mapping = aes(x=FoldEnrichment, #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = Count,
                                                      color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 10) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(1,50),breaks=c(5,25,45),range = c(2,12)) + # set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Reds"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("P24/P30/P90")+
  xlim(1,6)
  
plot4 <- ggplot2::ggplot(data = df2_DOWN, mapping = aes(x=FoldEnrichment, #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = Count,
                                                        color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_continuous(limits=c(1,50),breaks=c(5,25,45),range = c(2,12)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Blues"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  xlim(1,6)+
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free")

# grid.arrange(plot3, plot4, nrow = 3, layout_matrix = cbind(c(1,1,2))) #stack the two plots together, some additional cosmetics can be done in illustrator

grid.draw(rbind(ggplotGrob(plot3), ggplotGrob(plot4), size="last"))

####################
##### P24 only
####################
file <- here("Data","DavidOutput","DAVID_Results_Aug2022_ForPlotting_P24only.txt")
df3 <- read.table(file, sep="\t", header = TRUE)

df3_UP <- df3[df3$Direction==1,] 
df3_DOWN <- df3[df3$Direction==-1,]

data <- rep(c(" "), times = c(3))
groups <- matrix(data, ncol = 1, byrow = TRUE)
df3_UP$groups <- groups

data <- rep(c("Cluster 1"," "), times = c(2,6))
groups <- matrix(data, ncol = 1, byrow = TRUE)

df3_DOWN$groups <- groups

plot5 <- ggplot2::ggplot(data = df3_UP, mapping = aes(x=FoldEnrichment, #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = Count,
                                                      color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 10) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(1,80),breaks=c(20,40,60),range = c(2,17)) + # set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Reds"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("P24only")+
  xlim(1,8)

plot6 <- ggplot2::ggplot(data = df3_DOWN, mapping = aes(x=FoldEnrichment, #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = Count,
                                                        color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_continuous(limits=c(1,80),breaks=c(20,40,60),range = c(2,17)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Blues"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  xlim(1,8)+
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free")

# grid.arrange(plot5, plot6, nrow = 3, layout_matrix = cbind(c(1,2,2))) #stack the two plots together, some additional cosmetics can be done in illustrator

grid.draw(rbind(ggplotGrob(plot5), ggplotGrob(plot6), size="last"))


#####################
##### P90 only
##################
file <- here("Data","DavidOutput","DAVID_Results_Aug2022_ForPlotting_P90only.txt")
df5 <- read.table(file, sep="\t", header = TRUE)

df5_UP <- df5[df5$Direction==1,] 
df5_DOWN <- df5[df5$Direction==-1,]

data <- rep(c("Cluster 1", "unclustered"), times = c(17,5))
groups <- matrix(data, ncol = 1, byrow = TRUE)
df5_UP$groups <- groups

data <- rep(c("Cluster 2","Cluster 3","Cluster 4","Cluster 5", "Cluster 6","unclustered"), times = c(10,5,7,3,4,3))
groups <- matrix(data, ncol = 1, byrow = TRUE)

df5_DOWN$groups <- groups

plot9 <- ggplot2::ggplot(data = df5_UP, mapping = aes(x=FoldEnrichment, #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = Count,
                                                      color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 10) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(1,200),breaks=c(50,100,150),range = c(1,13)) + # set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Reds"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ .,space="free", scales = "free") +
  ggtitle("P90only")+
  xlim(0,10.75)


plot10 <- ggplot2::ggplot(data = df5_DOWN, mapping = aes(x=FoldEnrichment, #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = Count,
                                                        color= PValue)) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_continuous(limits=c(1,200),breaks=c(50,100,150),range = c(1,13)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Blues"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  xlim(0,10.75)+
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., space ="free", scales = "free")


# grid.arrange(plot9, plot10, nrow = 3, layout_matrix = cbind(c(1,2,2)))

grid.draw(rbind(ggplotGrob(plot9), ggplotGrob(plot10), size="last"))


