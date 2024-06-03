##########
# Wetlands
##########

library(factoextra)
library(gridExtra)
library(tidyverse)
library(cluster)
library(useful)
library(vegan)
library(labdsv)
library(MASS)
library(MVA)
library(optpart)
library(picante)
library(stats)
library(vegan)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read in environmental data and scale:
comm <- read.csv("Pooled_abiotic.csv", header = TRUE, row.names = 1) 
comm1 <- comm[,9:16]
comm1s <- scale(comm1, center=TRUE, scale=TRUE)
head(comm1s)

# Reminders: conductivity, TDS removed due to high correlations; 
# also removed basin fill and water temperature; veg is
# arcsine sqrt transformed and others are log10-transformed.
# P = PO4

# k-means clustering, for k = 2:12; number of starts can make a big difference
E2 <- eclust(comm1s, "kmeans", k = 2, nstart=25)
E3 <- eclust(comm1s, "kmeans", k = 3, nstart=25)
E4 <- eclust(comm1s, "kmeans", k = 4, nstart=25) # e.g. playa, urban, salina, former salina
# Those four classes do not parse out well.
E5 <- eclust(comm1s, "kmeans", k = 5, nstart=25) 
# Color-blind-friendly palette:
p1 <- fviz_cluster(E2, geom = "point", data = df, ggtheme=theme_minimal()) + ggtitle("k = 2") + scale_color_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73")) + 
  scale_fill_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73"))
p2 <- fviz_cluster(E3, geom = "point",  data = df, ggtheme=theme_minimal()) + ggtitle("k = 3")+ scale_color_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73")) + 
  scale_fill_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73"))
p3 <- fviz_cluster(E4, geom = "point",  data = df, ggtheme=theme_minimal()) + ggtitle("k = 4")+ scale_color_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73")) + 
  scale_fill_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73"))
p4 <- fviz_cluster(E5, geom = "point",  data = df, ggtheme=theme_minimal()) + ggtitle("k = 5")+ scale_color_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73")) + 
  scale_fill_manual(values=c("#999999", "#56B4E9", "#cc79A7", "#E69F00", "#009E73"))
grid.arrange(p1, p2, p3, p4, nrow = 2)
# Results:
# So there are 3 clusters of wetlands based on environmental variables.
# Cluster 1: F2, F3, F4, S3, P10, U1, U3, U4, U10, U11, U13, U16, U19, U20, U21, U23
# Cluster 2: F1, S1, S2, S4, P4, P6, P7, P9, U2, U6, U7, U8, U9, U12, U14, U15, U17, U18, U22, U24, U25, U26
# Cluster 3: P1, P2, P3, P5, P8, P11, P12, P13, U5

# None of these clusters is "homogeneous" in terms of site types (land use) or hydroperiod. 

#PCA:
pca.1 <- princomp(comm1) #centering and scaling are default in princomp
summary(pca.1)
# Examine the factor loadings:
pca.1$loadings
# It takes 3 axes to explain at least 50% of variation in the data. The
# first three axes explain 64.5% of variation in data.
# PC1 strongly positively associated with vegetation, PO4, NTU, and
# negatively associated with SO2, NO2, and salinity
# PC2  strongly positively associated with turbidity, negatively with salinity,
# vegetation, PO4
# PC3  associated positively with NO3, SO2, PO4; negatively with turbidity and 
# salinity
biplot(pca.1$scores, pca.1$loadings, cex=0.5)
fviz_cos2(pca.1, choice="var", axes=1)
fviz_cos2(pca.1, choice="var", axes=2)
fviz_cos2(pca.1, choice="var", axes=3)
# Scree plot:
fviz_eig(pca.1, addlabels=TRUE)

##########
# Odonates
##########

ode.biol <- read.csv("all_years_no_name.csv", header=TRUE, row.names=1)
odes <- ode.biol[,3:39]
head(odes)

# k-means clustering:
odesK3 <- kmeans(x=odes, center = 3, nstart=25) 
# Note: eclust method doesn't work on binary data with lots of 0s
# (cannot rescale a constant/zero column to unit variance)
odesK3
plot(odesK3, data=odes)
#more overlap (increases w/ increasing k)
# So there may be 3 clusters of odonates. Do these clusters 
# correspond to the 3 wetland clusters?
# 3 odonate clusters, sites separate as:
# Cluster 1 - S1, S2
# Cluster 2 - U1, U3, U4, U5, U6, U7, U8, U9, U10, U13, U14, U16, U20, U21, U24, 
# S3, F2, S4, F3, F4, P4, P11, P12
# Cluster 3 - U2, U11, U12, U15, U17, U18, U19, U22, U23, U25, U26, F3, P1, P2, P3,
# P5, P6, P7, P8, P9, P10, P13
# So, NOT the same as the wetland clusters

plot1 <- plot(E3, data=comm1s, xlab="Dim1 (25.3%)", ylab="Dim2 (18.2%)", title = "Wetlands")
plot2 <- plot(odesK3, data=odes, xlab="Dim1 (28.5%)", ylab="Dim2 (13.3%)", title = "Odonates")
grid.arrange(plot1, plot2, ncol=2)

#####################
# Wetlands + odonates
#####################

# RDA:

# Perform Hellinger transformation on odonate presence/absence data:
odes.hell <- hellinger(odes)

RDA1 <- rda(odes.hell~., comm1)
print(RDA1)
# 0.0897/0.4639 = 19.3% of odonate variation explained by environmental variables
# 0.030068/0.0897 = 33.52% of constrained variability is explained by the first axis alone;
# this value is greater than the sum of the variability explained by the first two
# unconstrained axes:
0.06249/0.3742
0.05087/0.3742
0.1669963+0.1359433
# Constrained proportion of variance explained = 19.33%
# Unconstrained proportion of variance explained = 80.66%
# Thus, the environmental variables are not explaining the majority of the variation
# seen in the odonate assemblage data.

# Significance testing for environmental variables:
RDA_A <- anova(RDA1, by="terms")
RDA_A
# SO2, Veg significant (others not) 
# The above is the same as a PERMANOVA:
odes.dist <- vegdist(odes.hell, method = "euclidean")
odes.div <- adonis2(odes.hell~., comm1, permutations = 999, method = "euclidean")
odes.div
# Redo RDA with only significant variables:# Redo RDA with only significant variables:NULL
RDA2 <- rda(odes.hell ~ SO2 + Veg, data = comm1)
RDA2

plot(RDA1, type = "p", main = "RDA, all variables")
plot(RDA2, type = "p", main = "RDA, significant variables only") 
legend('topleft', col=1:2, pch=c(20, 3), legend=c("sites", "species"))

# BIOENV:
bioenv(odes.hell, comm1s, method = "spearman", index = "euclidean")
# Vegetation is abiotic variable with best model (corr = 0.1350)

# Examining odonate patterns by groupings/natural breaks:

hp <- comm$Hydropd
RDA.hp <- rda(odes.hell ~ hp)
RDA.hp
type_num <- as.numeric(as.factor(comm$Hydropd))
ordiplot(RDA.hp, display=c('si', 'cn'), type = 'n', main="Hydroperiod")
points(RDA.hp, display='si', col = type_num, pch = 20)
legend('bottomright', col=1:2, pch=20, legend=levels(as.factor(comm$Hydropd)))
# axes based on RDA1 coordinates and first unconstrained eigenvalue PC1 (because
# only two groupings)
# There is some separation in the odonate assemblage based on hydroperiod.

vegetation <- comm$Vegetation
RDA_veg <- rda(odes.hell ~ vegetation)
RDA_veg
type_num2 <- as.numeric(as.factor(comm$Vegetation))
ordiplot(RDA_veg, display=c('si', 'cn'), type = 'n', main="Vegetation")
points(RDA_veg, display='si', col = type_num2, pch = 20)
legend('bottomright', col=1:3, pch=20, legend=levels(as.factor(comm$Vegetation)))
# Some separation, but much overlap

sulphate <- comm$SO2
RDA_SO2 <- rda(odes.hell ~ sulphate)
RDA_SO2
type_num3 <- as.numeric(as.factor(comm$Sulfate))
ordiplot(RDA_SO2, display=c('si', 'cn'), type = 'n', main="Sulphate Concentration")
points(RDA_SO2, display='si', col = type_num3, pch = 20)
legend('bottomright', col=1:3, pch=20, legend=levels(as.factor(comm$Sulfate)))
# No separation

par(mfrow = c(1, 3))
ordiplot(RDA.hp, display=c('si', 'cn'), type = 'n', main="Hydroperiod")
points(RDA.hp, display='si', col = type_num, pch=20, cex=2)
legend('bottomright', col=1:2, pch=20, cex=1.5, legend=levels(as.factor(comm$Hydropd)))
ordiplot(RDA_veg, display=c('si', 'cn'), type = 'n', main="Vegetation")
points(RDA_veg, display='si', col = type_num2, pch = 20, cex=2)
legend('bottomright', col=1:3, pch=20, cex=1.5, legend=c("some", "bare", "vegetated"))
ordiplot(RDA_SO2, display=c('si', 'cn'), type = 'n', main="Sulphate Concentration")
points(RDA_SO2, display='si', col = type_num3, pch = 20, cex=2)
legend('bottomright', col=1:3, pch=20, cex=1.5, legend=c("high", "low", "medium"))
par(mfrow = c(1,1))


