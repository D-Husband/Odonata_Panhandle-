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

plot1 <- plot(wetlandsK3, data=comm1s, xlab="Dim1 (25.3%)", ylab="Dim2 (18.2%)", title = "Wetlands")
plot2 <- plot(odesK3, data=odes, xlab="Dim1 (28.5%)", ylab="Dim2 (13.3%)", title = "Odonates")
grid.arrange(plot1, plot2, ncol=2)

#####################
# Wetlands + odonates
#####################

# CCA: 
ccamodel <- cca(odes~.,comm1)
ccamodel
plot(ccamodel, type = "p")
# Permutation tests:
# for model:
A <- anova(ccamodel, test = "F")
A
# Model is significant.
# for environmental variables:
B <- anova(ccamodel, by="terms")
B
# SO2, NTU, Veg, Na significant (others not) 
# Redo CCA with only significant variables:
cca2 <- cca(odes ~ SO2 + NTU + Veg + Na, data = comm1)
summary(cca2)
cca2
#% variance explained by first 3 eigenvalues:
#CCA1:
0.25641/0.4773
#CCA2:
0.12134/0.4773
#CCA3:
0.06213/0.4773
# First 2 account for most
plot(ccamodel, type = "p", main = "CCA, all variables")
plot(cca2, type = "p", main = "CCA, significant variables only") 
legend('topleft', col=1:2, pch=c(20, 3), legend=c("sites", "species"))
# envfit:
plot(ccamodel, type = "p")
fit <- envfit(ccamodel, comm1, perm = 999, display = "lc")
plot(fit, p.max = 0.05, col = "red")
# Reiterates previous findings and plots

# BIOENV:
bioenv(odes, comm1s, method = "spearman", index = "euclidean")
# Salinity is abiotic variable with best model (corr = 0.2988)

# Examining odonate patterns by groupings:

hp <- comm$Hydropd
CCA <- cca(odes ~ hp)
CCA
type_num <- as.numeric(as.factor(comm$Hydropd))
ordiplot(CCA, display=c('si', 'cn'), type = 'n', main="Hydroperiod")
points(CCA, display='si', col = type_num, pch = 20)
legend('bottomright', col=1:2, pch=20, legend=levels(as.factor(comm$Hydropd)))
# axes based on CCA1 coordinates and first unconstrained eigenvalue CA1
# There is separation in the odonate assemblage based on hydroperiod.

sal <- comm$Salinity
CCA_sal <- cca(odes ~ sal)
CCA_sal
type_num3 <- as.numeric(as.factor(comm$Salinity))
ordiplot(CCA_sal, display=c('si', 'cn'), type = 'n', main = "Salinity")
points(CCA_sal, display='si', col = type_num3, pch = 20)
legend('bottomright', col=1:3, pch=20, legend=levels(as.factor(comm$Salinity)))
# Some separation

vegetation <- comm$Vegetation
CCA_veg <- cca(odes ~ vegetation)
CCA_veg
type_num6 <- as.numeric(as.factor(comm$Vegetation))
ordiplot(CCA_veg, display=c('si', 'cn'), type = 'n', main="Vegetation")
points(CCA_veg, display='si', col = type_num6, pch = 20)
legend('bottomright', col=1:3, pch=20, legend=levels(as.factor(comm$Vegetation)))
# Some separation

turb <- comm$Turbid
CCA_turb <- cca(odes ~ turb)
CCA_turb
type_num5 <- as.numeric(as.factor(comm$Turbid))
ordiplot(CCA_turb, display=c('si', 'cn'), type = 'n', main="Turbidity")
points(CCA_turb, display='si', col = type_num5, pch = 20)
legend('bottomright', col=1:3, pch=20, legend=levels(as.factor(comm$Turbid)))
# No separation

par(mfrow = c(2, 2))
ordiplot(CCA, display=c('si', 'cn'), type = 'n', main="Hydroperiod")
points(CCA, display='si', col = type_num, pch = 20)
legend('bottomright', col=c(1,2), pch=20, legend=levels(as.factor(comm$Hydropd)))
ordiplot(CCA_sal, display=c('si', 'cn'), type = 'n', main = "Salinity")
points(CCA_sal, display='si', col = c(1,2,4), pch = 20)
legend('topright', col=c(1,2,4), pch=20, legend=c("medium", "low", "high"))
ordiplot(CCA_veg, display=c('si', 'cn'), type = 'n', main="Vegetation")
points(CCA_veg, display='si', col = c(1,2,4), pch = 20)
legend('topright', col=c(1,2,4), pch=20, legend=c("some", "unvegetated", "vegetated"))
ordiplot(CCA_turb, display=c('si', 'cn'), type = 'n', main="Turbidity")
points(CCA_turb, display='si', col = c(1,2,4), pch = 20)
legend('topleft', col=c(1,2,4), pch=20, legend=c("clear", "medium", "turbid"))

