#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Biodiversity analysis
# Author: Decolle Alicia & Grégoire Noël
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RStudio.Version()
rm(list=ls())

#1. Alpha diversity : (ANOVA + correlation test)
#_______________________________________________

# Packages installation & preparation 
#------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("phyloseq")

#install.packages("lme4")
#install.packages("multcompView")
#install.packages("lme4")
#install.packages("ggeffects")
#install.packages("stargazer")
#install.packages("ggfortify")
#install.packages("lsmeans")
#install.packages("msm")
#install.packages("phylosmith")
#devtools::install_github('schuyler-smith/phylosmith')
#install.packages(c("RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))

library(lme4)
library(multcompView)
library(readxl)
library(ggplot2)
library(reshape2)  # for melt() function
library(corrplot)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(broom)
library(ggeffects)
library(stargazer)
library(ggfortify)
library(lsmeans)
library(emmeans)
library(car)
library(Hmisc)
library(vegan)
library(ggsci)
library(rstatix)
library(ggpubr)
library(sandwich)
library(msm)
library(nnspat)
library(phyloseq)
library(microbiome)
library(mgcv)
#library(phylosmith)
#library(phyloseq)
#packageVersion("phyloseq")

#######

#FUNGI SPECIES RICHNESS AND STATISTIC ANALYSIS
#---------------------------------------------
#####
setwd(dir = "C://Users//grego/OneDrive/Pollen_Metabarcoding/Alicia Decolle/Decolle_Alicia_2021/R_analysis/data")

dir()

#Import data

data <- read.csv2("fungi_abs_pres2_final.csv", row.names = 1)
t_data<-as.matrix(t(data))

chao_site <- estimateR(t_data)
chao_site


tax <- read.csv2("tax_fungi2.csv",row.names = 1)
samples_df <- read.csv2("samples_fungi2.csv", row.names = 1)

data <- as.matrix(data)
tax <- as.matrix(tax)

str(data)

#Modify the name of site in the data : Canada.x into Canada-x: 
colnames(data) <- c("Canada-2","Canada-3","Canada-4","Canada-5","Canada-6","Canada-7", "Canada-8","Canada-9","Canada-10","Canada-11","Canada-12","Canada-13","Canada-14","Canada-15","Canada-16","Canada-17","Canada-18","Canada-19","Canada-20","Canada-21","Canada-22","Canada-23","Canada-24","Canada-25","Canada-26","Canada-27","Canada-28","Canada-29","Canada-30","Canada-31","Canada-33","Canada-34","Canada-35","Canada-36","Canada-37","Canada-38","Canada-39","Canada-40","Canada-41","Canada-42","Canada-43","Canada-44")

#Transform to phyloseq objects

species = otu_table(data, taxa_are_rows = TRUE)
TAX = tax_table(tax)
samples = sample_data(samples_df)

taxa_names(species)
taxa_names(TAX)
sample_names(species)

# Reorder the months for the graph.
samples$month <- factor(samples$month, levels = c("May", "June", "July", "August", "September"))

#Creation of the carbom
carbom <- phyloseq(species, TAX, samples)

#Visualization of the data characteristics

sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)


####Rarecurve
otu_table <- as(otu_table(carbom), "matrix")
otu_table <- t(otu_table)  # Transpose if necessary for vegan
rare_curve_data <- rarecurve(otu_table, step = 100, cex = 0.5, label = FALSE)




#########################
###Species composition####
#########################
############################FUNGI
###Species pruning
##By site
carbom2 <- aggregate_taxa(carbom, "Species")
carbom2<-merge_samples(carbom2, "month")


#filter <- phyloseq::genefilter_sample(carbom2)
carbom2_filtered = filter_taxa(carbom2, function(x) sum(x) > .005, TRUE)
#carbom2_filtered <- prune_taxa(filter, carbom2)
carbom2_filtered

#FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)
#carbom2 = filter_taxa(carbom2, function(x) x > 1000, TRUE)
#carbom2_filtered = filter_taxa(carbom2_filtered, function(x) sum(x) > .005, TRUE)
trans <- transform_sample_counts(carbom2_filtered, function(x) x / sum(x))

mdf = psmelt(trans)
#plot_heatmap(trans, sample.label="site", low="#66CCFF", high="#000033", na.value="white")
#Plot heatmap
mdf$Sample <- factor(mdf$Sample, levels=c("March", "April", "May","June","July","August","September"))
#Plot heatmap
dk1<-ggplot(filter(mdf, Abundance >= 0.05), 
           aes(x = Sample, y = reorder(OTU, Abundance), fill = Abundance)) +
  geom_tile(color = "gray40") +
  scale_fill_gradient2(low = "white", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,0.65), space = "Lab",
                       name="Proportional occurence")+
  theme_bw(9) +
  ylab("Fungi species") +
  xlab("Months") +
  labs(fill = "Proportional\nabundance") +
  theme(axis.text.y = element_text(face = "italic", size = 9)) +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1))
dk1
ggsave(filename = "Month Heatmap Fungi Figure 3.tiff",dk1, width = 15, height = 10,
       units = "cm",dpi = 500 )


###Species pruning
##By Month
carbom2 <- aggregate_taxa(carbom, "Species")
carbom2<-merge_samples(carbom2, "Map_site")


#filter <- phyloseq::genefilter_sample(carbom2)
carbom2_filtered = filter_taxa(carbom2, function(x) sum(x) > .005, TRUE)
#carbom2_filtered <- prune_taxa(filter, carbom2)
carbom2_filtered

#FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)
#carbom2 = filter_taxa(carbom2, function(x) x > 1000, TRUE)
#carbom2_filtered = filter_taxa(carbom2_filtered, function(x) sum(x) > .005, TRUE)
trans <- transform_sample_counts(carbom2_filtered, function(x) x / sum(x))

mdf = psmelt(trans)
#plot_heatmap(trans, sample.label="site", low="#66CCFF", high="#000033", na.value="white")

#Plot heatmap
dk2<-ggplot(filter(mdf, Abundance >= 0.10), 
           aes(x = Sample, y = reorder(OTU, Abundance), fill = Abundance)) +
  geom_tile(color = "gray40") +
  scale_fill_gradient2(low = "white", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,0.90), space = "Lab",
                       name="Proportional occurence")+
  theme_bw(9) +
  ylab("Fungi species") +
  xlab("Sites") +
  labs(fill = "Proportional\nabundance") +
  theme(axis.text.y = element_text(face = "italic", size = 9)) +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1))
dk2
ggsave(filename = "Site Heatmap Fungi Figure 3.tiff",dk2, width = 15, height = 10,
       units = "cm",dpi = 500 )

########################
#Alpha diversity graphs##
########################
#Import data

data <- read.csv2("fungi_abs_pres2_final.csv", row.names = 1)
tax <- read.csv2("tax_fungi2.csv",row.names = 1)
samples_df <- read.csv2("samples_fungi2.csv", header=T, sep=",", dec =".",row.names = 1)

data <- as.matrix(data)
tax <- as.matrix(tax)

str(data)

#Modify the name of site in the data : Canada.x into Canada-x: 
colnames(data) <- c("Canada-2","Canada-3","Canada-4","Canada-5","Canada-6","Canada-7", "Canada-8","Canada-9","Canada-10","Canada-11","Canada-12","Canada-13","Canada-14","Canada-15","Canada-16","Canada-17","Canada-18","Canada-19","Canada-20","Canada-21","Canada-22","Canada-23","Canada-24","Canada-25","Canada-26","Canada-27","Canada-28","Canada-29","Canada-30","Canada-31","Canada-33","Canada-34","Canada-35","Canada-36","Canada-37","Canada-38","Canada-39","Canada-40","Canada-41","Canada-42","Canada-43","Canada-44")

#Transform to phyloseq objects

species = otu_table(data, taxa_are_rows = TRUE)
TAX = tax_table(tax)
samples = sample_data(samples_df)

taxa_names(species)
taxa_names(TAX)
sample_names(species)

# Reorder the months for the graph.
samples$month <- factor(samples$month, levels = c("May", "June", "July", "August", "September"))

#Creation of the carbom
carbom <- phyloseq(species, TAX, samples)

#Visualization of the data characteristics

sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)


####Diversity_Greg
#richness testing
results = estimate_richness(carbom, measures = 'Shannon')
d = sample_data(carbom)
d

#Non phylogenetic Diversity

hmp.div <- microbiome::alpha(carbom, index = c("observed", "diversity_shannon", "chao1"))

# get the metadata out as seprate object
hmp.meta <- meta(carbom)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)

# Add the rownames to diversity table
hmp.div$sam_name <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")

# check the tables
colnames(div.df)

library(ggpubr)
pd.plot3 <- ggboxplot(div.df,
                      x = "Map_site",
                      y = "observed",
                      fill = "Map_site",
                      palette = "jco",
                      ylab = "Species richness")
pd.plot3 <- pd.plot3 + rotate_x_text()+ annotate(geom="text", x=1.5, y= 20, label="pvalue = 0.51",
                                                 color="black")
#ggsave("Fig 3.svg", plot = pd.plot)
pd.plot3

pd.plot2 <- ggboxplot(div.df,
                      x = "Map_site",
                      y = "diversity_shannon",
                      fill = "Map_site",
                      palette = "jco",
                      ylab = "Shannon index")
pd.plot2 <- pd.plot2 + rotate_x_text()
pd.plot2


write.csv2(hmp.div,"Fungi_Diversity_Measurebysamples_2.csv")

#####

p<-plot_richness(carbom, measures=c("Observed", "Shannon"), x = "urbanised_area", color="month")+   
  geom_smooth(method = lm, formula = y ~ x, colour="black")+
  xlab("Rural-Urban Gradient")+
  labs(color = "Months")
p+theme_bw(base_size = 11)

#Other graphs I didn't run (to run alpha_diversity_graph, you need phylosmith packages)
#plot_richness(carbom, measures=c("Observed", "Simpson", "Shannon"), x = "month", color="urbanised_area")
#alpha_diversity_graph(carbom, index = 'shannon', treatment = 'urbanised_area', subset = NULL, colors = 'default')
#alpha_diversity_graph(carbom, index = 'shannon', treatment = 'month', subset = NULL, colors = 'white')

#####
#STATISTIC ANALYSIS
#------------------
#####

#Import data :

data_fungi <- read.csv("Variables_PCA_fungi.csv", header=T, sep=",", dec =".")

str(data)

# Reorder the months for the graph.
data_fungi$month <- factor(data_fungi$month, levels = c("May", "June", "July", "August", "September"))

# GLM analysis
modglm_fungi <- glm(data_fungi$nbr_species_fungi ~ data_fungi$urbanised_area + data_fungi$month, family = gaussian)
modglm_fungi
summary(modglm_fungi)

plot(modglm_fungi)

#shapiro wilk test normally distributed
shapiro.test(modglm_fungi$residuals)

# Reorder the months for the graph.
data_fungi$month <- factor(data_fungi$month, levels = c("May", "June", "July", "August", "September"))

# Boxplot graph observed species according to the month
par(cex.lab=1.5) # is for y-axis
par(cex.axis=1) 
a <-boxplot(nbr_species_fungi ~ month , data=data_fungi, ylab="Number of fungi species (n)",xlab = "Month" , main="")

###GLM
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())



Div_Kingdom<-read.csv2("Div_Kingdoms.csv")

pollen_plant <- read.csv2("Species_PCoA_plants_2.csv")
str(pollen_plant)










pollen_plant_tomerge<-pollen_plant
colnames(pollen_plant_tomerge)[colnames(pollen_plant_tomerge) == "sample_ID"] <- "Sample"
Div_Kingdom<-merge(Div_Kingdom,pollen_plant_tomerge, by = "Sample",all.x = T)
Div_Kingdom<-Div_Kingdom[,-c(14:91)]


gg_obs<-ggplot(data=Div_Kingdom, aes(y = observed_fungi, x = observed_plant)) +
  geom_point(shape = 19) +
  geom_smooth(method = "glm",se = T, color = "black", method.args = list(family = "poisson"))+
  BioR.theme + xlab ("Plant species richness") + ylab("Fungi species richness")
gg_obs
ggsave("gg_obs.png",gg_obs)

gg_shan<-ggplot(data=Div_Kingdom, aes(y = shannon_fungi, x = shannon_plant)) +
  geom_point(shape = 19) +
  geom_smooth(method = "glm",se = T, color = "black", method.args = list(family = "poisson"))+
  BioR.theme + xlab ("Plant Shannon's index") + ylab("Fungi Shannon's index")
gg_shan
ggsave("gg_shan.png",gg_shan)

gg_chao<-ggplot(data=Div_Kingdom, aes(y = chao1_fungi, x = chao1_plant)) +
  geom_point(shape = 19) +
  geom_smooth(method = "glm",se = T, color = "black", method.args = list(family = "poisson"))+
  BioR.theme + xlab ("Plant Chao1 estimation") + ylab("Fungi Chao1 estimation")
gg_chao
ggsave("gg_chao.png",gg_chao)

library(lme4)
# GLM analysis
modglm_corr_obs <- glmer(observed_fungi ~ observed_plant + (1|site), data = Div_Kingdom, family = poisson)

summary(modglm_corr_obs)

plot(modglm_fungi)

modglm_corr_shan <- lmer(shannon_fungi ~ shannon_plant + (1|site), data = Div_Kingdom)

summary(modglm_corr_shan)
Anova.lmer(modglm_corr_shan)




# GLM analysis alpha
library(glmmTMB)

modglm_plant_obs <- glmmTMB(observed_plant ~ month2 + (1|site), data = Div_Kingdom, family = nbinom2)
#modglm_plant_obs<- glmer(observed_plant ~ month2 + (1|site), data = Div_Kingdom, family = poisson)
summary(modglm_plant_obs)

modglm_fungi_obs <- glmmTMB(observed_fungi ~ month2 + (1|site), data = Div_Kingdom, family = nbinom2)

summary(modglm_fungi_obs)

plot(modglm_fungi)

modglm_plant_shan <- lmer(shannon_plant ~ month2 + (1|site), data = Div_Kingdom)

summary(modglm_plant_shan)
Anova.lmer(modglm_plant_shan)


modglm_fungi_shan <- lmer(shannon_fungi ~ month2 + (1|site), data = Div_Kingdom)

summary(modglm_fungi_shan)
Anova.lmer(modglm_fungi_shan)


#GAMM
modglm_plant_obs_gamm <- gam(observed_plant ~ month2 , data = Div_Kingdom, family = poisson)

summary(modglm_plant_obs_gamm $gam) # For gamm
summary(modglm_plant_obs_gamm )     # For bam

plot(modglm_plant_obs_gamm)

library(DHARMa)
#Petite explication de DHARMa
#1. Simuler des résidus
res_simus = simulateResiduals(modglm_plant_obs, n = 1000)
#2. Plotter tes résidus
plot(res_simus)
plotResiduals(res_simus)
#3.Tester par exemple l'overdispersion
testDispersion(res_simus)
testZeroInflation(res_simus)



gg_plant_month_obs<-ggplot(data=Div_Kingdom, aes(y = observed_plant, x = month2)) +
  geom_point(shape = 19) +
  geom_smooth(method = "glm",se = T, color = "black", method.args = list(family = "nbinom2"))+
  BioR.theme + xlab ("Months") + ylab("Plant species richness")
gg_plant_month_obs
ggsave("gg_plant_month_obs.png",gg_plant_month_obs)

gg_fungi_month_obs<-ggplot(data=Div_Kingdom, aes(y = observed_fungi, x = month2)) +
  geom_point(shape = 19) +
  geom_smooth(method = "glm",se = T, color = "black", method.args = list(family = "nbinom2"))+
  BioR.theme + xlab ("Months") + ylab("Fungi species richness")
gg_fungi_month_obs
ggsave("gg_fungi_month_obs.png",gg_fungi_month_obs)

gg_fungi_month_shannon<-ggplot(data=Div_Kingdom, aes(y = shannon_fungi, x = month2)) +
  geom_point(shape = 19) +
  geom_smooth(method = "lm",se = T, color = "black")+
  BioR.theme + xlab ("Months") + ylab("Fungi Shannon's index")
gg_fungi_month_shannon
ggsave("gg_fungi_month_shannon.png",gg_fungi_month_shannon)

gg_plant_month_shannon<-ggplot(data=Div_Kingdom, aes(y = shannon_plant, x = month2)) +
  geom_point(shape = 19) +
  geom_smooth(method = "lm",se = T, color = "black")+
  BioR.theme + xlab ("Months") + ylab("Plant Shannon's index")
gg_plant_month_shannon
ggsave("gg_plant_month_shannon.png",gg_plant_month_shannon)


library(ggpubr)

dc<-ggarrange(gg_plant_month_obs, gg_fungi_month_obs, gg_plant_month_shannon ,gg_fungi_month_shannon, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
ggsave("FigureGLMM_MOnth.png",dc)

## linear mixed models - reference values from older code

#####

#PLANT SPECIES RICHNESS AND STATISTIC ANALYSIS
#---------------------------------------------
#####

#Import data

data_plants <- read.csv2("plants_abs_pres_2.csv", row.names = 1)
tax_plants <- read.csv2("Tax_plants.csv", row.names = 1)
samples_df_plants <- read.csv2("samples_plants.csv", row.names = 1)

data_plants <- as.matrix(data_plants)
tax_plants <- as.matrix(tax_plants)

#Modify the name of site in the data : Canada.x into Canada-x: 
colnames(data_plants) <- c("Canada-2","Canada-3","Canada-4","Canada-5","Canada-6","Canada-7", "Canada-8","Canada-9","Canada-10","Canada-11","Canada-12","Canada-13","Canada-14","Canada-15","Canada-16","Canada-17","Canada-18","Canada-19","Canada-20","Canada-21","Canada-22","Canada-23","Canada-24","Canada-25","Canada-26","Canada-27","Canada-28","Canada-29","Canada-30","Canada-31","Canada-33","Canada-34","Canada-35","Canada-36","Canada-37","Canada-38","Canada-39","Canada-40","Canada-41","Canada-42","Canada-43","Canada-44")

#Transform to phyloseq objects

species_plants = otu_table(data_plants, taxa_are_rows = TRUE)
TAX_plants = tax_table(tax_plants)
samples_plants = sample_data(samples_df_plants)

taxa_names(species_plants)
taxa_names(TAX_plants)
sample_names(species)

# Reorder the months for the graph.
samples_plants$month <- factor(samples_plants$month, levels = c("May", "June", "July", "August", "September"))

#Creation of carbom_plants
carbom_plants <- phyloseq(species_plants, TAX_plants, samples_plants)





#########################
###Species composition####
#########################
###########################PLANTS
###Species pruning
##By month
carbom2 <- aggregate_taxa(carbom_plants, "Species")
carbom2<-merge_samples(carbom2, "month")


#filter <- phyloseq::genefilter_sample(carbom2)
carbom2_filtered = filter_taxa(carbom2, function(x) sum(x) > .005, TRUE)
#carbom2_filtered <- prune_taxa(filter, carbom2)
carbom2_filtered

#FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)
#carbom2 = filter_taxa(carbom2, function(x) x > 1000, TRUE)
#carbom2_filtered = filter_taxa(carbom2_filtered, function(x) sum(x) > .005, TRUE)
trans <- transform_sample_counts(carbom2_filtered, function(x) x / sum(x))

mdf = psmelt(trans)
#plot_heatmap(trans, sample.label="site", low="#66CCFF", high="#000033", na.value="white")
#Plot heatmap
mdf$Sample <- factor(mdf$Sample, levels=c("March", "April", "May","June","July","August","September"))
#Plot heatmap
dk3<-ggplot(filter(mdf, Abundance >= 0.05), 
           aes(x = Sample, y = reorder(OTU, Abundance), fill = Abundance)) +
  geom_tile(color = "gray40") +
  scale_fill_gradient2(low = "white", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,0.5), space = "Lab",
                       name="Proportional occurence")+
  theme_bw(9) +
  ylab("Plant species") +
  xlab("Month") +
  labs(fill = "Proportional\nabundance") +
  theme(axis.text.y = element_text(face = "italic", size = 9)) +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1))
dk3
ggsave(filename = "Month Heatmap Plant Figure 3.tiff",dk3, width = 15, height = 10,
       units = "cm",dpi = 500 )


###Species pruning
##By Month
carbom2 <- aggregate_taxa(carbom_plants, "Species")
carbom2<-merge_samples(carbom2, "Map_site")


#filter <- phyloseq::genefilter_sample(carbom2)
carbom2_filtered = filter_taxa(carbom2, function(x) sum(x) > .005, TRUE)
#carbom2_filtered <- prune_taxa(filter, carbom2)
carbom2_filtered

#FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)
#carbom2 = filter_taxa(carbom2, function(x) x > 1000, TRUE)
#carbom2_filtered = filter_taxa(carbom2_filtered, function(x) sum(x) > .005, TRUE)
trans <- transform_sample_counts(carbom2_filtered, function(x) x / sum(x))

mdf = psmelt(trans)
#plot_heatmap(trans, sample.label="site", low="#66CCFF", high="#000033", na.value="white")

#Plot heatmap
dk4<-ggplot(filter(mdf, Abundance >= 0.1), 
           aes(x = Sample, y = reorder(OTU, Abundance), fill = Abundance)) +
  geom_tile(color = "gray40") +
  scale_fill_gradient2(low = "white", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,0.85), space = "Lab",
                       name="Proportional occurence")+
  theme_bw(9) +
  ylab("Plant species") +
  xlab("Sites") +
  labs(fill = "Proportional\nabundance") +
  theme(axis.text.y = element_text(face = "italic", size = 8)) +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1))
dk4
ggsave(filename = "Site Heatmap Plant Figure 3.tiff",dk4, width = 15, height = 10,
       units = "cm",dpi = 500 )


#Import data
data_plants <- read.csv2("plants_abs_pres.csv", row.names = 1)
tax_plants <- read.csv2("Tax_plants.csv", row.names = 1)
samples_df_plants <- read.csv2("samples_plants.csv", row.names = 1)
data_plants <- as.matrix(data_plants)
tax_plants <- as.matrix(tax_plants)

#Modify the name of site in the data : Canada.x into Canada-x: 
colnames(data_plants) <- c("Canada-2","Canada-3","Canada-4","Canada-5","Canada-6","Canada-7", "Canada-8","Canada-9","Canada-10","Canada-11","Canada-12","Canada-13","Canada-14","Canada-15","Canada-16","Canada-17","Canada-18","Canada-19","Canada-20","Canada-21","Canada-22","Canada-23","Canada-24","Canada-25","Canada-26","Canada-27","Canada-28","Canada-29","Canada-30","Canada-31","Canada-33","Canada-34","Canada-35","Canada-36","Canada-37","Canada-38","Canada-39","Canada-40","Canada-41","Canada-42","Canada-43","Canada-44")

#Transform to phyloseq objects
species_plants = otu_table(data_plants, taxa_are_rows = TRUE)
TAX_plants = tax_table(tax_plants)
samples_plants = sample_data(samples_df_plants)

taxa_names(species_plants)
taxa_names(TAX_plants)
sample_names(species)

# Reorder the months for the graph.
samples_plants$month <- factor(samples_plants$month, levels = c("May", "June", "July", "August", "September"))

#Creation of carbom_plants
carbom_plants <- phyloseq(species_plants, TAX_plants, samples_plants)

####Diversity_Greg
###richness testing
results = estimate_richness(carbom_plants, measures = 'Shannon')
d = sample_data(carbom)


#Non phylogenetic Diversity

hmp.div <- alpha(carbom_plants, index = "all")

# get the metadata out as seprate object
hmp.meta <- meta(carbom_plants)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)

# Add the rownames to diversity table
hmp.div$sam_name <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")

# check the tables
colnames(div.df)

library(ggpubr)
pd.plot3 <- ggboxplot(div.df,
                      x = "Map_site",
                      y = "observed",
                      fill = "Map_site",
                      palette = "jco",
                      ylab = "Species richness")
pd.plot3 <- pd.plot3 + rotate_x_text()+ annotate(geom="text", x=1.5, y= 20, label="pvalue = 0.51",
                                                 color="black")
#ggsave("Fig 3.svg", plot = pd.plot)
pd.plot3

pd.plot2 <- ggboxplot(div.df,
                      x = "Map_site",
                      y = "diversity_shannon",
                      fill = "Map_site",
                      palette = "jco",
                      ylab = "Shannon index")
pd.plot2 <- pd.plot2 + rotate_x_text()
pd.plot2


write.csv2(hmp.div,"Plants_Diversity_Measurebysamples2.csv")




#Alpha diversity graphs
plot_richness(carbom_plants, measures=c("Observed", "Simpson", "Shannon"), x = "urbanised_area", color="month")+   
  geom_smooth(method = glm, formula = y ~ x, colour="black")+
  xlab("Rural-Urban Gradient")+
  labs(color = "Months")

#Other graphs I didn't run (to run alpha_diversity_graph, you need phylosmith packages)
#plot_richness(carbom_plants, measures=c("Observed", "Simpson", "Shannon"), x = "month",  color="urbanised_area")  
#alpha_diversity_graph(carbom_plants, index = 'shannon', treatment = 'urbanised_area', subset = NULL, colors = 'default')
#alpha_diversity_graph(carbom_plants, index = 'shannon', treatment = 'month', subset = NULL, colors = 'default')

#####
#STATISTIC ANALYSIS
#------------------
#####

#Import data

data_plants <- read.csv("Variables_PCA_plants.csv", header=T, sep=",", dec =".")

str(data)

# Reorder the months for the graph.
data_plants$month <- factor(data_plants$month, levels = c("May", "June", "July", "August", "September"))

# GLM analysis
modglm_plants <- glm(data_plants$nbr_species_plant ~ data_plants$urbanised_area + data_plants$month, family = gaussian)
modglm_plants
summary(modglm_plants)

plot(modglm_plants)

#shapiro wilk test normally distributed
shapiro.test(modglm_plants$residuals)

# Reorder the months for the graph.
data_plants$month <- factor(data_plants$month, levels = c("May", "June", "July", "August", "September"))

# Boxplot graph observed species according to the month
par(cex.lab=1.5) # is for y-axis
par(cex.axis=1) 
a <-boxplot(nbr_species_plant ~ month , data=data_plants, ylab="Number of plant species (n)",xlab = "Month" , main="")

#############
#VENN Diagram



pollen_plant <- read.csv2("Species_PCoA_plants_2.csv")
str(pollen_plant)

library(venn)
library(UpSetR)
library(ggpolypath)


filter <-filter(pollen_plant[-c(2:7)],month==c("May"))
r<-colSums(filter[-c(1)])

filter <-filter(pollen_plant[-c(2:7)],month==c("June"))
s<-colSums(filter[-c(1)])


filter <-filter(pollen_plant[-c(2:7)],month==c("July"))
u<-colSums(filter[-c(1)])


filter <-filter(pollen_plant[-c(2:7)],month==c("August"))
uc<-colSums(filter[-c(1)])


filter <-filter(pollen_plant[-c(2:7)],month==c("September"))
se<-colSums(filter[-c(1)])


m<-as.matrix(data.frame(May=c(r),June=c(s),July=c(u),August=c(uc), September=c(se)))
colSums(m)
m<-ifelse(m>0, 1, 0)
venn(as.data.frame(m),opacity=0.5,plotsize=10,col=c(May='red',June='blue',July='green',August='purple', September='orange'),zcolor=c('#E31A1C','#1F78B4','#33A02C','#CAB2D6','#ffa500'), box=FALSE,ilcs = 1.3,sncs=1.4,par=FALSE)


##########

#2. Beta diversity : (NMDS + Permanova)
#______________________________________

# Packages installation & preparation 
#------------------------------------

#install.packages("RVAideMemoire")
#install.packages("rcompanion")
#install.packages("devtools")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(RVAideMemoire)
library(rcompanion)
library(devtools)         #to run pairwise adonis
library(pairwiseAdonis)   #to run pairwise adonis
library(vegan)
library(ggplot2)
#####

Dir <- paste0("C://Users/Abeille/OneDrive/ISHS2021_Pouilloux/Alicia Decolle/Decolle_Alicia_2021/R_analysis/data")
setwd(Dir)

#FUNGI SPECIES COMPOSITION - NMDS & PERMANOVA (+ post hoc) 
#---------------------------------------------------------

#import data

pollen_fungi <- read.csv("Species_PCoA_fungi.csv")
str(pollen_fungi)

#define the factors (urbanised_area is not define as a factor to conserve the variable as continuous)
pollen_fungi$season<-factor(pollen_fungi$season,levels=c("Spring","Summer","Fall"))
pollen_fungi$month <- factor(pollen_fungi$month, levels=c("May","June","July", "August", "September"))

#jaccard index for presence/absence matrix
pollen.matrix_fungi<-as.matrix(pollen_fungi[,-c(1:7)])
pollen_jac_fungi <- vegdist(pollen.matrix_fungi,method = "jaccard", binary = T)

### NMDS dimension, stress should be the smallest possible

pollen.mds_fungi <- metaMDS(pollen_jac_fungi, distance = "jaccard",k=4, autotransform = FALSE,maxit = 999,trymax = 250)
pollen.spp.fit_fungi <- envfit(pollen.mds_fungi, pollen_fungi [-c(1:7)], permutations = 999)
site.scrs_fungi <- as.data.frame(scores(pollen.mds_fungi, display = "sites"))
spp.scrs_fungi <- as.data.frame(scores(pollen.spp.fit_fungi, display = "vectors"))
spp.scrs_fungi <- cbind(spp.scrs_fungi, Species = rownames(spp.scrs_fungi))
spp.scrs_fungi <- cbind(spp.scrs_fungi, pval = pollen.spp.fit_fungi$vectors$pvals)
sig.spp.scrs_fungi <- subset(spp.scrs_fungi, pval<=0.05)

pollen.mean_fungi=aggregate(pollen.mds_fungi$points[ ,1:2],list(group=pollen_fungi$season),mean)
NMDS_fungi = data.frame(NMDS1 = pollen.mds_fungi$points[,1], NMDS2 = pollen.mds_fungi$points[,2],group=pollen_fungi$season)

###plot to determine the goodness of the fit
goodness(pollen.mds_fungi)
stressplot(pollen.mds_fungi)
pollen.mds_fungi$stress

# properties of the graphs (shapes for months and colors for Urban-Rural gradient)
shapes_fungi = c(21,22,23,24) 
shapes <- shapes_fungi[as.numeric(pollen_fungi$season)]
veganCovEllipse <-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

vector_fungi =pollen_fungi$urbanised_area
vector_fungi2=pollen_fungi$month2

col_vector_fungi <- factor(pollen_fungi$urbanised_area, levels=unique(pollen_fungi$urbanised_area)) #suite pas possible car trop d'?l?ments ? mettre en couleur (max 8 ?l?ments or ici on a 12)
col_vector_fungi
col_vector_fungi2 <- factor(pollen_fungi$season, levels=c("Spring","Summer","Fall"))
col_vector_fungi2


jColors <-
  with(pollen_fungi,
       data.frame(season = levels(season),
                  color = "black"))

(jDatColor <- merge(pollen_fungi, jColors))

#function for the ellipse
ord_fungi<-ordiellipse(pollen.mds_fungi,groups=col_vector_fungi2,conf = 0.8, label=TRUE)

df_ell <- data.frame()
for(g in levels(NMDS_fungi$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_fungi[NMDS_fungi$group==g,],
                                                   veganCovEllipse(ord_fungi[[g]]$cov,ord_fungi[[g]]$center,ord_fungi[[g]]$scale)))
                                ,group=g))
}

#NMDS ordination graph

ggplot(site.scrs_fungi, aes(x=NMDS1, y=NMDS2,colour=pollen_fungi$urbanised_area))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2,shape=pollen_fungi$month,fill=pollen_fungi$urbanised_area),size=3)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14),axis.text=element_text(size=12),axis.title=element_text(size=14))  +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2),col="black", linetype=1)+
  scale_color_continuous(name="Rural-Urban Gradient",type = "viridis", guide = guide_colorbar(nbin = 5))+
  scale_fill_continuous(name="Rural-Urban Gradient",type = "viridis", guide = guide_colorbar(nbin = 5))+
  #scale_color_gradient2(name="Rural-Urban Gradient",low="yellow", mid="red",high="brown", midpoint = mean(pollen_fungi$urbanised_area),guide = guide_colorbar(nbin = 2))+
  scale_shape_manual(name="Month",values=c(23, 22, 25,24,21))+
  annotate("text",x=c(0.38,0.07,0.10),y=c(0.28,-0.36,0.45),label=c("Spring","Summer","Fall"),size=5)+
  xlab(label="NMDS1")+
  ylab(lab="NMDS2") + theme_bw(base_size = 12)

#PERMANOVA to determine the level of significance between the different plant communities depending on the factors

adonis(pollen_fungi[-c(1:7)] ~ month+urbanised_area,
       method = "jaccard",
       data = pollen_fungi,
       permutations = 999, by = "terms")

# Post-Hoc test : Pairwise Adonis:

pad_fungi<- pairwise.adonis(x=pollen_fungi[-c(1:7)],factors=pollen_fungi$month, sim.function='vegdist', sim.method='jaccard',p.adjust.m='holm')
summary(pad_fungi)

library(stringr)
#attribute letter to significant difference
pad_fungi$pairs<-str_replace_all(pad_fungi$pairs, fixed(" vs "), "-")
cldList(p.adjusted ~ pairs,
        data = pad_fungi,
        threshold  = 0.05)

#Discriminating species between two groups using Bray-Curtis dissimilarities 
#r<-simper(pollen_fungi[-c(1:7)], pollen_fungi$month, permutations=100)
#summary(r)

#####

#PLANT SPECIES COMPOSITION - NMDS & PERMANOVA (+ post hoc) 
#---------------------------------------------------------

#import data

pollen_plant <- read.csv2("Species_PCoA_plants_2.csv")
str(pollen_plant)

#define the factors (urbanised_area is not define as a factor to conserve the variable as continuous)
pollen_plant$season<-factor(pollen_plant$season,levels=c("Spring","Summer","Fall"))
pollen_plant$month <- factor(pollen_plant$month, levels=c("May","June","July", "August", "September"))

#jaccard index for presence/absence matrix
pollen.matrix_plant<-as.matrix(pollen_plant[,-c(1:7)])
pollen_bray_plant <- vegdist(pollen.matrix_plant,method = "bray")

####################################
##Multivariate Analysis - RDA Plant##
####################################




library(tidyverse)
library(psych)
library(adespatial)

# Set seed
set.seed(123)


#--------------#
#
# Multicollinearity checks# Check here https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
#
#--------------#

# Remove correlated variables
env.data = subset(pollen_plant, select = -c(month,season))
#pairs.panels(env.data[,c(3,10:16)], scale = TRUE)

envidata<-env.data[,c(1:5)]

fungi_ccoa_hell<-decostand(fungi_ccoa, method = "hellinger")


vegou<-vegdist(pollen.matrix_plant, method = "hellinger")
spe.hell<-decostand(pollen.matrix_plant, method = "hellinger")
spe.hell<-vegdist(spe.hell,method="bray")


pollen_plant$site<-as.factor(pollen_plant$site)
# Perform RDA with all variables
rda1 = dbrda(spe.hell ~ as.factor(pollen_plant$month2), scale = TRUE)
#rda1$concont

# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(rda1)

# Model summaries

vif.cca(rda1) # variance inflation factor (<10 OK)
anova.cca(rda1, permutations = 1000) # full model
anova.cca(rda1, permutations = 1000, by="margin") # per variable 

library(RVAideMemoire)
pairwise.perm.manova(spe.hell, as.factor(pollen_plant$month2), R2 = T, p.method = "holm")

ordiplot(rda1, scaling = 1, type = "text")
ordiplot(rda1, scaling = 2, type = "text")
ordiplot(rda1)
         

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "unconstrained"))
summary(eigenvals(rda1, model = "constrained"))
rda1

screeplot(rda1)


###New plot
# Extract scores
site_scores <- scores(rda1, display = "sites", scaling = 2)  # Site scores
species_scores <- scores(rda1, display = "species", scaling = 2)  # Species scores
biplot_scores <- scores(rda1, display = "bp", scaling = 2)  # Environmental variables

# Convert to data frames for ggplot2
site_scores_df <- as.data.frame(site_scores)
species_scores_df <- as.data.frame(species_scores)
biplot_scores_df <- as.data.frame(biplot_scores)

# Add the factorial variable from dune.env (e.g., "Management")
site_scores_df$Month <- pollen_plant$month
site_scores_df$ID <- pollen_plant$site
site_scores_df$ID<-as.factor(site_scores_df$ID)

#plot
str(centroid.long1)
centroid.long1$axis1c<-as.numeric(centroid.long1$axis1c)
centroid.long1$axis2c<-as.numeric(centroid.long1$axis2c)
centroid.long1$Centroid<-as.character(centroid.long1$axis2c)



p <- ggplot(site_scores_df, aes(x = dbRDA1, y = dbRDA2,color = Month)) +
  geom_point(size = 3)  + 
  labs(title = "", x = "dbRDA1 (15.19%)", y = "dbRDA2 (5.96%)") +
  stat_ellipse(type = "t", size = 1,level = 0.95)

p
#+
#  geom_point(data=centroid.long1, 
#              aes(x=axis1c, y=axis2c), 
#              alpha=0.7, size=5, colour="orange", shape = 18) +
#   geom_text(data = centroid.long1, 
#             aes(x = axis1c, y = axis2c, label = Centroid), vjust = -1, size = 4,
#             alpha = 0.7,  color = "orange", fontface = "bold")

p<-p + theme_bw(base_size = 12)

p

nick<-betadisper(vegou, pollen_plant$month, type = c("centroid"), bias.adjust = FALSE,
           sqrt.dist = FALSE, add = FALSE)


plot(nick)


## Permutation test for F
permutest(nick, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(rda1))
plot(mod.HSD)

ggsave(filename = "RDA_Plant.svg",p,dpi = 1500)

p <- ggplot(site_scores_df, aes(x = dbRDA1, y = MDS1)) +
  geom_point(size = 3)  + 
  geom_text(data = site_scores_df,
            aes(x = dbRDA1, y = MDS1, label = pollen_plant$site), vjust = -1, size = 4)+
  labs(title = "", x = "dbRDA1 (5.49%)", y = "MDS1 (19.32%)") +
  geom_point(data=centroid.long1, 
             aes(x=axis1c, y=axis2c), 
             alpha=0.7, size=5, colour="orange", shape = 18) +
  geom_text(data = centroid.long1, 
            aes(x = axis1c, y = axis2c, label = Centroid), vjust = -1, size = 4,
            alpha = 0.7,  color = "orange", fontface = "bold")
p

p<-p + theme_bw(base_size = 12)

p

ggsave(filename = "Plant_RDA.svg",p,dpi = 1500)

####FUNGI####

Month2<-read.csv2("Month2.csv")
Month2$Month<-as.factor(Month2$Month)

# Perform RDA with all variables
rda2 = dbrda(fungi_ccoa_bray ~ Month2$Month, scale = TRUE)
#rda1$concont

# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(rda2)

# Model summaries

vif.cca(rda2) # variance inflation factor (<10 OK)
anova.cca(rda2, permutations = 1000) # full model
anova.cca(rda2, permutations = 1000, by="margin") # per variable 
pairwise.perm.manova(fungi_ccoa_bray, Month2$Month, R2 = T, p.method = "holm")

ordiplot(rda1, scaling = 1, type = "text")
ordiplot(rda1, scaling = 2, type = "text")
ordiplot(rda1)


# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "unconstrained"))
summary(eigenvals(rda1, model = "constrained"))
rda1

screeplot(rda1)


###New plot
# Extract scores
site_scores <- scores(rda1, display = "sites", scaling = 2)  # Site scores
species_scores <- scores(rda1, display = "species", scaling = 2)  # Species scores
biplot_scores <- scores(rda1, display = "bp", scaling = 2)  # Environmental variables

# Convert to data frames for ggplot2
site_scores_df <- as.data.frame(site_scores)
species_scores_df <- as.data.frame(species_scores)
biplot_scores_df <- as.data.frame(biplot_scores)

# Add the factorial variable from dune.env (e.g., "Management")
site_scores_df$Month <- Month2$Month
site_scores_df$ID <- pollen_plant$site
site_scores_df$ID<-as.factor(site_scores_df$ID)

Month2$Month<-factor(Month2$Month, levels = c("May", "June", "July", "August", "September"))

p <- ggplot(site_scores_df, aes(x = dbRDA1, y = dbRDA2,color = Month2$Month)) +
  geom_point(size = 3)  + 
  labs(title = "", x = "dbRDA1 (5.13%)", y = "dbRDA2 (4.69%)", color = "Month") +
  stat_ellipse(type = "t", size = 1,level = 0.95)
p<-p + theme_bw(base_size = 12)

p
ggsave("RDA_Fungi.svg",p,dpi = 1500)


#+
#  geom_point(data=centroid.long1, 
#              aes(x=axis1c, y=axis2c), 
#              alpha=0.7, size=5, colour="orange", shape = 18) +
#   geom_text(data = centroid.long1, 
#             aes(x = axis1c, y = axis2c, label = Centroid), vjust = -1, size = 4,
#             alpha = 0.7,  color = "orange", fontface = "bold")









# Custom triplot code!

## extract % explained by the first 2 axes
perc <- round(100*(summary(rda1)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(rda1, display="sites", choices=c(1,2), scaling=2)
sc_sp <- scores(rda1, display="species", choices=c(1,2), scaling=2)
sc_bp <- scores(rda1, display="bp", choices=c(1, 2), scaling=2)

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(rda1,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - Plant",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 1)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 0.65, 
     font = 1)


##Plot

# Create colour scheme
# blue=#377EB8, green=#7FC97F, orange=#FDB462, red=#E31A1C
cols = c("#7FC97F","#377EB8","#FDB462")

plot(rda1, type="n", scaling = 3)


# SITES
points(rda1, display="sites", pch=21, scaling=3, cex=1.5, col="black") # sites



# PREDICTORS
text(rda1, display="bp", scaling=3, col="red1", cex=1, lwd=2)
# PREDICTORS
text(rda1, display="species", scaling=3, col="blue1", cex=1, lwd=2)




library(ggplot2)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)
library(concaveman)


##Code from https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html#4_Example_2:_Dune_meadow_data_with_constrained_ordination
plot2 <- ordiplot(rda1, choices=c(1,2))


sites.long2 <- sites.long(plot2, env.data=pollen_plant$site)
head(sites.long2)

species.long2 <- species.long(plot2)
species.long2


axis.long2 <- axis.long(rda1, choices=c(1, 2))
axis.long2

spec.envfit <- envfit(plot2, env=spe.hell, permutations=99)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)


species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long3 <- species.long2[species.long2$r >= 0.40, ]
species.long3$labels <- make.cepnames(species.long3$labels)
species.long3


vectors.envfit <- envfit(plot2, env=pollen_plant$month2)
vectors.long3 <- vectorfit.long(vectors.envfit)
vectors.long3


axis.long3 <- axis.long(rda1, choices=c(1, 2))
axis.long3

sites.long2<-as.data.frame(cbind(sites.long2,pollen_plant$month2))
colnames(sites.long2)[12] <- "Month"


centroid.long1<-centroids.long(sites.long2,pollen_plant$site, FUN = median, centroids.only=T)

##GGPLOT

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())


sites.long2$axis1<-as.numeric(sites.long2$axis1)
sites.long2$axis2<-as.numeric(sites.long2$axis2)

plotgg5 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long3[1, "label"]) +
  ylab(axis.long3[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour=pollen_plant$month, shape=pollen_plant$month), 
             alpha=0.7, size=3) +
  geom_point(data=species.long3, 
             aes(x=axis1, y=axis2), 
             alpha=0.7, size=5)+
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*.5, y=axis2*.5, label=labels),
                  colour="black", size = 5) +
  geom_segment(data=vectors.long3,
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="orange", size=1.7, arrow=arrow(), label = "Month") +
  geom_text_repel(data=subset(vectors.long3, vector %in% c("Month")), 
                  aes(x=axis1*1.3, y=axis2*1.3, label=vector),
                  colour="orange") +
  geom_point(data=centroid.long1, 
             aes(x=axis1c, y=axis2c,labels = Centroid), 
             alpha=0.7, size=5, colour="blue", shape = 19) + BioR.theme 
plotgg5 + BioR.theme 

ggsave(filename = "RDA_Season_A.svg",plotgg5,dpi = 500)


##GGplot2

plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour=Larvae_Data_Multi$Sites, shape=Larvae_Data_Multi$Sites), 
             size=5) +
  geom_point(data=species.long2, 
             aes(x=axis1, y=axis2)) +theme_ipsum() 
plotgg2

#######################################################################################################
######PCOA ANALYSIS
#PCOA FUngi
#######################################################################
#######################################################################################################

fungidata2<-read.csv2("Fungi_CCOA.csv",row.names=1)
fungi_ccoa<-fungidata2
plant_ccoa<-as.data.frame (pollen.matrix_plant)
head<-rownames(fungi_ccoa)
# Set the row names of the data frame to the values in the first column
rownames(plant_ccoa) <- head

#spe<-colnames(fungi_ccoa)
############PCOA

fungi_ccoa_hell<-decostand(fungi_ccoa, method = "hellinger")
fungi_ccoa_bray<-vegdist(fungi_ccoa,method="bray")

fungi.pcoa <- cmdscale(fungi_ccoa_bray, eig = TRUE) #Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis


eigendata<-fungi.pcoa $eig 
eigendata
valeurs_propres = eigendata/sum(eigendata)*100 
valeurs_propres



##Code from https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html#4_Example_2:_Dune_meadow_data_with_constrained_ordination
plot2 <- ordiplot(scores(fungi.pcoa, choices=c(1,2)))



ordiplot(scores(fungi.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add weighted average projection of species
#spe.wa <- wascores(fungi.pcoa$points[, 1:2], fungi_ccoa)
#text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")



sites.long2 <- sites.long(plot2)
head(sites.long2)


# Change the name of column 'B' to 'NewColumnName'
colnames(pollen_plant)[colnames(pollen_plant) == "sample_ID"] <- "labels"

sites.long2<-merge(sites.long2,pollen_plant, by = "labels", all.x = TRUE)
head(sites.long2)
sites.long2<-sites.long2[,-c(11:87)]


##GGPLOT

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())


sites.long2$axis1<-as.numeric(sites.long2$axis1)
sites.long2$axis2<-as.numeric(sites.long2$axis2)
colnames(sites.long2)[colnames(sites.long2) == "month"] <- "Month"

plotgg5 <- ggplot(data=sites.long2, 
                  aes(x=axis1, y=axis2, colour = Month, fill = Month)) + 
  stat_ellipse(level = 0.8)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("PCoA1 (22.58%)") +
  ylab("PCoA2 (15.84%)") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point( size=3)  +
 BioR.theme 
plotgg5

ggsave(filename = "PCOA_FUNGI.jpeg",plotgg5,dpi = 500)

#PERMANOVA to determine the level of significance between the different plant communities depending on the factors
labels<-row.names(fungi_ccoa)
fungi_ccoa_month<-cbind(labels, fungi_ccoa)
fungi_ccoa_month<-merge(fungi_ccoa_month, pollen_plant[,c(1:5)],by="labels",all.x = T)




perm<-adonis2(fungi_ccoa_hell ~ month,
       method = "bray",
       data = fungi_ccoa_month,
       permutations = 999, by = "terms")

perm

# Post-Hoc test : Pairwise Adonis
pad_fungi<- pairwise.adonis(fungi_ccoa_month[,c(2:46)],factors=pollen_plant$month, sim.function='vegdist', sim.method='jaccard',p.adjust.m='holm')
summary(pad_fungi)

#attribute letter to significant difference
pad_fungi$pairs<-str_replace_all(pad_fungi$pairs, fixed(" vs "), "-")
cldList(p.adjusted ~ pairs,
        data = pad_fungi,
        threshold  = 0.05)

#Discriminating species between two groups using Bray-Curtis dissimilarities
#r<-simper(pollen_plant[-c(1:7)], pollen_plant$month, permutations=100)
#summary(r)



################################################
############PCOA Plant###########################
#################################################

plant_ccoa_hell<-decostand(plant_ccoa, method = "hellinger")
plant_ccoa_bray<-vegdist(plant_ccoa,method="bray")

plant.pcoa <- cmdscale(plant_ccoa_bray, eig = TRUE) #Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis


eigendata<-plant.pcoa $eig 
eigendata
valeurs_propres = eigendata/sum(eigendata)*100 
valeurs_propres



##Code from https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html#4_Example_2:_Dune_meadow_data_with_constrained_ordination
plot2 <- ordiplot(scores(plant.pcoa, choices=c(1,2)))



ordiplot(scores(plant.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add weighted average projection of species
#spe.wa <- wascores(plant.pcoa$points[, 1:2], plant_ccoa)
#text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")



sites.long2 <- sites.long(plot2)
head(sites.long2)


# Change the name of column 'B' to 'NewColumnName'
colnames(pollen_plant)[colnames(pollen_plant) == "sample_ID"] <- "labels"

sites.long2<-merge(sites.long2,pollen_plant, by = "labels", all.x = TRUE)
sites.long2<-sites.long2[,-c(10:87)]

sites.long2$axis1<-as.numeric(sites.long2$axis1)
sites.long2$axis2<-as.numeric(sites.long2$axis2)
colnames(sites.long2)[colnames(sites.long2) == "month"] <- "Month"

plotgg5 <- ggplot(data=sites.long2, 
                  aes(x=axis1, y=axis2, colour = Month, fill = Month)) + 
  stat_ellipse(level = 0.8)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("PCoA1 (16.06%)") +
  ylab("PCoA2 (9.99%)") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point( size=3)  +
  BioR.theme 
plotgg5

ggsave(filename = "PCOA_Plant.jpeg",plotgg5,dpi = 500)

#PERMANOVA to determine the level of significance between the different plant communities depending on the factors
labels<-row.names(plant_ccoa)
plant_ccoa_month<-cbind(labels, plant_ccoa)
plant_ccoa_month<-merge(plant_ccoa_month, pollen_plant[,c(1:5)],by="labels",all.x = T)

perm<-adonis2(plant_ccoa_hell ~ month,
             method = "bray",
             data = plant_ccoa_month,
             permutations = 999, by = "terms")

perm

# Post-Hoc test : Pairwise Adonis
pad_plant<- pairwise.adonis(plant_ccoa_month[,c(2:79)],factors=pollen_plant$month, sim.function='vegdist', sim.method='jaccard',p.adjust.m='holm')
summary(pad_plant)

#attribute letter to significant difference
pad_plant$pairs<-str_replace_all(pad_plant$pairs, fixed(" vs "), "-")
cldList(p.adjusted ~ pairs,
        data = pad_plant,
        threshold  = 0.05)

#Discriminating species between two groups using Bray-Curtis dissimilarities
#r<-simper(pollen_plant[-c(1:7)], pollen_plant$month, permutations=100)
#summary(r)






