---
  output: html_document
editor_options: 
  chunk_output_type: console
---
  **Data Manipulation of greenhouse and climate data**
  ===
  
  *Load libraries* IMPORTANT --> DO NOT SKIP
---
  ```{r, include=F}
library(dplyr)
library(tidyr)
library(plyr)
library(vegan)
library(usdm)
library(ggplot2)
#install.packages("multcompView")
library(multcompView)
#install.packages("lsmeans")
library(emmeans) 

#install.packages("ggpubr")
library(ggpubr)
library(ggfortify)
library(nlme)
library(lmerTest)
library(lme4)
#library(mvpart) # for MRT, package needs to be installed from github using devtools::install_github("cran/mvpart")
library(visreg)
library(gridExtra)
library(multcomp)
library(phia)
library(Hmisc) # for rcorr
#install.packages("corrplot")
library(corrplot)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

```

*Load data* IMPORTANT --> DO NOT SKIP

```{r}
#load individual common garden data

phenraw <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/Garden.data.compiled.long.txt", header = T)

#set factors
phenraw$Treatment <- factor(phenraw$Treatment)
phenraw$Pop.ID <- factor(phenraw$Pop.ID)
phenraw$Mom <- factor(phenraw$Mom)
phenraw$Range  <- factor(phenraw$Range)
phenraw$gender <- factor(phenraw$gender)
phenraw$t4.stem <- as.numeric(phenraw$t4.stem) 
phenraw$group <- as.factor(paste(phenraw$Treatment,phenraw$Range, sep="_") )
phenraw$Ind.ID <-as.factor(paste(phenraw$Treatment, phenraw$Mom, sep="_"))
#drop shoot1 as not variable and highly skewed
phenraw <- phenraw[ ,-which(names(phenraw) %in% c("t0.shoots"))]


```

*Aggregate by Treatment and Population ID to get averages* IMPORTANT --> DO NOT SKIP

---
  ```{r}
popID<-phenraw[,2:3]
popID<-unique(popID)
phen<-cbind(phenraw[,1:2], phenraw[,6:dim(phenraw)[2]])
phen$Range<-ifelse(phenraw$Range=="native", 1,0)
colnames(phen)[1:2] <- c("Treatment", "Pop.ID")
agg_phen_mean <- aggregate(phen, by= list(phen$Treatment, phen$Pop.ID), FUN=mean, na.rm=T)
colnames(agg_phen_mean)[1:2] <- c("Treatment", "Pop.ID")
agg_phen_mean <- agg_phen_mean[,-c(3:4)] # remove empty columns
agg_phen_mean <- agg_phen_mean[,-c(71:72)] # remove empty columns

agg_phen_mean<- agg_phen_mean[ ,-which(names(agg_phen_mean) %in% c("t4t3.leaves"))]
head(agg_phen_mean[,28:71])


#get corresponding matrix of bioclim variables and their PCA

*Calculate PC scores for environmental variables*
  
  ```{r}
#with standarisation --> better solution for climate, without latitude
bioclim<-read.table("~/Dropbox/Documents/Canada_Thistle_CG/Bioclim.txt", header = T)

bioclim<-bioclim[bioclim$Pop.ID != "150908-1",]
bioclim<-merge(bioclim, popID, by="Pop.ID")
climnolatnolonst_pca <- prcomp((bioclim[,5:23]), scale=T)
climlatst_pca <- prcomp((bioclim[,c(3:23)]), scale=T)

autoplot(climnolatnolonst_pca,	type="obs")
pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_nolatlon.pdf", height=5, width=8)
autoplot(climnolatnolonst_pca, data=bioclim, colour="Range",label.label=bioclim$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#autoplot(climnolatnolonst_pca, data=bioclim, colour="Range",label.label="Pop",loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
dev.off()
pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges.pdf", height=5, width=8)
autoplot(climlatst_pca, data=bioclim, colour="Range",label.label="Pop",loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#autoplot(climlatst_pca, data=bioclim, colour="Range",label.label=bioclim$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
dev.off()
aload <- abs(climlatst_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(climnolatnolonst_pca)
bioclim$CLIMPC1 <- climnolatnolonst_pca$x[,1]
bioclim$CLIMPC2 <- climnolatnolonst_pca$x[,2]
bioclim$CLIMPC3 <- climnolatnolonst_pca$x[,3]

#bioclim$CLIMPC1 <- climlatst_pca$x[,1]
#bioclim$CLIMPC2 <- climlatst_pca$x[,2]
#bioclim$CLIMPC3 <- climlatst_pca$x[,3]


##plot correlation matrix of bioclim variables
bioclim2<- bioclim[ ,-which(names(bioclim) %in% c("Pop.acronym"))]
corrs<-rcorr(as.matrix(bioclim2[,3:dim(bioclim2)[2]]), type="spearman")
corrs_r <- as.matrix(corrs$r)

#pdf("~/Dropbox/Documents/Canada_Thistle_CG/BIO_cor.pdf", height=5, width=8)
corrplot(corrs_r)
#dev.off()



#redo with all treatments but shading

#redo with all treatments
#reduce the trait number by removing traits with NA #remove herbivory leaf length   
agg_phen_mean_all <- agg_phen_mean[ agg_phen_mean$Treatment != "Shading",-which(names(agg_phen_mean) %in% c("flowering_day", "Shading.death", "flowers", "Ind.ID","t1.LLL", "t0.leaves", "t1.LLW", "t1.ellipse","t0.ellipse", "t0.LLL", "t0.LLW", "t2.LLL", "t2.LLW", "t2.ellipse", "t3.LLL", "t3.LLW", "t3.ellipse", "t4t3.leaves" , "Above",       "t4.height"))]
#agg_phen_mean_all <- agg_phen_mean[ agg_phen_mean$Treatment=="Nutrient",-which(names(agg_phen_mean) %in% c("flowering_day", "Shading.death", "flowers", "Ind.ID","t1.LLL", "t0.leaves", "t1.LLW", "t1.ellipse","t0.ellipse", "t0.LLL", "t0.LLW", "t2.LLL", "t2.LLW", "t2.ellipse", "t3.LLL", "t3.LLW", "t3.ellipse", "t4t3.leaves", "t4t2.leaves", "t4.leaves", "t4.stem", "t4.height","t4.shoots","t4t0.leaves", "t4t1.leaves"))]
#agg_phen_mean_all <- agg_phen_mean[ ,-which(names(agg_phen_mean) %in% c("flowering_day", "Shading.death", "flowers", "Ind.ID"))]
agg_phen_mean_all<-merge(agg_phen_mean_all, popID, by="Pop.ID")
pc_mean_sub_all <- agg_phen_mean_all[,30:dim(agg_phen_mean_all)[2]-2]
corrs <- rcorr(as.matrix(pc_mean_sub_all), type="spearman")
corrs_r <- data.frame(corrs$r)
corrs_p <- data.frame(corrs$P)
corrs_pbonf <- data.frame(corrs$p_bonf)
corr <- flattenCorrMatrix(corrs_r,corrs_p)
corr$p_bonf <- p.adjust(corr$p, method= "bonferroni")
corr[corr$cor > 0.70,]

#remove correlated variables
variable_cor <- corrs_r
variable_cor[upper.tri(corrs_r)] <- 0
diag(variable_cor) <- 0
variable.new <- na.omit(colnames(pc_mean_sub_all)[!apply(variable_cor,2,function(x) any(abs(x) > 0.7))])
agg_phen_mean_sub_all <- agg_phen_mean_all[ ,which(names(agg_phen_mean_all) %in% variable.new)]
#Control "punch"       "t1t0.leaves" "t1.height"   "t3.shoots"   "t3t0.leaves" "t3.height"   "Above"       "Below"
#Herbivory "punch"       "t3.shoots"   "t3t0.leaves" "t3.height"   "Above"       "Below"
#Shading  "bolting_day" "punch"       "t2.shoots"   "t3.shoots"   "t3t2.leaves" "t3.height"   "Below"
#Nutrient "bolting_day" "punch"       "t3.shoots"   "t3t0.leaves" "t3.height"   "Above"       "Below"      
#agg_phen_mean_sub_all <- agg_phen_mean_all[ ,which(names(agg_phen_mean_all) %in% c("punch", "t3.height", "Below"))]
#agg_phen_mean_sub_all <- agg_phen_mean_all[ ,which(names(agg_phen_mean_all) %in% c("punch", "t1t0.leaves", "t1.height","t3.shoots", "t3t0.leaves", "t3.height", "Above", "Below"))]

#check to make sure they are gone
corrs <- rcorr(as.matrix(agg_phen_mean_sub_all), type="spearman")
corrs_r <- data.frame(corrs$r)
corrs_p <- data.frame(corrs$P)
corrs_pbonf <- data.frame(corrs$p_bonf)
corr <- flattenCorrMatrix(corrs_r,corrs_p)
corr$p_bonf <- p.adjust(corr$p, method= "bonferroni")
corr[corr$cor > 0.7,]

agg_phen_mean_all2<-na.omit(cbind(agg_phen_mean_all$Range, agg_phen_mean_all[, 1:2],  agg_phen_mean_all[,dim(agg_phen_mean_all)[2]], agg_phen_mean_all[, 30:dim(agg_phen_mean_all)[2]-2]))
colnames(agg_phen_mean_all2)[1] <- c("Range")
agg_phen_mean_all2$Range<-ifelse(agg_phen_mean_all2$Range==1,"native", "invasive" )
colnames(agg_phen_mean_all2)[4] <- c("Pop.acronym")
agg_phen_mean_all2$group<-paste(agg_phen_mean_all2$Treatment, agg_phen_mean_all2$Range, sep="_")
pc_mean_sub<-agg_phen_mean_all2[,6:dim(agg_phen_mean_all2)[2]-1]
trait_pca <- prcomp((pc_mean_sub), scale=TRUE)
agg_phen_mean_all2$Range<-as.factor(agg_phen_mean_all2$Range)
autoplot(trait_pca,	type="obs")
pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_trait_all_NO_SHADE.pdf", height=5, width=8)
autoplot(trait_pca, data=agg_phen_mean_all2, colour="group",label.label=agg_phen_mean_all2$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
dev.off()
aload <- abs(trait_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(trait_pca)
agg_phen_mean_all$TRAITPC1 <- trait_pca$x[,1]
agg_phen_mean_all$TRAITPC2 <- trait_pca$x[,2]
agg_phen_mean_all$TRAITPC3 <- trait_pca$x[,3]

pc_mean_sub$TRAITPC1 <- trait_pca$x[,1]
pc_mean_sub$TRAITPC2 <- trait_pca$x[,2]
pc_mean_sub$TRAITPC3 <- trait_pca$x[,3]

corrs<-rcorr(as.matrix(pc_mean_sub, type="spearman"))
corrs_r <- as.matrix(corrs$r)


pdf("~/Dropbox/Documents/Canada_Thistle_CG/TRAIT_all_cor_NO_SHADE.pdf", height=5, width=8)
corrplot(corrs_r)
dev.off()

#PCA reduced corrlation traits
agg_phen_mean_sub<-na.omit(cbind(agg_phen_mean_all$Range, agg_phen_mean_all$Pop.acronym, agg_phen_mean_all[, 1:2], agg_phen_mean_all[,which(names(agg_phen_mean_all) %in% variable.new)]))
colnames(agg_phen_mean_sub)[1] <- c("Range")
colnames(agg_phen_mean_sub)[2] <- c("Pop.acronym")
agg_phen_mean_sub$Range<-ifelse(agg_phen_mean_sub$Range==1,"native", "invasive" )
agg_phen_mean_sub$group<-paste(agg_phen_mean_sub$Treatment,agg_phen_mean_sub$Range, sep="_")
pc_mean_sub<-agg_phen_mean_sub[,6:dim(agg_phen_mean_sub)[2]-1]
trait_pca <- prcomp((pc_mean_sub), scale=TRUE)
agg_phen_mean_sub$Range<-as.factor(agg_phen_mean_sub$Range)
autoplot(trait_pca,	type="obs")
pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_trait_sub_NO_SHADE.pdf", height=5, width=8)
autoplot(trait_pca, data=agg_phen_mean_sub, colour="group",label.label=agg_phen_mean_all2$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
dev.off()
aload <- abs(trait_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(trait_pca)
agg_phen_mean_sub$TRAITPC1 <- trait_pca$x[,1]
agg_phen_mean_sub$TRAITPC2 <- trait_pca$x[,2]
agg_phen_mean_sub$TRAITPC3 <- trait_pca$x[,3]

pc_mean_sub$TRAITPC1 <- trait_pca$x[,1]
pc_mean_sub$TRAITPC2 <- trait_pca$x[,2]
pc_mean_sub$TRAITPC3 <- trait_pca$x[,3]

corrs<-rcorr(as.matrix(pc_mean_sub, type="spearman"))
corrs_r <- as.matrix(corrs$r)
pdf("~/Dropbox/Documents/Canada_Thistle_CG/TRAIT_sub_co_NO_SHADE.pdf", height=5, width=8)
corrplot(corrs_r)
dev.off()


*Test the trait differences between Ranges in MANOVA*
  ---
  
  Do this without any rho>0.70 traits
```{r}

#make sure remaining variables are normal "punch"       "t4.shoots"   "t4t2.leaves" "t4.leaves"   "t4.stem"     "Below"
# punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below
#pdf("~/Dropbox/Documents/Canada_Thistle_CG/hist.pdf", height=5, width=8)
hist(agg_phen_mean_sub$punch[agg_phen_mean_sub$Treatment=="Control"])
hist(agg_phen_mean_sub$punch[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(agg_phen_mean_sub$punch[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(agg_phen_mean_sub$t4.shoots[agg_phen_mean_sub$Treatment=="Control"])
hist(agg_phen_mean_sub$t4.shoots[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(agg_phen_mean_sub$t4.shoots[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(sqrt(agg_phen_mean_sub$t4.shoots)[agg_phen_mean_sub$Treatment=="Control"])
hist(sqrt(agg_phen_mean_sub$t4.shoots)[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(sqrt(agg_phen_mean_sub$t4.shoots)[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(sqrt(agg_phen_mean_sub$t4t2.leaves)[agg_phen_mean_sub$Treatment=="Control"])
hist(sqrt(agg_phen_mean_sub$t4t2.leaves)[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(sqrt(agg_phen_mean_sub$t4t2.leaves)[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(agg_phen_mean_all$t4t2.leaves[agg_phen_mean_all$Treatment=="Control"])
hist(agg_phen_mean_all$t4t2.leaves[agg_phen_mean_all$Treatment=="Herbivory"])
hist((agg_phen_mean_all$t4t2.leaves)[agg_phen_mean_all$Treatment=="Nutrient"])

hist(agg_phen_mean_sub$t4.leaves[agg_phen_mean_sub$Treatment=="Control"])
hist(agg_phen_mean_sub$t4.leaves[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(agg_phen_mean_sub$t4.leaves[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(sqrt(agg_phen_mean_sub$t4.leaves)[agg_phen_mean_sub$Treatment=="Control"])
hist(sqrt(agg_phen_mean_sub$t4.leaves)[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(sqrt(agg_phen_mean_sub$t4.leaves)[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(agg_phen_mean_sub$t4.stem[agg_phen_mean_sub$Treatment=="Control"])
hist(agg_phen_mean_sub$t4.stem[agg_phen_mean_sub$Treatment=="Herbivory"])
hist((agg_phen_mean_sub$t4.stem)[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(agg_phen_mean_sub$Below[agg_phen_mean_sub$Treatment=="Control"])
hist(agg_phen_mean_sub$Below[agg_phen_mean_sub$Treatment=="Herbivory"])
hist((agg_phen_mean_sub$Below)[agg_phen_mean_sub$Treatment=="Nutrient"])

hist(log(agg_phen_mean_sub$Below)[agg_phen_mean_sub$Treatment=="Control"])
hist(log(agg_phen_mean_sub$Below)[agg_phen_mean_sub$Treatment=="Herbivory"])
hist(log(agg_phen_mean_sub$Below)[agg_phen_mean_sub$Treatment=="Nutrient"])


pc_mean_sub<-na.omit(agg_phen_mean_sub[,5:(dim(agg_phen_mean_sub)[2]-4)])
agg_phen_mean_sub2<-na.omit(cbind(agg_phen_mean_sub[,1:3], agg_phen_mean_sub[,4:(dim(agg_phen_mean_sub)[2]-4)]))
agg_phen_mean_sub2$group<-paste(agg_phen_mean_sub2$Treatment, agg_phen_mean_sub2$Range, sep="_")
agg_phen_mean_sub2$group<-as.factor(agg_phen_mean_sub2$group)
trait_pca <- prcomp((pc_mean_sub), scale=TRUE)
agg_phen_mean_sub2$Range<-as.factor(agg_phen_mean_sub2$Range)
autoplot(trait_pca,	type="obs")
#pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_all_trait.pdf", height=5, width=8)
autoplot(trait_pca, data=agg_phen_mean_sub2, colour="group",label.label=agg_phen_mean_sub2$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#dev.off()
aload <- abs(trait_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(trait_pca)
agg_phen_mean_sub$TRAITPC1 <- trait_pca$x[,1]
agg_phen_mean_sub$TRAITPC2 <- trait_pca$x[,2]
agg_phen_mean_sub$TRAITPC3 <- trait_pca$x[,3]

##ALL_mannova<-merge(bioclim, agg_phen_mean_sub2, by="Pop.ID",keep=T)
ALL_mannova<-merge(bioclim, agg_phen_mean_all, by="Pop.ID",keep=T)

"punch"       "t4.shoots"   "t4t2.leaves" "t4.leaves"   "t4.stem"     "Below"     

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x*Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Treatment*Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Treatment*Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Treatment+Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x+Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")

#univariate analyis

#univariate analysis, pop means Range, Treatment
#punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below 
#punch not sig range or interaction
Anova(lm(punch~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(punch~ Treatment+Range.x , data = ALL_mannova)
testInteractions(m, pairwise = "Treatment") 

#t4.shoots not sig range or interaction
Anova(lm(sqrt(t4.shoots)~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(sqrt(t4.shoots)~ Treatment+Range.x , data = ALL_mannova)
testInteractions(m, pairwise = "Treatment") 

#t4t2.leaves  - range sig
Anova(lm(sqrt(t4t2.leaves)~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(t4t2.leaves~ Treatment*Range.x , data = ALL_mannova)
testInteractions(m, pairwise = "Treatment") 


marginal=lsmeans(m,~Range.x:Treatment, type = "response")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 
#plot lsmeans and se
A<-ggplot(CLD,
          aes(x     = Treatment,
              y     = lsmean,
              color= Range.x,
              label = .group)) +
  
  geom_point(shape  = 15,
             size   = 4,
             position = pd, show.legend = FALSE) +
  
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7,
                position = pd, show.legend = FALSE) +
  
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("change in leaf number") +
  xlab("treatment")

#t4.leaves range significant
Anova(lm(sqrt(t4.leaves)~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(t4.leaves~ Treatment*Range.x , data = ALL_mannova)
marginal=lsmeans(m,~Range.x:Treatment, type = "response")

testInteractions(m, pairwise = "Treatment") 

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 
#plot lsmeans and se
B<-ggplot(CLD,
          aes(x     = Treatment,
              y     = lsmean,
              color= Range.x,
              label = .group)) +
  
  geom_point(shape  = 15,
             size   = 4,
             position = pd, show.legend = FALSE) +
  
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7,
                position = pd, show.legend = FALSE) +
  
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("leaf number") +
  xlab("treatment") 

#t4.stem
Anova(lm(t4.stem~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(t4.stem~ Treatment*Range.x , data = ALL_mannova)
testInteractions(m, pairwise = "Treatment") 
marginal=lsmeans(m,~Range.x:Treatment, type = "response")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

C<-ggplot(CLD,
          aes(x     = Treatment,
              y     = lsmean,
              color= Range.x,
              label = .group)) +
  
  geom_point(shape  = 15,
             size   = 4,
             position = pd, show.legend = FALSE) +
  
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7,
                position = pd, show.legend = FALSE) +
  
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("stem diameter") +
  xlab("treatment") 


#Below range and interaction not significant
Anova(lm(Below~ Treatment*Range.x , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x , data = ALL_mannova), type="III")
m<-lm(Below~ Treatment*Range.x , data = ALL_mannova)
testInteractions(m, pairwise = "Treatment") 
marginal=lsmeans(m,~Range.x:Treatment, type = "response")

pdf("~/Dropbox/Documents/Canada_Thistle_CG/NoShade_trait_lat.pdf", height=5, width=8)
ggarrange(A, B, C,   
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()
#PC1 PC2 latitude shoot0

#pop means Range, Treatment, CLIMPC1 MANOVA - Range and CLIMPC1 highly correlated
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC1*Range.x*Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x*CLIMPC1*Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC1*Treatment*Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC1+Range.x+Treatment+Range.x:Treatment+Range.x:CLIMPC1+Treatment:CLIMPC1, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x+CLIMPC1+Treatment+Range.x:Treatment+Range.x:CLIMPC1+Treatment:CLIMPC1, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC1+Treatment+Range.x++Range.x:Treatment+Range.x:CLIMPC1+Treatment:CLIMPC1, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC1+Treatment+Range.x++Range.x:Treatment+Range.x:CLIMPC1, data=ALL_mannova)
summary(m.man, test="Wilks")

#neither climate of range significant nor interaction
Anova(lm(punch~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC1 , data = ALL_mannova), type="III")

#t4.shoots - only treatment significant
Anova(lm(sqrt(t4.shoots)~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC1 , data = ALL_mannova), type="III")


#t4t2.leaves - range pc1 interaction
Anova(lm(sqrt(t4t2.leaves)~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = ALL_mannova), type="III")

m<-lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1, data = ALL_mannova)
interactionMeans(m, factor="Range.x", slope = "CLIMPC1")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC1")
testInteractions(m)

testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(ALL_mannova$CLIMPC1[which(ALL_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(ALL_mannova$CLIMPC1[which(ALL_mannova$Range.x=="native")]))) 



Native<- ALL_mannova[ALL_mannova$Range.x=='native',]
Anova(lm(sqrt(t4t2.leaves)~ Treatment+CLIMPC1 , data = Native), type="III")
summary(lm(sqrt(t4t2.leaves)~ Treatment+CLIMPC1 , data = Native))


Invasive<- ALL_mannova[ALL_mannova$Range.x=='invasive',]
Anova(lm(sqrt(t4t2.leaves)~ Treatment+CLIMPC1 , data = Invasive), type="III")
summary(lm(sqrt(t4t2.leaves)~ Treatment+CLIMPC1 , data = Invasive))



#t4 leaves - range pc1 interaction
Anova(lm(sqrt(t4.leaves)~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC1+ Treatment:Range.x + Treatment:CLIMPC1 + Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = ALL_mannova), type="III")
m<-lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC1+ Range.x:CLIMPC1, data = ALL_mannova)
interactionMeans(m, factor="Range.x", slope = "CLIMPC1")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC1")
testInteractions(m)

Native<- ALL_mannova[ALL_mannova$Range.x=='native',]
Anova(lm(sqrt(t4.leaves)~ Treatment+CLIMPC1 , data = Native), type="III")
summary(lm(sqrt(t4.leaves)~ Treatment+CLIMPC1 , data = Native))
#test if most northern Australian plants differ significantly from the native range in size
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(ALL_mannova$CLIMPC1[which(ALL_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(ALL_mannova$CLIMPC1[which(ALL_mannova$Range.x=="native")]))) 




Invasive<- ALL_mannova[ALL_mannova$Range.x=='invasive',]
Anova(lm(sqrt(t4.leaves)~ Treatment+CLIMPC1 , data = Invasive), type="III")
summary(lm(sqrt(t4.leaves)~ Treatment+CLIMPC1 , data = Invasive))



#t4 stem  - NO effect of range or PC1 - cancel out
Anova(lm(t4.stem~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC1, data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x, data = ALL_mannova), type="III")

#Below  - NO effect of range or pc1 with biomass
Anova(lm(Below~ Treatment*Range.x*CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Treatment:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC1+Treatment:Range.x+Range.x:CLIMPC1 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC1, data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x, data = ALL_mannova), type="III")






#pop means Range, Treatment, CLIMPC2 MANOVA - Range and CLIMPC1 highly correlated
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC2*Range.x*Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x*CLIMPC2*Treatment, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC2*Treatment*Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC2+Range.x+Treatment+Range.x:Treatment+Range.x:CLIMPC2+Treatment:CLIMPC2, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x+CLIMPC2+Treatment+Range.x:Treatment+Range.x:CLIMPC2+Treatment:CLIMPC2, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC2+Treatment+Range.x++Range.x:Treatment+Range.x:CLIMPC2+Treatment:CLIMPC2, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ CLIMPC2+Treatment+Range.x+Range.x:Treatment+Range.x:CLIMPC2, data=ALL_mannova)
summary(m.man, test="Wilks")

testInteractions(m.man, pairwise = "Range.x", test="Wilks")
testInteractions(m.man, pairwise = "Range.x",slope="CLIMPC2", test="Wilks")


#univariate tests

#punch range*pc2 significant
Anova(lm(punch~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC2+ Treatment:Range.x + Treatment:CLIMPC2 + Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova), type="III")
m<-lm(punch~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova)
interactionMeans(m, factor="Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m)

testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="native")]))) 

Native<- ALL_mannova[ALL_mannova$Range.x=='native',]
Anova(lm(punch~ CLIMPC2 , data = Native), type="III")
summary(lm(punch~ CLIMPC2 , data = Native))

Invasive<- ALL_mannova[ALL_mannova$Range.x=='invasive',]
Anova(lm(punch~ CLIMPC2 , data = Invasive), type="III")
summary(lm(punch~ CLIMPC2 , data = Invasive))


#t4.shoots range and pc2 not significant
Anova(lm(sqrt(t4.shoots)~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC2+ Treatment:Range.x + Treatment:CLIMPC2 + Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x+CLIMPC2, data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment+Range.x, data = ALL_mannova), type="III")

#t4t2.leaves - range pc2 interaction
Anova(lm(sqrt(t4t2.leaves)~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC2+ Treatment:Range.x + Treatment:CLIMPC2 + Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova), type="III")

m<-lm(sqrt(t4t2.leaves)~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova)
interactionMeans(m, factor="Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m)

testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="native")]))) 

Native<- ALL_mannova[ALL_mannova$Range.x=='native',]
Anova(lm(sqrt(t4t2.leaves)~ CLIMPC2 , data = Native), type="III")
summary(lm(sqrt(t4t2.leaves)~ CLIMPC2 , data = Native))

Invasive<- ALL_mannova[ALL_mannova$Range.x=='invasive',]
Anova(lm(sqrt(t4t2.leaves)~ CLIMPC2 , data = Invasive), type="III")
summary(lm(sqrt(t4t2.leaves)~ CLIMPC2 , data = Invasive))


#t4.leaves - range pc2 interaction
Anova(lm(sqrt(t4.leaves)~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC2+ Treatment:Range.x + Treatment:CLIMPC2 + Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova), type="III")

m<-lm(sqrt(t4.leaves)~ Treatment+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = ALL_mannova)
interactionMeans(m, factor="Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m)

testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(ALL_mannova$CLIMPC2[which(ALL_mannova$Range.x=="native")]))) 

Native<- ALL_mannova[ALL_mannova$Range.x=='native',]
Anova(lm(sqrt(t4.leaves)~ CLIMPC2 , data = Native), type="III")
summary(lm(sqrt(t4.leaves)~ CLIMPC2 , data = Native))

Invasive<- ALL_mannova[ALL_mannova$Range.x=='invasive',]
Anova(lm(sqrt(t4.leaves)~ CLIMPC2 , data = Invasive), type="III")
summary(lm(sqrt(t4.leaves)~ CLIMPC2 , data = Invasive))



#t4 stem  - NO effect interactions or PC2 
Anova(lm(t4.stem~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Treatment:CLIMPC2+Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Treatment:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC2+Range.x:CLIMPC2, data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment+Range.x+CLIMPC2, data = ALL_mannova), type="III")


#Below  - NO effect of range or pc2 with below
Anova(lm(Below~ Treatment*Range.x*CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Treatment:CLIMPC2+Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Treatment:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC2+Treatment:Range.x+Range.x:CLIMPC2 , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x+CLIMPC2, data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment+Range.x, data = ALL_mannova), type="III")



#covariate 
covar <- agg_phen_mean[ agg_phen_mean$Treatment != "Shading",which(names(agg_phen_mean) %in% c("t0.leaves", "Treatment", "Pop.ID"))]
ALL_mannova$ID<-as.factor(paste(ALL_mannova$Treatment,ALL_mannova$Pop.ID, sep="_") )
covar$ID<-as.factor(paste(covar$Treatment,covar$Pop.ID, sep="_") )
ALL_mannova<-merge(ALL_mannova, covar, by="ID")


m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves*Range.x*Treatment.y, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ Range.x*t0.leaves*Treatment.y, data=ALL_mannova)
summary(m.man, test="Wilks")
m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves*Treatment.y*Range.x, data=ALL_mannova)
summary(m.man, test="Wilks")


m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves+Range.x+Treatment.y+Range.x:Treatment.y+Range.x:t0.leaves+Treatment.y:t0.leaves, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves+Range.x+Treatment.y+Range.x:Treatment.y+Treatment.y:t0.leaves, data=ALL_mannova)
summary(m.man, test="Wilks")

m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves+Range.x+Treatment.y+Range.x:Treatment.y+Treatment.y:t0.leaves, data=ALL_mannova)
summary(m.man, test="Wilks")


m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves+Range.x+Treatment.y+Treatment.y:t0.leaves, data=ALL_mannova)
summary(m.man, test="Wilks")


m.man<-manova(cbind(punch, sqrt(t4.shoots), sqrt(t4t2.leaves), sqrt(t4.leaves), t4.stem, Below) ~ t0.leaves+Range.x+Treatment.y, data=ALL_mannova)
summary(m.man, test="Wilks")


Anova(lm(log(Below)~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(log(Below)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")


# punch  - NO effect of range with or interaction
Anova(lm(punch~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(punch~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x , data = ALL_mannova), type="III")

# t4.shoots - NO effect of range with biomass
Anova(lm(sqrt(t4.shoots)~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment.y+Range.x+t0.leaves+Treatment.y:t0.leaves, data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.shoots)~ Treatment.y+Range.x+t0.leaves, data = ALL_mannova), type="III")

Anova(lm(sqrt(t4t2.leaves)~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x, data = ALL_mannova), type="III")
Anova(lm(sqrt(t4t2.leaves)~ Treatment.y+Range.x+t0.leaves, data = ALL_mannova), type="III")

Anova(lm(sqrt(t4.leaves)~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x, data = ALL_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Treatment.y+Range.x+t0.leaves, data = ALL_mannova), type="III")

#Below  - NO effect of range with biomass
Anova(lm(Below~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment.y+Range.x+t0.leaves+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(Below~ Treatment.y+Range.x+t0.leaves, data = ALL_mannova), type="III")

Anova(lm(t4.stem~ Treatment.y*Range.x*t0.leaves , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves+Range.x:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment.y+Range.x+t0.leaves+Treatment.y:Range.x+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment.y+Range.x+t0.leaves+Treatment.y:t0.leaves , data = ALL_mannova), type="III")
Anova(lm(t4.stem~ Treatment.y+Range.x+t0.leaves, data = ALL_mannova), type="III")



##PC1
##plot t4t2.leaves
pc_meandummy<- ALL_mannova[ ,which(names(ALL_mannova) %in% c("t4t2.leaves","CLIMPC1","Range.x","Treatment.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
pc_meandummy$group<-paste(pc_meandummy$Treatment.x, pc_meandummy$Range.x, sep="_")
mod<- lm(t4t2.leaves~Range.x*Treatment.x*CLIMPC1, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

A<-  ggplot(pc_meandummy,aes(x=CLIMPC1, y=t4t2.leaves, color=factor(group), shape=factor(group))) +
  geom_point(size=2, aes(fill=factor(group))) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+  
  theme_bw()+  ylab("t4t2.leaves")+ xlab("CLIMPC1") +
  theme(panel.grid.minor = element_blank())+ 
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(group)), alpha=0.3, lty=0) 

pc_meandummy$CLIMPC1<-as.numeric(pc_meandummy$CLIMPC1)
pc_meandummy$group<-as.factor(pc_meandummy$group)
m<- lm(t4t2.leaves~group*CLIMPC1, data = pc_meandummy)
interactionMeans(m, factor="group", slope= "CLIMPC1")
testInteractions(m, pairwise = "group", slope = "CLIMPC1")
testInteractions(m, pairwise = "group", covariates = c(CLIMPC1=max(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat

pc_meandummy$CLIMPC1<-as.numeric(pc_meandummy$CLIMPC1)
pc_meandummy$Range.x<-as.factor(pc_meandummy$Range.x)
m<- lm(sqrt(t4t2.leaves)~ Treatment.x+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = pc_meandummy)
interactionMeans(m, factor="Range.x", slope= "CLIMPC1")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC1")
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat


##plot t4.leaves
 
  pc_meandummy<- ALL_mannova[ ,which(names(ALL_mannova) %in% c("t4.leaves","CLIMPC1","Range.x","Treatment.x"))]
  pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
  pc_meandummy$group<-paste(pc_meandummy$Treatment.x, pc_meandummy$Range.x, sep="_")
  mod<- lm(t4.leaves~Range.x*Treatment.x*CLIMPC1, data = pc_meandummy)
  fitted <- predict(mod, interval="confidence")
  pc_meandummy<-cbind(pc_meandummy, fitted)
  
 B<- ggplot(pc_meandummy,aes(x=CLIMPC1, y=t4.leaves, color=factor(group), shape=factor(group))) +
    geom_point(size=2, aes(fill=factor(group))) +
    scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+  
    theme_bw()+  ylab("t4.leaves")+ xlab("CLIMPC1") +
    theme(panel.grid.minor = element_blank())+ 
    geom_line(data=pc_meandummy, aes(y=fit)) + 
    geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(group)), alpha=0.3, lty=0)  

 pc_meandummy$CLIMPC1<-as.numeric(pc_meandummy$CLIMPC1)
 pc_meandummy$Range.x<-as.factor(pc_meandummy$Range.x)
 m<- lm(sqrt(t4.leaves)~ Treatment.x+Range.x+CLIMPC1+ Range.x:CLIMPC1 , data = pc_meandummy)
 interactionMeans(m, factor="Range.x", slope= "CLIMPC1")
 testInteractions(m, pairwise = "Range.x", slope = "CLIMPC1")
 testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
 #testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
 testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=min(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
 #testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC1=max(pc_meandummy$CLIMPC1[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
 
 
pdf("~/Dropbox/Documents/Canada_Thistle_CG/TREAT_NOSHADE_trait_pc1.pdf", height=5, width=8)
ggarrange(A, B,  
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()





##PC2
##plot t4t2.leaves
pc_meandummy<- ALL_mannova[ ,which(names(ALL_mannova) %in% c("t4t2.leaves","CLIMPC2","Range.x","Treatment.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
pc_meandummy$group<-paste(pc_meandummy$Treatment.x, pc_meandummy$Range.x, sep="_")
mod<- lm(t4t2.leaves~Range.x*Treatment.x*CLIMPC2, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)


pc_meandummy$CLIMPC2<-as.numeric(pc_meandummy$CLIMPC2)
pc_meandummy$Range.x<-as.factor(pc_meandummy$Range.x)
m<- lm(sqrt(t4t2.leaves)~ Treatment.x+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = pc_meandummy)
interactionMeans(m, factor="Range.x", slope= "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat



A<-  ggplot(pc_meandummy,aes(x=CLIMPC2, y=t4t2.leaves, color=factor(group), shape=factor(group))) +
  geom_point(size=2, aes(fill=factor(group))) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+  
  theme_bw()+  ylab("t4t2.leaves")+ xlab("CLIMPC2") +
  theme(panel.grid.minor = element_blank())+ 
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(group)), alpha=0.3, lty=0) 

##plot t4.leaves

pc_meandummy<- ALL_mannova[ ,which(names(ALL_mannova) %in% c("t4.leaves","CLIMPC2","Range.x","Treatment.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
pc_meandummy$group<-paste(pc_meandummy$Treatment.x, pc_meandummy$Range.x, sep="_")
mod<- lm(t4.leaves~Range.x*Treatment.x*CLIMPC2, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

pc_meandummy$CLIMPC2<-as.numeric(pc_meandummy$CLIMPC2)
pc_meandummy$Range.x<-as.factor(pc_meandummy$Range.x)
m<- lm(t4.leaves~ Treatment.x+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = pc_meandummy)
interactionMeans(m, factor="Range.x", slope= "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat


B<- ggplot(pc_meandummy,aes(x=CLIMPC2, y=t4.leaves, color=factor(group), shape=factor(group))) +
  geom_point(size=2, aes(fill=factor(group))) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+  
  theme_bw()+  ylab("t4.leaves")+ xlab("CLIMPC2") +
  theme(panel.grid.minor = element_blank())+ 
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(group)), alpha=0.3, lty=0)  

pc_meandummy<- ALL_mannova[ ,which(names(ALL_mannova) %in% c("punch","CLIMPC2","Range.x","Treatment.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
pc_meandummy$group<-paste(pc_meandummy$Treatment.x, pc_meandummy$Range.x, sep="_")
mod<- lm(punch~Range.x*Treatment.x*CLIMPC2, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

pc_meandummy$CLIMPC2<-as.numeric(pc_meandummy$CLIMPC2)
pc_meandummy$Range.x<-as.factor(pc_meandummy$Range.x)
m<- lm(punch~ Treatment.x+Range.x+CLIMPC2+ Range.x:CLIMPC2 , data = pc_meandummy)
interactionMeans(m, factor="Range.x", slope= "CLIMPC2")
testInteractions(m, pairwise = "Range.x", slope = "CLIMPC2")
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="native")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=min(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat
testInteractions(m, pairwise = "Range.x", covariates = c(CLIMPC2=max(pc_meandummy$CLIMPC2[which(pc_meandummy$Range.x=="invasive")]))) #Lowest AU lat



C<- ggplot(pc_meandummy,aes(x=CLIMPC2, y=punch, color=factor(group), shape=factor(group))) +
  geom_point(size=2, aes(fill=factor(group))) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+  
  theme_bw()+  ylab("punch")+ xlab("CLIMPC2") +
  theme(panel.grid.minor = element_blank())+ 
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(group)), alpha=0.3, lty=0)  


pdf("~/Dropbox/Documents/Canada_Thistle_CG/TREAT_NOSHADE_trait_pc2.pdf", height=5, width=8)
ggarrange(A, B,  
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 1)
dev.off()

#Test interaction at min and max PC1 values Range difference
#Slope different from 0 and different from eachother

#survival analysis Shade