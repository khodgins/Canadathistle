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

#get the control treatment and reduce the trait number by removing highly correlated traits
pc_mean_sub_C <- agg_phen_mean[agg_phen_mean$Treatment=="Control",28:68]
corrs <- rcorr(as.matrix(pc_mean_sub_C), type="spearman")
corrs_r <- data.frame(corrs$r)
corrs_p <- data.frame(corrs$P)
corrs_pbonf <- data.frame(corrs$p_bonf)
corr <- flattenCorrMatrix(corrs_r,corrs_p)

corr$p_bonf <- p.adjust(corr$p, method= "bonferroni")
corr[corr$cor < 0.70,]

#remove correlated variables
variable_cor <- corrs_r
variable_cor[upper.tri(corrs_r)] <- 0
diag(variable_cor) <- 0
variable.new <- colnames(pc_mean_sub_C)[!apply(variable_cor,2,function(x) any(abs(x) > 0.7))]
agg_phen_mean_sub_C <- agg_phen_mean[ agg_phen_mean$Treatment=="Control",which(names(agg_phen_mean) %in% variable.new)]

#check to make sure they are gone
corrs <- rcorr(as.matrix(agg_phen_mean_sub_C), type="spearman")
corrs_r <- data.frame(corrs$r)
corrs_p <- data.frame(corrs$P)
corrs_pbonf <- data.frame(corrs$p_bonf)
corr <- flattenCorrMatrix(corrs_r,corrs_p)
corr$p_bonf <- p.adjust(corr$p, method= "bonferroni")
corr[corr$cor > 0.7,]
```
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

agg_phen_mean_all_C<-na.omit(cbind(agg_phen_mean$Range[agg_phen_mean$Treatment=="Control"], agg_phen_mean[agg_phen_mean$Treatment=="Control", 1:2], agg_phen_mean[agg_phen_mean$Treatment=="Control", 31:dim(agg_phen_mean)[2]-3]))
colnames(agg_phen_mean_all_C)[1] <- c("Range")
agg_phen_mean_all_C<-merge(agg_phen_mean_all_C, popID, by="Pop.ID")
agg_phen_mean_all_C$group<-paste(agg_phen_mean_all_C$Range,agg_phen_mean_all_C$Treatment, sep="_")
pc_mean_sub<-agg_phen_mean_all_C[,6:dim(agg_phen_mean_all_C)[2]-2]
trait_pca <- prcomp((pc_mean_sub), scale=TRUE)
agg_phen_mean_all_C$Range<-as.factor(agg_phen_mean_all_C$Range)
autoplot(trait_pca,	type="obs")
pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_control_trait_all.pdf", height=5, width=8)
autoplot(trait_pca, data=agg_phen_mean_all_C, colour="Range",label.label=agg_phen_mean_all_C$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
dev.off()
aload <- abs(trait_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(trait_pca)
agg_phen_mean_all_C$TRAITcPC1 <- trait_pca$x[,1]
agg_phen_mean_all_C$TRAITcPC2 <- trait_pca$x[,2]
agg_phen_mean_all_C$TRAITcPC3 <- trait_pca$x[,3]

pc_mean_subc$TRAITPC1 <- trait_pca$x[,1]
pc_mean_subc$TRAITPC2 <- trait_pca$x[,2]
pc_mean_subc$TRAITPC3 <- trait_pca$x[,3]

corrs<-rcorr(as.matrix(pc_mean_sub, type="spearman"))
corrs_r <- as.matrix(corrs$r)


pdf("~/Dropbox/Documents/Canada_Thistle_CG/TRAIT_all_C_cor.pdf", height=5, width=8)
corrplot(corrs_r)
dev.off()


#PCA reduced corrlation traits

agg_phen_mean_sub_C<-na.omit(cbind(agg_phen_mean$Range[agg_phen_mean$Treatment=="Control"], agg_phen_mean[agg_phen_mean$Treatment=="Control", 1:2], agg_phen_mean[ agg_phen_mean$Treatment=="Control",which(names(agg_phen_mean) %in% variable.new)]))
colnames(agg_phen_mean_sub_C)[1] <- c("Range")
agg_phen_mean_sub_C<-merge(agg_phen_mean_sub_C, popID, by="Pop.ID")
agg_phen_mean_sub_C$Range<-ifelse(agg_phen_mean_sub_C$Range==1,"native", "invasive" )
agg_phen_mean_sub_C$group<-paste(agg_phen_mean_sub_C$Range,agg_phen_mean_sub_C$Treatment, sep="_")
pc_mean_sub<-agg_phen_mean_sub_C[,6:dim(agg_phen_mean_sub_C)[2]-2]
trait_pca <- prcomp((pc_mean_sub), scale=TRUE)
agg_phen_mean_sub_C$Range<-as.factor(agg_phen_mean_sub_C$Range)
autoplot(trait_pca,	type="obs")
#pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_control_trait_sub.pdf", height=5, width=8)
autoplot(trait_pca, data=agg_phen_mean_sub_C, colour="Range", label.label= agg_phen_mean_all_C$Pop.acronym, loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#dev.off()
aload <- abs(trait_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(trait_pca)
agg_phen_mean_sub_C$TRAITcPC1 <- trait_pca$x[,1]
agg_phen_mean_sub_C$TRAITcPC2 <- trait_pca$x[,2]
agg_phen_mean_sub_C$TRAITcPC3 <- trait_pca$x[,3]

pc_mean_sub$TRAITcPC1 <- trait_pca$x[,1]
pc_mean_sub$TRAITcPC2 <- trait_pca$x[,2]
pc_mean_sub$TRAITcPC3 <- trait_pca$x[,3]

corrs<-rcorr(as.matrix(pc_mean_sub, type="spearman"))
corrs_r <- as.matrix(corrs$r)
pdf("~/Dropbox/Documents/Canada_Thistle_CG/TRAIT_sub_C_cor.pdf", height=5, width=8)
corrplot(corrs_r)
dev.off()

Control_mannova<-merge(bioclim, agg_phen_mean_sub_C, by="Pop.ID")
hist(Control_mannova$flowering_day)
hist(sqrt(Control_mannova$flowering_day))
hist(Control_mannova$punch)
hist(Control_mannova$t1t0.leaves)
hist(log(Control_mannova$t1t0.leaves))
hist(Control_mannova$t1.ellipse)
hist(log(Control_mannova$t1.height))
hist(Control_mannova$t3.ellipse)
hist(Control_mannova$flowers)
hist(Control_mannova$t4.shoots)
hist(log(Control_mannova$t4.height))
hist(Control_mannova$Above)
hist(log(Control_mannova$Below))
hist(sqrt(Control_mannova$t4.leaves))
hist(sqrt(Control_mannova$t4.stem))

#range and pc1 - interaction not significant so remove
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Range.x*CLIMPC1, data=Control_mannova)

summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

#range and pc1 (range only significant depending on the order with pc1)
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Range.x+CLIMPC1, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

#range and pc1
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ CLIMPC1+Range.x, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

##do not use pc1 as range and pc1 completely confounded
#test if PC1 sig within each range
invasive<-Control_mannova[Control_mannova$Range.x=="invasive",]
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~  CLIMPC1, data=invasive)
summary(m.man, test="Wilks")
native<-Control_mannova[Control_mannova$Range.x=="native",]
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~  CLIMPC1, data=native)
summary(m.man, test="Wilks")


#range and pc2 (range significant not pc2 or interaction )
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ CLIMPC2*Range.x, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

#range signficant even when the order in the model was varied
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Range.x+CLIMPC2, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ CLIMPC2+Range.x, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

#range -  significant
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Range.x, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") #every trait effect

#Lat and range effects significant both ways (interaction not significant and removed)
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem))~ Latitude*Range.x, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") 
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Range.x+Latitude, data=Control_mannova)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") 
native<-Control_mannova[Control_mannova$Range.x=="native",]
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Latitude, data=native)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") 
invasive<-Control_mannova[Control_mannova$Range.x=="invasive",]
m.man<-manova(cbind(sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Latitude, data=invasive)
summary(m.man, test="Wilks")
summary.aov(m.man, test="Wilks") 


#test trait differences between ranges
Anova(lm(TRAITPC1~ Range.x , data = Control_mannova), type="III")
Anova(lm(TRAITPC2~ Range.x , data = Control_mannova), type="III")
m<-lm(TRAITPC2~ Range.x , data = Control_mannova)
marginal=lsmeans(m,~Range.x, type = "response")

sqrt(flowering_day),  punch, log(t1t0.leaves), t1.ellipse, log(t1.height), t3.ellipse, flowers, t4.shoots, log(t4.height), Above, log(Below), sqrt(t4.leaves), sqrt(t4.stem)) ~ Latitude, data=invasive)

Anova(lm(sqrt(flowering_day)~ Range.x , data = Control_mannova), type="III")
m<-lm(sqrt(flowering_day)~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(punch~ Range.x , data = Control_mannova), type="III")
m<-lm(punch~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(log(t1t0.leaves)~ Range.x , data = Control_mannova), type="III")
m<-lm(log(t1t0.leaves)~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(t1.ellipse~ Range.x , data = Control_mannova), type="III")
m<-lm(t1.ellipse~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(log(t1.height)~ Range.x , data = Control_mannova), type="III")
m<-lm(log(t1.height)~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(t3.ellipse~ Range.x , data = Control_mannova), type="III")
m<-lm(t3.ellipse~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(flowers~ Range.x , data = Control_mannova), type="III")
m<-lm(flowers~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(t4.shoots~ Range.x , data = Control_mannova), type="III")
m<-lm(t4.shoots~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(log(t4.height)~ Range.x , data = Control_mannova), type="III")
m<-lm(log(t4.height)~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(Above~ Range.x , data = Control_mannova), type="III")
m<-lm(Above~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(log(Below)~ Range.x , data = Control_mannova), type="III")
m<-lm(log(Below)~ Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(t1.ellipse~  Range.x , data = Control_mannova), type="III")
m<-lm(t1.ellipse~  Range.x , data = Control_mannova)
lsmeans(m,~Range.x, type = "response")

Anova(lm(sqrt(t4.leaves)~ Range.x , data = Control_mannova), type="III") 
m<-lm(sqrt(t4.leaves)~ Range.x , data = Control_mannova)
marginal=lsmeans(m,~Range.x, type = "response")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)

#plot lsmeans and se
ggplot(CLD,
       aes(x     = Range.x,
           y     = response,
           label = NA)) +
  
  geom_point(shape  = 15,
             size   = 4) +
  
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7) +
  
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("t4 leaf number") +
  xlab("range") +
  ## ylim(15,18)+
  geom_text(nudge_x = c(0, 0, 0, 0),
            nudge_y = c(4, 4, 4, 4),
            color   = "black")



Anova(lm(sqrt(t4.stem)~ Range.x , data = Control_mannova), type="III") 
m<-lm(sqrt(t4.stem)~ Range.x , data = Control_mannova)
marginal=lsmeans(m,~Range.x, type = "response")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)

#plot lsmeans and se
ggplot(CLD,
       aes(x     = Range.x,
           y     = response,
           label = NA)) +
  
  geom_point(shape  = 15,
             size   = 4) +
  
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7) +
  
  theme_bw() +
 theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
       plot.caption = element_text(hjust = 0)) +
  ylab("stem diameter (mm)") +
  xlab("range") +
 ## ylim(15,18)+
  geom_text(nudge_x = c(0, 0, 0, 0),
            nudge_y = c(4, 4, 4, 4),
            color   = "black")



#NS
Anova(lm(sqrt(flowering_day)~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(sqrt(flowering_day)~ Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(punch~  Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(punch~  Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(log(t1t0.leaves)~  Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(log(t1t0.leaves)~  Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(log(t1.height)~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(log(t1.height)~ Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(t3.ellipse~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(t3.ellipse~ Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(t4.shoots~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(t4.shoots~ Range.x+Latitude , data = Control_mannova), type="III")

Anova(lm(log(t4.height)~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(log(t4.height)~ Range.x+Latitude , data = Control_mannova), type="III")


#Lat sig
Anova(lm(flowers~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(flowers~ Range.x+Latitude , data = Control_mannova), type="III")
summary(lm(flowers~ Range.x*Latitude , data = Control_mannova))
summary(lm(flowers~ Range.x+Latitude , data = Control_mannova))
m<-lm(flowers~ Range.x+Latitude , data = Control_mannova)
interactionMeans(m, factor="Range.x", slope = "Latitude")


Anova(lm(Above~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(Above~ Range.x+Latitude , data = Control_mannova), type="III")

#Int sig
Anova(lm(log(Below)~ Range.x*Latitude , data = Control_mannova), type="III")

summary(lm(log(Below)~ Range.x*Latitude , data = Control_mannova))
summary(lm(log(Below)~ Range.x+Latitude , data = Control_mannova))
m<-lm(log(Below)~ Range.x*Latitude , data = Control_mannova)
interactionMeans(m, factor="Range.x", slope = "Latitude")
testInteractions(m, pairwise = "Range.x", slope = "Latitude")
testInteractions(m)
testInteractions(m, pairwise = "Range.x", covariates = c(Latitude=max(Control_mannova$Latitude[which(Control_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(Latitude=min(Control_mannova$Latitude[which(Control_mannova$Range.x=="invasive")]))) 


Native<- Control_mannova[Control_mannova$Range.x=='native',]
Anova(lm(log(Below)~ Latitude , data = Native), type="III")
summary(lm(log(Below)~ Latitude , data = Native))


Invasive<- Control_mannova[Control_mannova$Range.x=='invasive',]
Anova(lm(log(Below)~ Latitude , data = Invasive), type="III")
summary(lm(log(Below)~ Latitude , data = Invasive))

Anova(lm(t1.ellipse~  Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(t1.ellipse~  Range.x+Latitude , data = Control_mannova), type="III")

#Both sig
Anova(lm(sqrt(t4.stem)~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(sqrt(t4.stem)~ Range.x+Latitude , data = Control_mannova), type="III")

m<-lm(sqrt(t4.stem)~ Range.x+Latitude , data = Control_mannova)
interactionMeans(m, factor="Range.x", slope = "Latitude")
testInteractions(m, pairwise = "Range.x", slope = "Latitude")
testInteractions(m)

testInteractions(m, pairwise = "Range.x", covariates = c(Latitude=max(Control_mannova$Latitude[which(Control_mannova$Range.x=="invasive")]))) 
testInteractions(m, pairwise = "Range.x", covariates = c(Latitude=min(Control_mannova$Latitude[which(Control_mannova$Range.x=="invasive")]))) 


Anova(lm(sqrt(t4.leaves)~ Range.x*Latitude , data = Control_mannova), type="III")
Anova(lm(sqrt(t4.leaves)~ Range.x+Latitude , data = Control_mannova), type="III")
m<-lm(sqrt(t4.leaves)~ Range.x+Latitude , data = Control_mannova)
interactionMeans(m, factor="Range.x", slope = "Latitude")
testInteractions(m, pairwise = "Range.x", slope = "Latitude")
testInteractions(m)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#Stem
pc_meandummy<- Control_mannova[ ,which(names(Control_mannova) %in% c("t4.stem","Latitude","Range.x"))]
colnames(pc_meandummy)[1]<-"Range"
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
mod<- lm(sqrt(t4.stem)~Range+Latitude, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

A<-ggplot(pc_meandummy,aes(x=Latitude, y=sqrt(t4.stem), color=Range, shape=Range)) +
  geom_point(size=2, aes(fill=Range)) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+ 
  #scale_color_manual(values=c(12,12,12,12,22,22,22,22))+
  theme_bw()+  ylab("Stem diameter")+ xlab("Latitude") +
  theme(panel.grid.minor = element_blank(), legend.position="top")+
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(Range)), alpha=0.3, lty=0) 

pc_meandummy<- Control_mannova[ ,which(names(Control_mannova) %in% c("t4.leaves","Latitude","Range.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
mod<- lm(sqrt(t4.leaves)~Range.x+Latitude, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

B<-ggplot(pc_meandummy,aes(x=Latitude, y=sqrt(t4.leaves), color=factor(Range.x), shape=factor(Range.x))) +
  geom_point(size=2, aes(fill=factor(Range.x)),show.legend = FALSE) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+ 
  #scale_color_manual(values=c(12,12,12,12,22,22,22,22))+
  theme_bw()+  ylab("Leaf number")+ xlab("Latitude") +
  theme(panel.grid.minor = element_blank())+
  geom_line(data=pc_meandummy, aes(y=fit),show.legend = FALSE) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(Range.x)), alpha=0.3, lty=0,show.legend = FALSE) 



pc_meandummy<- Control_mannova[ ,which(names(Control_mannova) %in% c("flowers","Latitude","Range.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
mod<- lm(flowers~Range.x+Latitude, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

C<-ggplot(pc_meandummy,aes(x=Latitude, y=flowers, color=factor(Range.x), shape=factor(Range.x))) +
  geom_point(size=2, aes(fill=factor(Range.x)),show.legend = FALSE) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+ 
  #scale_color_manual(values=c(12,12,12,12,22,22,22,22))+
  theme_bw()+  ylab("Flower number")+ xlab("Latitude") +
  theme(panel.grid.minor = element_blank())+
  geom_line(data=pc_meandummy, aes(y=fit),show.legend = FALSE) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(Range.x)), alpha=0.3, lty=0,show.legend = FALSE) 

pc_meandummy<- Control_mannova[ ,which(names(Control_mannova) %in% c("Below","Latitude","Range.x"))]
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
mod<- lm(log(Below)~Range.x*Latitude, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

D<-ggplot(pc_meandummy,aes(x=Latitude, y=log(Below), color=factor(Range.x), shape=factor(Range.x))) +
  geom_point(size=2, aes(fill=factor(Range.x)),show.legend = FALSE) +
  scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+ 
  #scale_color_manual(values=c(12,12,12,12,22,22,22,22))+
  theme_bw()+  ylab("Below ground biomass")+ xlab("Latitude") +
  theme(panel.grid.minor = element_blank())+
  geom_line(data=pc_meandummy, aes(y=fit),show.legend = FALSE) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(Range.x)), alpha=0.3, lty=0,show.legend = FALSE) 

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Control_trait_lat.pdf", height=5, width=8)
ggarrange(A, B, C, D ,  
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()














