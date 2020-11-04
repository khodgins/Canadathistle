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

#delete shading, repeat with shading
#no mixed models on all maternal and clone common traits and see if there are any significant range effects
#combine experiments for common traits and see if there are experiment*range interactions for any traits
#take pop*experiment means and see if significant correlations within each treatment

#read in maternal effects data
germraw <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/germination.txt", header = T)
germraw$Pop.ID <- factor(germraw$Pop.ID)
germraw$Mom.ID <- factor(germraw$Mom.ID)
germraw$Range  <- factor(germraw$Range)
germraw$rate <- as.numeric(germraw$rate) 
Anova(lmer(asin(sqrt(rate))~ Range  + (1|Pop.ID) , data=germraw), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#pop means
germ<-germraw
germ$Range<-ifelse(germraw$Range=="native", 1,0)
germ_mean <- aggregate(germ, by= list(germ$Pop.ID), FUN=mean, na.rm=T)
germ_mean <- germ_mean[,-c(3:4)] # remove empty columns
colnames(germ_mean)[1] <- c("Pop.acronym")
germ_mean$Range<-ifelse(germ_mean$Range==1, "native","invasive")

#still significant
Anova(lm(asin(sqrt(rate))~ Range  , data=germ_mean))  

m<-lm(asin(sqrt(rate))~ Range  , data=germ_mean)
marginal=lsmeans(m,~Range, type = "response")
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

A<-ggplot(CLD,
          aes(x     = c("Introduced", "Native"),
              y     = response,
              color= Range,
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
  ylab("Germination rate") +
  xlab("Range") 

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Germination.pdf", height=5, width=8)
ggarrange(A,    
          labels = c("A"),
          ncol = 1, nrow = 1)
dev.off()



#
bioclim<-read.table("~/Dropbox/Documents/Canada_Thistle_CG/Bioclim.txt", header = T)

bioclim<-bioclim[bioclim$Pop.ID != "150908-1",]
bioclim<-merge(bioclim, popID, by="Pop.ID")
climnolatnolonst_pca <- prcomp((bioclim[,5:23]), scale=T)
climlatst_pca <- prcomp((bioclim[,c(3:23)]), scale=T)

autoplot(climnolatnolonst_pca,	type="obs")
#pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges_nolatlon.pdf", height=5, width=8)
autoplot(climnolatnolonst_pca, data=bioclim, colour="Range",label.label=bioclim$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#autoplot(climnolatnolonst_pca, data=bioclim, colour="Range",label.label="Pop",loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,loadings.label.repel=T,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#dev.off()
#pdf("~/Dropbox/Documents/Canada_Thistle_CG/PCA_ranges.pdf", height=5, width=8)
autoplot(climlatst_pca, data=bioclim, colour="Range",label.label="Pop",loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#autoplot(climlatst_pca, data=bioclim, colour="Range",label.label=bioclim$Pop.acronym,loadings=T,label.size=3, loadings.label=T,loadings.colour=1,loadings.label.colour=1,loadings.label.angle=-10,loadings.label.hjust=-.1,shape=F,  loadings.label.size=3,frame = TRUE, frame.alpha=.05)+theme_bw()
#dev.off()
aload <- abs(climlatst_pca$rotation)
sweep(aload, 2, colSums(aload), "/")
summary(climnolatnolonst_pca)
bioclim$CLIMPC1 <- climnolatnolonst_pca$x[,1]
bioclim$CLIMPC2 <- climnolatnolonst_pca$x[,2]
bioclim$CLIMPC3 <- climnolatnolonst_pca$x[,3]

germ_mean_C<-merge(germ_mean, bioclim, by="Pop.acronym")
germ_mean_C$CLIMPC1<-as.numeric(germ_mean_C$CLIMPC1)

#not correlated with latitude or PC2 but correlated with PC1
Anova(lm(asin(sqrt(rate))~ Range.x *CLIMPC1 , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ Range.x + CLIMPC1 , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ CLIMPC1 , data=germ_mean_C))  
m<-lm(asin(sqrt(rate))~ Range.x + CLIMPC1 , data=germ_mean_C)
m

pc_meandummy<- germ_mean_C[ ,which(names(germ_mean_C) %in% c("rate","CLIMPC1","Range.x"))]
colnames(pc_meandummy)[1]<-"Range"
pc_meandummy<-pc_meandummy[complete.cases(pc_meandummy), ]
mod<- lm(asin(sqrt(rate))~Range+CLIMPC1, data = pc_meandummy)
fitted <- predict(mod, interval="confidence")
pc_meandummy<-cbind(pc_meandummy, fitted)

ggplot(pc_meandummy,aes(x=CLIMPC1, y=asin(sqrt(rate)), color=Range, shape=Range)) +
  geom_point(size=2, aes(fill=Range)) +
  #scale_shape_manual(values=c(21,22,24,25,21,22,24,25))+ 
  #scale_color_manual(values=c(12,12,12,12,22,22,22,22))+
  theme_bw()+  ylab("Germination rate")+ xlab("CLIMPC1") +
  theme(panel.grid.minor = element_blank(), legend.position="top")+
  geom_line(data=pc_meandummy, aes(y=fit)) + 
  geom_ribbon(data=pc_meandummy,aes(ymin=lwr,ymax=upr, fill=factor(Range)), alpha=0.3, lty=0) 

Anova(lm(asin(sqrt(rate))~ Range.x *CLIMPC2 , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ Range.x + CLIMPC2 , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ CLIMPC2 , data=germ_mean_C))

Anova(lm(asin(sqrt(rate))~ Range.x *Latitude   , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ Range.x + Latitude , data=germ_mean_C), type="III")
Anova(lm(asin(sqrt(rate))~ Latitude , data=germ_mean_C))  
