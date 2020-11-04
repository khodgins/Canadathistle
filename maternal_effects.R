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
maternalraw <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/Maternal.effects_long.txt", header = T)

##mixed models Treatment*Range

clone<-maternalraw[maternalraw$Generation=="Clone" & maternalraw$Treatment !="Shading",]

#INT SIG
step(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.leaves
step(lmer(t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.ellipse
step(lmer(t1.ellipse~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.ellipse~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.hheight
step(lmer(t1.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.height~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.shoots 
step(lmer(t2.shoots~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.shoots~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.leaves 
step(lmer(t2.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t1.leaves 
step(lmer(t2t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t0.leaves
step(lmer(t2t0.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t0.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#ellipse 
step(lmer(ellipse~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(ellipse~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.height 
step(lmer(t2.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.height~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.shoots
step(lmer(t3.shoots~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.shoots~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.leaves 
step(lmer(t3.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.flowers 
step(lmer(t3.flowers~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.flowers~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.stem
step(lmer(t3.stem~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.stem~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")


#t3.height
step(lmer(t3.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.height~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#Above 
step(lmer(Above~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Above~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")


#Below - Range ns
step(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=clone), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#no significant effects of range for any of the clones or a treatment*range interaction


#get same moms from main experiment
mom<-unique(clone$Mom)
#see if an interaction involving experiment for any trait
#experiment * range interaction or  experiment*treatment 

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

mat<-phenraw[phenraw$Mom %in% mom,]
#mat<-mat[ mat$Treatment == "Control" | mat$Treatment == "Nutrient" |  mat$Treatment == "Shading",]
mat<-mat[ mat$Treatment == "Control" | mat$Treatment == "Nutrient" ,]

#stil patterns in the significantly different traits? NO

step(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.leaves - RANGE SIG and RANGE*INT MARG SIG
step(lmer(t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.ellipse
step(lmer(t1.ellipse~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.ellipse~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.height INT SIG RANGE NOT
step(lmer(t1.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.shoots 
step(lmer(t2.shoots~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.shoots~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.leaves 
step(lmer(t2.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t1.leaves  RANGE AND INT SIGN
step(lmer(t2t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t0.leaves  INT SIGN
step(lmer(t2t0.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t0.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

####ellipse 
step(lmer(t2.ellipse~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.ellipse~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.height 
step(lmer(t2.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.height~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.shoots
step(lmer(t3.shoots~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.shoots~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.leaves 
step(lmer(t3.leaves~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.leaves~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.flowers 
step(lmer(flowers~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(flowers~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.stem
step(lmer(t4.stem~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t4.stem~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.height
step(lmer(t3.height~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.height~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#Above 
step(lmer(Above~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Above~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#Below - Range ns
step(lmer(Below~ Treatment*Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below~ Treatment+Range+ (1|Pop.ID) +(1|Mom), data=mat), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")



#make combined dataset
library(gtools)
mat2$Generation<-"Mat"
mat2$Generation <- factor(mat2$Generation)
clone$Generation <- factor(clone$Generation)
colnames(mat2)[69]<-"stem"
colnames(clone)[30]<-"flowers"
colnames(clone)[25]<-"t2.ellipse"
colnames(clone)[31]<-"stem"
combo<-smartbind(mat2,clone)


#t1.leaves - RANGE SIG Range:Generation Treatment:Generation
step(lmer(t1.leaves~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.leaves ~ Treatment + Range + Generation + Treatment:Range + Treatment:Generation + Range:Generation +  (1|Pop.ID) + (1|Mom) , data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.ellipse Range:Generation Treatment:Generation
step(lmer(t1.ellipse~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.ellipse~ Treatment + Range + Generation + (1 | Pop.ID) + (1 | Mom) + Treatment:Generation + Range:Generation, data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t1.height G*R*T R*G T*G
step(lmer(t1.height~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.height~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.shoots NO INT OR RANGE
step(lmer(t2.shoots~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.shoots~ Treatment+Range+ Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.leaves  NO INT OR RANGE
step(lmer(t2.leaves~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.leaves~ Treatment+Range+Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t1.leaves  G*R*T R*G T*G
step(lmer(t2t1.leaves~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~ Treatment*Range*Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2t0.leaves  R*T NO interaction with GEN
step(lmer(t2t0.leaves~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t0.leaves~ Treatment+Range+Generation+Treatment:Range+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.ellipse T*G G T; NO RANGE
step(lmer(t2.ellipse~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.ellipse~ Treatment+Range+Generation + Treatment:Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t2.height T; NO RANGE OR GEN
step(lmer(t2.height~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2.height~ Treatment+Range+Generation +(1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.shoots T G; NO RANGE OR INT
step(lmer(t3.shoots~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.shoots~ Treatment+Range+Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.leaves 
step(lmer(t3.leaves~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.leaves~ Treatment+Range+Generation + Treatment:Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#flowers T G; NO RANGE OR INT
step(lmer(flowers~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(flowers~ Treatment+Range+Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#stem T G NO R or INT
step(lmer(stem~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(stem~ Treatment+Range+Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#t3.height NO INT T NO G
step(lmer(t3.height~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t3.height~ Treatment+Range+Generation + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#Above  T*G G T NO R
step(lmer(Above~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Above~ Treatment+Range+ Generation+Treatment:Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#Below - T*G T*G G T - NO R
step(lmer(Below~ Treatment*Range*Generation+ (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below~ Treatment+Range+ Generation+Treatment:Generation +Treatment:Range + (1|Pop.ID) +(1|Mom), data=combo), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")


########
#take pop*experiment means and see if significant correlations within each treatment

#merge by maternal line paste generation, treatment, line
mat2$group <- paste(mat2$Treatment, mat2$Mom, sep="_")
clone$group <- paste(clone$Treatment, clone$Mom, sep="_")
combo2<-merge(mat2,clone, by="group")
control<-combo2[combo2$Treatment.x=="Control",]
nutrient<-combo2[combo2$Treatment.x=="Nutrient",]

#Below
Anova(lmer(Below.y~ Below.x*Treatment.x*Range.x  + (1|Pop.ID.x) , data=combo2), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

#sign trait 
Anova(lmer(Below.y~ Below.x*Range.x  + (1|Pop.ID.x) , data=control), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below.y~ Below.x+Range.x  + (1|Pop.ID.x) , data=control), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
plot(control$Below.y~ control$Below.x)
Anova(lm(control$Below.y~ control$Below.x))

#NS trait; R*Trait interaction
Anova(lmer(Below.y~ Below.x*Range.x  + (1|Pop.ID.x) , data=nutrient), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below.y~ Below.x+Range.x  + (1|Pop.ID.x) , data=nutrient), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
plot(nutrient$Below.y~ nutrient$Below.x)
Anova(lm(nutrient$Below.y~ nutrient$Below.x))

#Stem - not correlated with anything
step(Anova(lmer(stem.y~ stem.x*Treatment.x*Range.x  + (1|Pop.ID.x) , data=combo2), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F"))

Anova(lmer(stem.y~ stem.x*Range.x  + (1|Pop.ID.x) , data=control), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(stem.y~ stem.x+Range.x  + (1|Pop.ID.x) , data=control), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(stem.y~ stem.x + (1|Pop.ID.x) , data=control), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

plot(control$stem.y~ control$stem.x)
Anova(lm(control$stem.y~ control$stem.x))

#NS trait; R*Trait interaction
Anova(lmer(Below.y~ Below.x*Range.x  + (1|Pop.ID.x) , data=nutrient), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(Below.y~ Below.x+Range.x  + (1|Pop.ID.x) , data=nutrient), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
plot(nutrient$Below.y~ nutrient$Below.x)
Anova(lm(nutrient$Below.y~ nutrient$Below.x))


#t1.leaves - not correlated except for treatment

step(lmer(t1.leaves.y~ t1.leaves.x*Treatment.x*Range.x  + (1|Pop.ID.x) , data=combo2), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t1.leaves.y~ t1.leaves.x+Treatment.x + (1|Pop.ID.x) , data=combo2), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

plot(control$t1.leaves.y~ control$t1.leaves.x)
Anova(lm(control$t1.leaves.y~ control$t1.leaves.x))

plot(nutrient$t1.leaves.y~ nutrient$t1.leaves.x)
Anova(lm(nutrient$t1.leaves.y~ nutrient$t1.leaves.x))



#take population means
combo2_pop<- aggregate(combo2, by= list(combo2$Treatment.x, combo2$Pop.ID.x), FUN=mean, na.rm=T)
plot(combo2_pop$Below.y~ combo2_pop$Below.x)

Anova(lm(Below.y~ Below.x*Group.1 , data=combo2_pop), type="III")



#take pop means of combo data beloe - no interactions
#below
combo_pop<- aggregate(combo, by= list(combo$Treatment, combo$Generation, combo$Pop.ID), FUN=mean, na.rm=T)
colnames(combo_pop)[1]<-"treatment"
colnames(combo_pop)[2]<-"generation"
colnames(combo_pop)[3]<-"range"

#below - G T
Anova(lm(Below~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range , data=combo_pop), type="III")

#t1.leaves         T G R T*G G*R  - Sig Range difference in maternal not in Clone (I<N)
#same as mixed model
Anova(lm(t1.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t1.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t1.leaves~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
m<-lm((t1.leaves)~  treatment+generation+range +treatment:generation +generation:range, data = combo_pop)
marginal=lsmeans(m,~generation:range, type = "response")
 
#t1.ellipse T G T*G G*R NO R - NO differences in range but direction reversed mat N>I clone N<I; Con > Nut (sig in mat and clone)
Anova(lm(t1.ellipse~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t1.ellipse~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t1.ellipse~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
m<-lm((t1.ellipse)~  treatment+generation+range +treatment:generation +generation:range, data = combo_pop)
marginal=lsmeans(m,~generation:range, type = "response")
lsmeans(m,~generation:treatment, type = "response")

#t1.height T G T*G NO R Con > Nut (sig in maternal line not clone)
Anova(lm(t1.height~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t1.height~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t1.height~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t1.height~ treatment+generation+range +treatment:generation    , data=combo_pop), type="III")
m<-lm((t1.height)~  treatment+generation+range +treatment:generation , data = combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")
 
#t2.shoots T G NO R
Anova(lm(t2.shoots~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t2.shoots~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t2.shoots~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t2.shoots~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(t2.shoots~ treatment+generation+range , data=combo_pop), type="III")

#t2.leaves T G NO R
Anova(lm(t2.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t2.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t2.leaves~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t2.leaves~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(t2.leaves~ treatment+generation+range , data=combo_pop), type="III")

#t2t1.leaves T*G NO T G R ; Con > Nut (Sig in clone not mat)
Anova(lm(t2t1.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
m<-lm((t2t1.leaves)~   treatment+generation+range +treatment:generation  , data = combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

#below - G T
Anova(lm(Below~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range , data=combo_pop), type="III")

#above G T T*G
Anova(lm(Above~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range , data=combo_pop), type="III")

#t2t0.leaves G T NO R (Con < Nut)
Anova(lm(t2t0.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t2t0.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t2t0.leaves~ treatment+generation+range +treatment:generation +treatment:range  , data=combo_pop), type="III")
Anova(lm(t2t0.leaves~ treatment+generation+range +treatment:range   , data=combo_pop), type="III")
Anova(lm(t2t0.leaves~ treatment+generation+range   , data=combo_pop), type="III")
m<-lm((t2t0.leaves)~   treatment+generation+range  , data = combo_pop)
marginal=lsmeans(m,~treatment, type = "response")

#above G T T*G
Anova(lm(Stem~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range , data=combo_pop), type="III")

#########CORR
combo_popx <- combo_pop[ ,which(names(combo_pop) %in% c("treatment", "range", "generation", "t1.leaves", "t1.ellipse", "t1.height", "t2.shoots", "t2.leaves", "t2t1.leaves", "t2t0.leaves", "t2.height", "t2.ellipse", "t3.shoots", "t3.leaves", "flowers", "stem", "t3.height", "Above", "Below"))]

#get the control treatment and reduce the trait number by removing highly correlated traits
combo_popx_sub <- combo_popx[,4:19]
corrs <- rcorr(as.matrix(combo_popx_sub), type="spearman")
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
variable.new <- colnames(combo_popx_sub)[!apply(variable_cor,2,function(x) any(abs(x) > 0.7))]
combo_popx_sub_R <- combo_popx_sub[ ,which(names(combo_popx_sub) %in% variable.new)]

#check to make sure they are gone
corrs <- rcorr(as.matrix(combo_popx_sub_R), type="spearman")
corrs_r <- data.frame(corrs$r)
corrs_p <- data.frame(corrs$P)
corrs_pbonf <- data.frame(corrs$p_bonf)
corr <- flattenCorrMatrix(corrs_r,corrs_p)
corr$p_bonf <- p.adjust(corr$p, method= "bonferroni")
corr[corr$cor > 0.7,]

hist(sqrt(combo_popx$t2t1.leaves[combo_popx$treatment=="Control"]))
hist(sqrt(combo_popx$t2t1.leaves[combo_popx$treatment=="Nutrient"]))
hist(sqrt(combo_popx$Above[combo_popx$treatment=="Control"]))
hist(sqrt(combo_popx$Above[combo_popx$treatment=="Nutrient"]))
hist(sqrt(combo_popx$Below[combo_popx$treatment=="Control"]))
hist(sqrt(combo_popx$Below[combo_popx$treatment=="Nutrient"]))

m.man<-manova(cbind(t2t1.leaves,  Above, Below) ~ range*treatment*generation, data=combo_popx )
m.man<-manova(cbind(t2t1.leaves,  Above, Below) ~ range+treatment+generation+range:treatment+range:generation+treatment:generation, data=combo_popx )
m.man<-manova(cbind(t2t1.leaves,  Above, Below) ~ range+treatment+generation+range:treatment+treatment:generation, data=combo_popx )
m.man<-manova(cbind(t2t1.leaves,  Above, Below) ~ range+treatment+generation+treatment:generation, data=combo_popx )
summary(m.man, test="Wilks")


##univariate analysis
#below - G T Control > Nut Mat > Clone
Anova(lm(Below~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Below~ treatment+generation+range , data=combo_pop), type="III")
m<-lm(Below~ treatment+generation+range , data=combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

A<-ggplot(CLD,
       aes(x     = treatment,
           y     = lsmean,
           color= generation,
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
  ylab("below ground biomass") +
  xlab("treatment") +
  geom_text(color   = "black") 



################

#above G T T*G - Con > Nut but both signficantly different
Anova(lm(Above~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
Anova(lm(Above~ treatment+generation+range , data=combo_pop), type="III")
m<-lm(Above~ treatment+generation+range , data=combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

B<-ggplot(CLD,
          aes(x     = treatment,
              y     = lsmean,
              color= generation,
              label=.group)) +
  
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
  ylab("above ground biomass") +
  xlab("treatment") +
  geom_text(color   = "black") 



#t2t1.leaves T*G NO T G R ; Con > Nut (Sig in clone not mat)
Anova(lm(t2t1.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t2t1.leaves~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
m<-lm((t2t1.leaves)~   treatment+generation+range +treatment:generation  , data = combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

C<-ggplot(CLD,
          aes(x     = treatment,
              y     = response,
              color= generation,
              label=.group)) +
  
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
  xlab("treatment") +
  geom_text(color   = "black") 


#stem
Anova(lm(stem~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(stem~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(stem~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(stem~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
m<-lm((stem)~   treatment+generation+range +treatment:generation  , data = combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

D<-ggplot(CLD,
          aes(x     = treatment,
              y     = response,
              color= generation,
              label=.group)) +
  
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
  xlab("treatment") +
  geom_text(color   = "black") 

#final leaves
Anova(lm(t3.leaves~ treatment*generation*range , data=combo_pop), type="III")
Anova(lm(t3.leaves~ treatment+generation+range +treatment:generation +treatment:range  +generation:range , data=combo_pop), type="III")
Anova(lm(t3.leaves~ treatment+generation+range +treatment:generation +generation:range   , data=combo_pop), type="III")
Anova(lm(t3.leaves~ treatment+generation+range +treatment:generation   , data=combo_pop), type="III")
m<-lm((t3.leaves)~   treatment+generation+range +treatment:generation  , data = combo_pop)
marginal=lsmeans(m,~generation:treatment, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

E<-ggplot(CLD,
          aes(x     = treatment,
              y     = response,
              color= generation,
              label=.group)) +
  
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
  ylab("final leaf number") +
  xlab("treatment") +
  geom_text(color   = "black") 

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Maternal_effects.pdf", height=5, width=8)
ggarrange(A, B, C, D ,E,  
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)
dev.off()

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Maternal_effects_nostem.pdf", height=5, width=8)
ggarrange(A, B, C,   
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

######################PLOT SIG EXP*RANGE

step(lmer(final.stem~ Treatment*Range*gen  +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.stem~ Treatment*Range*gen +(1|Pop.ID) , data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

m<-lmer(final.stem~ Treatment*Range*gen +(1|Pop.ID) , data=mat3)
marginal=lsmeans(m,~Treatment:Range:gen, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

ggplot(CLD,
          aes(x     = Treatment,
              y     = lsmean,
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
  ylab("stem diameter") +
  xlab("treatment") +
  geom_text(color   = "black") 










#########
step(lmer(final.leaves~ Treatment*Range*gen+ (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.leaves~  Treatment*Range*gen+ (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen + Range:gen + (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen + Range:gen + (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen +  (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(final.leaves~  Treatment+Range+gen+ Treatment:gen +  (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")


m<-lmer(final.leaves~ Treatment*Range*gen +(1|Pop.ID) , data=mat3)
marginal=lsmeans(m,~Treatment:Range:gen, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

ggplot(CLD,
       aes(x     = Treatment,
           y     = lsmean,
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
  ylab("stem diameter") +
  xlab("treatment") 



step(lmer(t2t1.leaves~ Treatment*Range*gen+ (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~ Treatment*Range*gen+ (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")


m<-lmer(t2t1.leaves~ Treatment*Range*gen +(1|Pop.ID) , data=mat3)
marginal=lsmeans(m,~Treatment:Range:gen, type = "response")

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters)
pd = position_dodge(0.4) 

ggplot(CLD,
       aes(x     = Treatment,
           y     = lsmean,
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
  ylab("stem diameter") +
  xlab("treatment") 





Anova(lmer(t2t1.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen + Range:gen + (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen + Range:gen + (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~  Treatment+Range+gen+ Treatment:Range + Treatment:gen +  (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")
Anova(lmer(t2t1.leaves~  Treatment+Range+gen+ Treatment:gen +  (1|Pop.ID) +(1|Mom), data=mat3), reduce.fixed = TRUE, reduce.random = FALSE, alpha.fixed = 0.05,type =3, ddf="Kenward-Roger",test.statistic = "F")

