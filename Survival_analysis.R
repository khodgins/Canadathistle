library(survival)
#install.packages("survminer")
library(survminer)
library(lubridate)
#install.packages("coxme") 
library(coxme) 
shading <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/Shading_survival.txt", header = T)
shading$Code<-as.numeric(shading$Code)
shading$Survival<-as.factor(shading$Shading.death)
shading$Day<-as.numeric(shading$Day)
shading$Range<-as.factor(shading$Range)
shading$Pop.ID<-as.factor(shading$Pop.ID)
shading$Shading.t0.leaves<-as.numeric(shading$Shading.t0.leaves)
shading$Survival<-as.numeric(as.character(shading$Survival)) 

fit <- coxph(Surv(Day, Code) ~ cluster(Pop.ID) , data=shading)
fit_full <- coxph(Surv(Day, Code) ~ Range+ cluster(Pop.ID), data=shading)

coxph(Surv(Day, Code) ~ Range+ cluster(Pop.ID), data=shading)

fit_full<-coxme(Surv(Day, Code) ~ Range+ Shading.t0.leaves + (1 | Pop.ID), data=shading)
fit<-coxme(Surv(Day, Code) ~ Shading.t0.leaves + (1 | Pop.ID), data=shading)
anova(fit, fit_full)
fit_full<-coxme(Surv(Day, Code) ~ Range + (1 | Pop.ID), data=shading)
fit<-coxme(Surv(Day, Code) ~ (1 | Pop.ID), data=shading)
anova(fit, fit_full)

summary(fit)
summary(fit)$table
fit<-survdiff(Surv(Day, Code) ~ Range, data=shading)
fit<-survfit(Surv(Day, Code) ~ Range, data=shading)

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Shading_survival.pdf", height=5, width=8)
ggsurvplot(fit,
           pval = FALSE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "Range", # Change risk table color by groups
           surv.median.line = "hv") # Specify median survival

dev.off()




mowing <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/Mowing_survival.txt", header = T)
mowing$Code<-as.numeric(mowing$Code)
mowing$Survival<-as.factor(mowing$Survival)
mowing$Day<-as.numeric(mowing$Day)
mowing$Range<-as.factor(mowing$Range)

mowing$Survival<-as.numeric(as.character(mowing$Survival)) 


fit <- coxph(Surv(Day, Code) ~ cluster(Pop.ID) , data=mowing)
fit_full <- coxph(Surv(Day, Code) ~ Range+ cluster(Pop.ID), data=mowing)

coxph(Surv(Day, Code) ~ Range+ cluster(Pop.ID), data=mowing)

fit_full<-coxme(Surv(Day, Code) ~ Range+ (1 | Pop.ID), data=mowing)
fit<-coxme(Surv(Day, Code) ~ (1 | Pop.ID), data=mowing)
anova(fit, fit_full)

summary(fit)
summary(fit)$table
fit<-survdiff(Surv(Day, Code) ~ Range, data=mowing)
fit<-survfit(Surv(Day, Code) ~ Range, data=mowing)

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Mowing_survival.pdf", height=5, width=8)
ggsurvplot(fit,
           pval = FALSE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "Range", # Change risk table color by groups
           surv.median.line = "hv") # Specify median survival
           
dev.off()

ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 200))

##
drought <-  read.table("~/Dropbox/Documents/Canada_Thistle_CG/Drought_survival.txt", header = T)
drought$Day.Wilt<-as.numeric(drought$Day.Wilt)
drought$First.Wilt<-as.numeric(drought$First.Wilt)
drought$Range<-as.factor(drought$Range)

fit <- survfit(Surv(Day.Wilt, First.Wilt) ~ Range, data=drought)


summary(fit)
summary(fit)$table
survdiff(Surv(Day.Wilt, First.Wilt) ~ Range, data=drought)

ggsurvplot(fit,
           pval = FALSE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "Range", # Change risk table color by groups
           surv.median.line = "hv")

ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 10))


drought$Day.Wilt2<-as.numeric(drought$Day.Wilt2)
drought$Wilt<-as.numeric(drought$Wilt)

fit <- survfit(Surv(Day.Wilt2, Wilt) ~ Range, data=drought)

summary(fit)
summary(fit)$table
survdiff(Surv(Day.Wilt, First.Wilt) ~ Range, data=drought)

ggsurvplot(fit,
           pval = FALSE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "Range", # Change risk table color by groups
           surv.median.line = "hv")

ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 20))

drought$Day.Death<-as.numeric(drought$Day.Death)
drought$Survival<-as.numeric(drought$Survival)

coxph(Surv(Day.Death, Survival) ~ Range+ cluster(Pop.ID), data=drought)

fit_full<-coxme(Surv(Day.Death, Survival) ~ Range+ (1 | Pop.ID), data=drought)
fit<-coxme(Surv(Day.Death, Survival) ~ (1 | Pop.ID), data=drought)
anova(fit, fit_full)


summary(fit)
summary(fit)$table
survdiff(Surv(Day.Death, Survival) ~ Range, data=drought)
fit<-survfit(Surv(Day.Death, Survival) ~ Range, data=drought)

pdf("~/Dropbox/Documents/Canada_Thistle_CG/Drought_survival.pdf", height=5, width=8)
ggsurvplot(fit,
           pval = FALSE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "Range", # Change risk table color by groups
           surv.median.line = "hv")
dev.off()

ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 25))

