





# check if required packages are installed
list.of.packages <- c("lme4", "lmerTest", "multcomp", "doBy", "car")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load packages
library(lme4)
library(lmerTest)
library(multcomp)
library(doBy)
library(car)



# read data
eeg_r <- read.csv('./data/XPLOWHIGH_r.csv')

eeg_r$ID <- as.factor(eeg_r$ID)
eeg_r$freq <- as.factor(eeg_r$freq)
eeg_r$tone <- ifelse(grepl('^H',eeg_r$condition),'high','low')
eeg_r$tone <- factor(eeg_r$tone, levels=c('high','low'))
eeg_r$rhythm <- ifelse(grepl('_unsyncopated',eeg_r$condition),'unsyncopated','syncopated')
eeg_r$rhythm <- factor(eeg_r$rhythm, levels=c('unsyncopated','syncopated'))
eeg_r$isShifted <- as.factor(ifelse(eeg_r$binshift==0,'no','yes'))



# fit mixed model
eeg_r$isShifted <- as.factor(ifelse(eeg_r$binshift==0,'no','yes'))
eeg_r4anova <- summaryBy(r~ID+tone+rhythm+freq+isShifted, data=eeg_r, FUN=mean, keep.names = 1)

fit_r_full <- lmer(r ~ rhythm*tone*freq*isShifted+(1|ID), data=eeg_r4anova)
Anova(fit_r_full, test='F')





# because there was a marginally significant Baseline x Tone interaction, 
# fit separate models for low and high tone and check that the 
# Frequency x Rhythm x Baseline interaction is significant in both. 

fit_r_LowTone <- lmer(r ~ rhythm*freq*isShifted+(1|ID), data=subset(eeg_r4anova,tone=='low'))
Anova(fit_r_LowTone, test='F')

fit_r_HighTone <- lmer(r ~ rhythm*freq*isShifted+(1|ID), data=subset(eeg_r4anova,tone=='high'))
Anova(fit_r_HighTone, test='F')


