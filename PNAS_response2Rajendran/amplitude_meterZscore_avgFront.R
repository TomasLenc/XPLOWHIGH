

# replicate the results with zscores at meter frequencies with frontal channels only






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




# load data
eeg_amps <- read.csv('./data/XPLOWHIGH_avgFront_amplitudes.csv')
eeg_amps$tone <- ifelse(grepl('^H',eeg_amps$condition),'high','low')
eeg_amps$tone <- factor(eeg_amps$tone, levels=c('high','low'))
eeg_amps$rhythm <- ifelse(grepl('_unsyncopated',eeg_amps$condition),'unsyncopated','syncopated')
eeg_amps$rhythm <- factor(eeg_amps$rhythm, levels=c('unsyncopated','syncopated'))

# get mean meter-related frequencies
eeg_amps_meter <- subset(eeg_amps, isMeterRel==1)
eeg_amps_meter <- summaryBy(zscore ~ ID+rhythm+tone, FUN=mean, data=eeg_amps_meter)

# fit a mixed model
m1_eeg_amps_meter <- lmer(zscore.mean ~ rhythm * tone + (1|ID), data=eeg_amps_meter)
Anova(m1_eeg_amps_meter,test="F")

# decompose the interaction
posthoc_eeg_amps_meter <- lmer(zscore.mean ~ rhythm + rhythm:tone -1 + (1|ID), data=eeg_amps_meter)
summary(glht(posthoc_eeg_amps_meter, linfct=c('rhythmunsyncopated:tonelow==0', 
                                              'rhythmsyncopated:tonelow==0')), test=adjusted('bonferroni'))
