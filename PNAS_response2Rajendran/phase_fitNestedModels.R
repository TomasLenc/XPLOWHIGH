


# check if required packages are installed
list.of.packages <- c("dplyr", "bpnreg", "pander")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load packages
library(bpnreg)
library(pander)
library(dplyr)




# load data
eeg_phase <- read.csv('./data/XPLOWHIGH_phases.csv')

eeg_phase$ID <- as.double(eeg_phase$ID)

eeg_phase$tone <- ifelse(grepl('^H',eeg_phase$condition),'high','low')
eeg_phase$tone <- factor(eeg_phase$tone, levels=c('high','low'))

eeg_phase$rhythm <- ifelse(grepl('_unsyncopated',eeg_phase$condition),'unsyncopated','syncopated')
eeg_phase$rhythm <- factor(eeg_phase$rhythm, levels=c('unsyncopated','syncopated'))

eeg_phase$freq <- as.factor(eeg_phase$freq)

eeg_phase <- as_tibble(eeg_phase)



# fit the models (!!! this takes a long time to run)
fit_eegPhase_interceptOnly <- bpnme(pred.I = phase ~ (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)

fit_eegPhase_mainFreq <- bpnme(pred.I = phase ~ freq + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_interceptOnly$model.fit[1,'WAIC']
fit_eegPhase_mainFreq$model.fit[1,'WAIC']

fit_eegPhase_mainFreq_mainRhythm <- bpnme(pred.I = phase ~ freq + rhythm + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_mainFreq$model.fit[1,'WAIC']
fit_eegPhase_mainFreq_mainRhythm$model.fit[1,'WAIC']

fit_eegPhase_mainFreq_mainRhythm_mainTone <- bpnme(pred.I = phase ~ freq + rhythm + tone + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_mainFreq_mainRhythm$model.fit[1,'WAIC']
fit_eegPhase_mainFreq_mainRhythm_mainTone$model.fit[1,'WAIC']


fit_eegPhase_FREQxRHYTHM <- bpnme(pred.I = phase ~ freq + rhythm + freq:rhythm + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_mainFreq_mainRhythm$model.fit[1,'WAIC']
fit_eegPhase_FREQxRHYTHM$model.fit[1,'WAIC']

fit_eegPhase_FREQxRHYTHM_FREQxTONE <- bpnme(pred.I = phase ~ freq + rhythm + freq:rhythm + freq:tone + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_FREQxRHYTHM$model.fit[1,'WAIC']
fit_eegPhase_FREQxRHYTHM_FREQxTONE$model.fit[1,'WAIC']

fit_eegPhase_FREQxRHYTHM_RHYTHMxTONE <- bpnme(pred.I = phase ~ freq + rhythm + freq:rhythm + rhythm:tone + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_FREQxRHYTHM$model.fit[1,'WAIC']
fit_eegPhase_FREQxRHYTHM_RHYTHMxTONE$model.fit[1,'WAIC']

fit_eegPhase_FULL <- bpnme(pred.I = phase ~ freq * rhythm * tone + (1|ID), data = eeg_phase, its = 10000, burn = 1000, n.lag = 3, seed = 101)
fit_eegPhase_FREQxRHYTHM$model.fit[1,'WAIC']
fit_eegPhase_FULL$model.fit[1,'WAIC']


# print a nice table
model_res <- matrix(c(fit_eegPhase_interceptOnly$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_mainFreq$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_mainFreq_mainRhythm$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_mainFreq_mainRhythm_mainTone$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_FREQxRHYTHM$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_FREQxRHYTHM_FREQxTONE$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_FREQxRHYTHM_RHYTHMxTONE$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")], 
                      fit_eegPhase_FULL$model.fit[1, c("DIC","DICalt","WAIC","WAIC2")]), 
                      nrow=4, byrow=F, dimnames=list(c("DIC","DICalt","WAIC","WAIC2"),c('intercept-only', 'frequency',
                                                                                      'frequency + rhythm', 
                                                                                      'frequency + rhythm + tone', 
                                                                                      'frequency + rhythm + frequency:rhythm', 
                                                                                      'frequency + rhythm + frequency:rhythm + frequency:tone', 
                                                                                      'frequency + rhythm + frequency:rhythm + rhythm:tone', 
                                                                                      'full model'
                      )))
panderOptions('round', 2)
pander(pandoc.table(model_res, split.cells = rep(1,dim(model_res)[2]), split.table=Inf))
panderOptions('round', 3)
