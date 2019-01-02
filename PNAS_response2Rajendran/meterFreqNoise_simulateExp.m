

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% This code shows that decreasing coherence of the response at meter-related frequencies decreases ITPC at the beat frequency, 
% stimulus-response phase locking at the beat frequency, but also meter-related z-scored amplitudes as revealed by frequency-tagging. 
% 
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all



% ------------- define parameters -------------

fs      = 500; 
pattern = [1 1 1 0 1 1 1 0 1 1 0 0]; 
IOI     = 0.2; 
rampon  = 0.01; 
rampoff = 0.05; 

n_cycles        = 24; 
n_beats         = n_cycles*3; 
beat_times      = IOI*4 * [1:n_beats]; 
pattern_times   = IOI*12 * [1:n_cycles]; 

n_experiments   = 50; 
n_partic        = 20; 
n_trials        = 10; 


frex = 1/(length(pattern)*IOI) * [1:12]; % 12 frequencies of interest
frex2jitter = 1/(length(pattern)*IOI) * [1,3,6,12]; % frequencies to jitter (meter-related frequencies)

different_phases_between_cond_BOOL = 1; % for each participant, generate random phases for the response separately for the standard and jittered condition

% frequency-jittering parameters
freq_jitter_mean = 0; % mean
freq_jitter_sd = [0.05, 0.1, 0.2]; % standard deviation

% slding window fft parameters
fft_win_overlap = 0.5; % 50% overlap
fft_win_dur = IOI*length(pattern)*4;  % 4 pattern cycles
fft_n_win = round(((n_cycles*length(pattern)*IOI) - (fft_win_dur*fft_win_overlap)) / (fft_win_dur*fft_win_overlap)); 




% ------------- create the stimulus and calculate its spectra -------------

t = [0 : round(n_cycles*length(pattern)*IOI*fs)-1]/fs; 
env_stimulus = zeros(1, length(t)); 
env_event = ones(1,round(IOI*fs)); 
env_event(1:round(rampon*fs)) = linspace(0,1,round(rampon*fs)); 
env_event(end-round(rampoff*fs)+1:end) = linspace(1,0,round(rampoff*fs)); 
for i=0:n_cycles*length(pattern)-1
    if pattern(mod(i,length(pattern))+1)        
        env_stimulus(round(i*IOI*fs)+1:round(i*IOI*fs+length(env_event))) = env_event; 
    end
end
N               = length(env_stimulus); 
hN              = floor(N/2)+1; 
X               = fft(env_stimulus);
Xshifted        = fftshift(X); 
idx5Hz          = round(5/fs*N)+1; 
Xshifted([1:hN-idx5Hz-1,hN+idx5Hz:end]) = 0; 
mX_original_s   = abs(Xshifted); 


% calculate phase angles with a sliding window for the stimulus at 12 frequencies of interest
N_win = fft_win_dur*fs; 
hN_win = floor(N_win/2)+1; 
fft_win_frex_idx = round(frex*N_win/fs)+1; 

aX_stim_slidingWin = zeros(fft_n_win, length(frex)); 
for wini=1:fft_n_win
    idx = (wini-1)*round(fft_win_dur*fft_win_overlap*fs); 
    x = env_stimulus(idx+1:idx+round(fft_win_dur*fs)); 
    aX = angle(fft(x)); 
    aX_stim_slidingWin(wini,:) = aX(fft_win_frex_idx); 
end






% ------------- allocate variables for the simulation -------------

mX_res = zeros(length(freq_jitter_sd),n_experiments,n_partic,2,N); 

H_zscore = nan(length(freq_jitter_sd),n_experiments); 
P_zscore = nan(length(freq_jitter_sd),n_experiments); 

H_phase_coherence = nan(length(freq_jitter_sd),n_experiments); 
P_phase_coherence = nan(length(freq_jitter_sd),n_experiments); 



% ------------- run the simulation -------------

if isempty(gcp('nocreate')); parpool(3); end; % use parallel computation to speed up the simulation

n_jitters = length(freq_jitter_sd); 
parfor jitteri=1:n_jitters
        disp(sprintf('\nfrequency jitter SD %d out of %d\n\n',jitteri,length(freq_jitter_sd)))
        
        
    for expi=1:n_experiments
        disp(sprintf('epxeriment %d',expi))

        % allocate variables for the current experiment
        response = zeros(n_partic, n_trials, 2, N); % simulated responses in the time domain
        phase_locking = zeros(n_partic, n_trials, 2, length(frex)); % phase-locking between stimulus and response for each participant, condition and trial
        
        
        for partici=1:n_partic

            % generate random phases for each participant (and condition if
            % different_phases_between_cond_BOOL == 1)
            phases_half = (rand(1,N/2-1)*2*pi-pi); 
            phases_s = [0, fliplr(-phases_half), 0, phases_half]; 
            if different_phases_between_cond_BOOL
                phases_half = (rand(1,N/2-1)*2*pi-pi); 
            end
            phases_s_jittered = [0, fliplr(-phases_half), 0, phases_half]; 

            
            for triali=1:n_trials

                % generate response in the standard condition
                frex_idx = find(mX_original_s); 
                frex_idx = frex_idx(13:end); 

                frex2use = 1/(length(pattern)*IOI) * [0:12]; 
                amps2use = mX_original_s(frex_idx); 
                phases2use_s = phases_s(frex_idx); 
                phases2use_s_jittered = phases_s_jittered(frex_idx); 

                resp_standard = zeros(size(t)); 
                for fi=1:length(amps2use)
                    resp_standard = resp_standard + amps2use(fi)*cos(2*pi*t*frex2use(fi) + phases2use_s(fi)); 
                end

                % generate response in the jittered condition
                resp_jittered = zeros(size(t)); 
                for fi=1:length(amps2use)        
                    if ismember(frex2use(fi),frex2jitter)
                        f_t = repmat(frex2use(fi), 1, N) + (randn(1,N)+freq_jitter_mean)*freq_jitter_sd(jitteri); 
                        resp_jittered = resp_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 
                    else
                        f_t = repmat(frex2use(fi), 1, N); 
                        resp_jittered = resp_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 

                    end
                end
                
                response(partici,triali,1,:) = resp_standard; 
                response(partici,triali,2,:) = resp_jittered; 

                
                % calculate stimulus-response phase locking for this trial
                aX_resp_standard_slidingWin = zeros(fft_n_win, length(frex)); 
                for wini=1:fft_n_win
                    idx = (wini-1)*round(fft_win_dur*fft_win_overlap*fs); 
                    x = resp_standard(idx+1:idx+round(fft_win_dur*fs)); 
                    aX = angle(fft(x)); 
                    aX_resp_standard_slidingWin(wini,:) = aX(fft_win_frex_idx); 
                end
                aX_resp_jittered_slidingWin = zeros(fft_n_win, length(frex)); 
                for wini=1:fft_n_win
                    idx = (wini-1)*round(fft_win_dur*fft_win_overlap*fs); 
                    x = resp_jittered(idx+1:idx+round(fft_win_dur*fs));                     
                    aX = angle(fft(x)); 
                    aX_resp_jittered_slidingWin(wini,:) = aX(fft_win_frex_idx); 
                end
                plv_standard = abs(mean(exp(1i*(aX_resp_standard_slidingWin-aX_stim_slidingWin)),1)); 
                plv_jittered = abs(mean(exp(1i*(aX_resp_jittered_slidingWin-aX_stim_slidingWin)),1)); 
                
                phase_locking(partici,triali,1,:) = plv_standard; 
                phase_locking(partici,triali,2,:) = plv_jittered; 
                
            end

        end
        
        
        % calculate meter zscores and do the ttest
        time_avg = squeeze(mean(response,2)); 

        frex_idx = round(frex*N/fs)+1; 
        idx_meterRel = [1,3,6,12]; 

        mX_simulated = abs(fft(time_avg,[],3))*2/N; 
        mX_res(jitteri,expi,:,:,:) = mX_simulated; 

        amps = mX_simulated(:,:,frex_idx); 
        z = zscore(amps,[],3); 

        z_meterRel_s = mean(z(:,1,idx_meterRel),3); 
        z_meterRel_s_jittered = mean(z(:,2,idx_meterRel),3); 

        meter_zscore_diff = z_meterRel_s - z_meterRel_s_jittered; 
        [H_zscore(jitteri,expi),P_zscore(jitteri,expi)] = ttest(z_meterRel_s,z_meterRel_s_jittered); 

         
        % do the ttest on phase-locking values
        phase_locking_mean = mean(phase_locking,2); 
        phase_locking_beat = squeeze(phase_locking_mean(:,:,:,3)); 
        [H_phase_coherence(jitteri,expi),P_phase_coherence(jitteri,expi)] = ttest(phase_locking_beat(:,1),phase_locking_beat(:,2)); 
        
    end

end







%% RESULTS

for jitteri=1:length(freq_jitter_sd)
    fprintf('\n-------------------------\nFrequency jitter = %.2f\n',freq_jitter_sd(jitteri))
    fprintf('Proportion of experiments with significant meter zscore difference = %.2f\n',  sum(H_zscore(jitteri,:))/n_experiments)
    fprintf('Proportion of experiments with significant stimulus-response coherence difference = %.2f\n',  sum(H_phase_coherence(jitteri,:))/n_experiments)
end














%% PLOT 2 EXAMPLE TRIALS


freq_jitter_sd2plot = 0.1; 

phases_half = (rand(1,N/2-1)*2*pi-pi); 
phases_s = [0, fliplr(-phases_half), 0, phases_half]; 
if different_phases_between_cond_BOOL
    phases_half = (rand(1,N/2-1)*2*pi-pi); 
end
phases_s_jittered = [0, fliplr(-phases_half), 0, phases_half]; 
for triali=1:n_trials

    % generate standard version
    frex_idx = find(mX_original_s); 
    frex_idx = frex_idx(13:end); 

    frex2use = 1/(length(pattern)*IOI) * [0:12]; 
    amps2use = mX_original_s(frex_idx); 
    phases2use_s = phases_s(frex_idx); 
    phases2use_s_jittered = phases_s_jittered(frex_idx); 

    resp_standard = zeros(size(t)); 
    for fi=1:length(amps2use)
        resp_standard = resp_standard + amps2use(fi)*cos(2*pi*t*frex2use(fi) + phases2use_s(fi)); 
    end

    % generate jittered version
    resp_jittered = zeros(size(t)); 
    for fi=1:length(amps2use)        
        if ismember(frex2use(fi),frex2jitter)
            f_t = repmat(frex2use(fi), 1, N) + (randn(1,N)+freq_jitter_mean)*freq_jitter_sd2plot; 
            resp_jittered = resp_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 
        else
            f_t = repmat(frex2use(fi), 1, N); 
            resp_jittered = resp_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 

        end
    end
end

% plot the original and jittered version 
figure('color','white')
subplot 211
plot(t,resp_standard,'k','LineWidth',1.7); 
hold on
plot(beat_times,max(resp_standard)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,(length(pattern)*IOI)*3])
box off
title('original')
subplot 212
plot(t,resp_jittered,'LineWidth',1.7,'color',[0.4902,0.1804,0.5608]); 
hold on
plot(beat_times,max(resp_jittered)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,(length(pattern)*IOI)*3])
box off
title('jittered frequency')
xlabel('time (s)')

% plot FFT spectra of the waveforms
figure('name', 'FFT of the two waveforms')
subplot 211
mXtmp = abs(fft(resp_standard)); 
plot([0:hN-1]/N*fs,mXtmp(1:hN),'b','LineWidth',1.7); 
xlim([0.1,6])
title('original')
subplot 212
mXtmp = abs(fft(resp_jittered)); 
plot([0:hN-1]/N*fs,mXtmp(1:hN),'r','LineWidth',1.7); 
xlim([0.1,6])
title('jittered frequency')
xlabel('frequency (Hz)'); 

        
    
    
    
    
