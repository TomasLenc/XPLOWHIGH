

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% This code shows that decreasing coherence of the response at meter-related frequencies decreases ITPC at the beat frequency, 
% but also meter-related z-scored amplitudes as revealed by frequency-tagging. 
% 
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all

fs = 500; 

pattern = [1 1 1 0 1 1 1 0 1 1 0 0]; 

IOI = 0.2; 
rampon = 0.01; 
rampoff = 0.05; 

ncycles = 24; 


% create the stimulus 
t = [0 : round(ncycles*length(pattern)*IOI*fs)-1]/fs; 
env = zeros(1, length(t)); 
env_event = ones(1,round(IOI*fs)); 
env_event(1:round(rampon*fs)) = linspace(0,1,round(rampon*fs)); 
env_event(end-round(rampoff*fs)+1:end) = linspace(1,0,round(rampoff*fs)); 
for i=0:ncycles*length(pattern)-1
    if pattern(mod(i,length(pattern))+1)        
        env(round(i*IOI*fs)+1:round(i*IOI*fs+length(env_event))) = env_event; 
    end
end
N = length(env); 
hN = floor(N/2)+1; 
X = fft(env);
Xshifted = fftshift(X); 
idx5Hz = round(5/fs*N)+1; 
Xshifted([1:hN-idx5Hz-1,hN+idx5Hz:end]) = 0; 
mX_original_s = abs(Xshifted); 

n_beats = length(t)/fs/(IOI*4); 
beat_times = IOI*4 * [1:n_beats]; 
pattern_times = IOI*12 * [1:ncycles]; 






%% SIMULATION



% ------------- DEFINE -------------

different_phases_between_cond_BOOL = 1; % for each participant, generate random phases for the response separately for the standard and jittered condition

n_experiments = 50; 
n_partic = 20; 
n_trials = 10; 

frex2jitter = 1/2.4 * [1,3,6,12]; % meter-related frequencies

% frequency-jittering parameters
freq_jitter_mean = 0; 
freq_jitter_sd = [0.05, 0.1, 0.2]; 




% ------------- ALLOCATE -------------

mX_res = zeros(length(freq_jitter_sd),n_experiments,n_partic,2,N); 

H_ITPC = nan(length(freq_jitter_sd),n_experiments); 
P_ITPC = nan(length(freq_jitter_sd),n_experiments); 

H_zscore = nan(length(freq_jitter_sd),n_experiments); 
P_zscore = nan(length(freq_jitter_sd),n_experiments); 






% ------------- RUN -------------

if isempty(gcp('nocreate')); parpool(3); end; 
jitter_length = length(freq_jitter_sd); 
parfor jitteri=1:jitter_length
        disp(sprintf('\nfrequency jitter SD %d out of %d\n\n',jitteri,length(freq_jitter_sd)))
        
    for expi=1:n_experiments
        disp(sprintf('epxeriment %d',expi))

        res = zeros(n_partic, n_trials, 2, N); 
        for partic=1:n_partic

            phases_half = (rand(1,N/2-1)*2*pi-pi); 
            phases_s = [0, fliplr(-phases_half), 0, phases_half]; 
            if different_phases_between_cond_BOOL
                phases_half = (rand(1,N/2-1)*2*pi-pi); 
            end
            phases_s_jittered = [0, fliplr(-phases_half), 0, phases_half]; 

            for trial=1:n_trials

                % generate standard version
                frex_idx = find(mX_original_s); 
                frex_idx = frex_idx(13:end); 

                frex2use = 1/2.4 * [0:12]; 
                amps2use = mX_original_s(frex_idx); 
                phases2use_s = phases_s(frex_idx); 
                phases2use_s_jittered = phases_s_jittered(frex_idx); 

                s = zeros(size(t)); 
                for fi=1:length(amps2use)
                    s = s + amps2use(fi)*cos(2*pi*t*frex2use(fi) + phases2use_s(fi)); 
                end

                % generate jittered version
                s_jittered = zeros(size(t)); 
                for fi=1:length(amps2use)        
                    if ismember(frex2use(fi),frex2jitter)
                        f_t = repmat(frex2use(fi), 1, N) + (randn(1,N)+freq_jitter_mean)*freq_jitter_sd(jitteri); 
                        s_jittered = s_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 
                    else
                        f_t = repmat(frex2use(fi), 1, N); 
                        s_jittered = s_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 

                    end
                end
                
                res(partic,trial,1,:) = s; 
                res(partic,trial,2,:) = s_jittered; 

            end

        end
        
        % calculate meter zscores
        time_avg = squeeze(mean(res,2)); 

        frex = 1/2.4*[1:12]; 
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

        
        % calculate ITPC
        aX = angle(fft(res,[],4)); 
        beat_angles = aX(:,:,:,frex_idx(3)); 
        ITPC = abs(mean(exp(1i*beat_angles),2)); 
        [H_ITPC(jitteri,expi),P_ITPC(jitteri,expi)] = ttest(ITPC(:,:,1),ITPC(:,:,2)); 

    end

end





%% RESULTS

for jitteri=1:length(freq_jitter_sd)
    fprintf('\n-------------------------\nFrequency jitter = %.2f\n',freq_jitter_sd(jitteri))
    fprintf('Proportion of experiments with significant meter zscore difference = %.2f\n',  sum(H_zscore(jitteri,:))/n_experiments)
    fprintf('Proportion of experiments with significant ITPC difference = %.2f\n',  sum(H_ITPC(jitteri,:))/n_experiments)
end











%% PLOT 2 EXAMPLE TRIALS


freq_jitter_sd2plot = 0.1; 

phases_half = (rand(1,N/2-1)*2*pi-pi); 
phases_s = [0, fliplr(-phases_half), 0, phases_half]; 
if different_phases_between_cond_BOOL
    phases_half = (rand(1,N/2-1)*2*pi-pi); 
end
phases_s_jittered = [0, fliplr(-phases_half), 0, phases_half]; 
for trial=1:n_trials

    % generate standard version
    frex_idx = find(mX_original_s); 
    frex_idx = frex_idx(13:end); 

    frex2use = 1/2.4 * [0:12]; 
    amps2use = mX_original_s(frex_idx); 
    phases2use_s = phases_s(frex_idx); 
    phases2use_s_jittered = phases_s_jittered(frex_idx); 

    s = zeros(size(t)); 
    for fi=1:length(amps2use)
        s = s + amps2use(fi)*cos(2*pi*t*frex2use(fi) + phases2use_s(fi)); 
    end

    % generate jittered version
    s_jittered = zeros(size(t)); 
    for fi=1:length(amps2use)        
        if ismember(frex2use(fi),frex2jitter)
            f_t = repmat(frex2use(fi), 1, N) + (randn(1,N)+freq_jitter_mean)*freq_jitter_sd2plot; 
            s_jittered = s_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 
        else
            f_t = repmat(frex2use(fi), 1, N); 
            s_jittered = s_jittered + amps2use(fi)*cos(2*pi*cumsum(f_t)/fs+phases2use_s_jittered(fi)); 

        end
    end
end

% plot the original and jittered version 
figure('color','white')
subplot 211
plot(t,s,'k','LineWidth',1.7); 
hold on
plot(beat_times,max(s)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,2.4*3])
box off
title('original')
subplot 212
plot(t,s_jittered,'LineWidth',1.7,'color',[0.4902,0.1804,0.5608]); 
hold on
plot(beat_times,max(s_jittered)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,2.4*3])
box off
title('jittered frequency')
xlabel('time (s)')

% plot FFT spectra of the waveforms
figure('name', 'FFT of the two waveforms')
subplot 211
mXtmp = abs(fft(s)); 
plot([0:hN-1]/N*fs,mXtmp(1:hN),'b','LineWidth',1.7); 
xlim([0.1,6])
title('original')
subplot 212
mXtmp = abs(fft(s_jittered)); 
plot([0:hN-1]/N*fs,mXtmp(1:hN),'r','LineWidth',1.7); 
xlim([0.1,6])
title('jittered frequency')
xlabel('frequency (Hz)'); 

        
    
    
    
    
