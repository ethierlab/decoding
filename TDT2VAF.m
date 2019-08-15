%% TDT2mat_et_waveculs_fonction.m
% Directory of TDT block location
BLOCKPATH = uigetdir('/Users/christianethier/Google Drive (Work)/Projects/Chronic Array and Mototrak/sample data/jados-19-04-30-am','Select TDT data block');

% Extracts data from all channels in 'spik' storage
spikes = TDTbin2mat(BLOCKPATH, 'TYPE', {'all'});

% Total number of channels
chan_size = size(spikes.streams.spik.data,1);

%Force stuff
force_data = spikes.streams.Fore.data; %data
fs{1,1} = spikes.streams.Fore.fs; % sampling rate

%LFP stuff
lfp_data = spikes.streams.Lfp1.data;
fs{1,3} = spikes.streams.Lfp1.fs;

%spik stuff
fs{1,4} = spikes.streams.spik.fs;

%% trial_extraction_main_function.m + trig2event.m
hold_time = 0.8;

events = trig2event(spikes.epocs.Bev_.onset);

session_time = [spikes.info.date(end-1:end) '-' spikes.info.date(end-5:end-3) '-' spikes.info.date(1:4) ' ' spikes.info.utcStartTime];
serial_session_time = datenum(session_time);


%% extraction_signals_during_trials2

k = 1;
for i = 1:size(events,1) %i: numéro d'événement
    if strfind(events{i,2},'start')
        trial_start_time(k) = events{i,1};
        if i+1 <= size(events,1)
            if strfind(events{i+1,2},'end') %si pas de trig reward (essai échoué)
                trial_end_time(k) = events{i+1,1};
                success(k) = 0;
                time2reward(k) = nan;
            else
                if i+2 <= size(events,1)
                    if strfind(events{i+2,2},'end') %si trig reward (essai réussi)
                        trial_end_time(k) = events{i+1,1};
                        time2reward(k) = (trial_end_time(k)-trial_start_time(k))-hold_time; %time2reward
                        %extrayons jusqu'à 500 ms après délivrance récompense (pas artefacts mastication)
                        trial_end_time(k) = trial_end_time(k)+.5;
                        success(k) = 1;
                    end
                else
                    trial_start_time(k) = [];
                end
            end
        else
            trial_start_time(k) = [];
        end
        
        if numel(trial_start_time) == k %ce start_time n'a pas été supprimé à défaut de end_time correspondant
            if strfind(events{i,2},'jackpot')
                jackpot(k) = 1;
            elseif strfind(events{i,2},'normal')
                jackpot(k) = 0;
            elseif strfind(events{i,2},'noreward')
                jackpot(k) = -1;
            end
        end
        k = k+1;
    end 
end
        
for i = 1:numel(trial_start_time)
    if trial_start_time(i) - 0.5 >= 0
        trial_start_time(i) = trial_start_time(i) - 0.5; 	
    end
end

for i = 1:numel(trial_start_time)
    data1_SSDP{1,i} = jackpot(i);
    data1_SSDP{2,i} = success(i);
    data1_SSDP{3,i} = time2reward(i);
end

trial_names = {};
for i = 1:numel(trial_start_time) %i: numéro d'essai de l'expérience
    data1_SSDP{4,i} = force_data(round(trial_start_time(i)*fs{1,1}):round(trial_end_time(i)*fs{1,1}));
    trial_names = [trial_names ['trial' num2str(i)]];
end

for i = 1:numel(trial_start_time) %i: numéro d'essai de l'expérience
    for j = 1:chan_size %numéro d'électrode
        data1_SSDP{7+j,i} = lfp_data(j,round(trial_start_time(i)*fs{1,3}):round(trial_end_time(i)*fs{1,3}));
    end
end

for v = 1:size(lfp_data,1)
    lfp_channels{v}=['LFP' num2str(v)];
end

cd([BLOCKPATH '\chandata'])

index = 8 + chan_size;
for i = 1:chan_size
    load_name = ['ch' num2str(i) '_snips.mat'];
    load (load_name, '-regexp', '^(?!fs$).')
    
    for n = 1:numel(trial_start_time)
    data1_SSDP{index,n} = intersect(...
        (round (trial_start_time(n)*fs{1,4}):round(trial_end_time(n)*fs{1,4})), ...
        round(ts*fs{1,4})) - (round(trial_start_time(n)*fs{1,4})-1);
    end
end
    
row = [{'trial_type','successful/unsuccessful','time2succes',...
        'FORCE','EMG1','EMG2','EMG3'},lfp_channels,'ch_mu'];   
data1_SSDP = cell2table(data1_SSDP, 'RowNames', row, 'VariableNames', trial_names);

col2={'FORCE','EMG','LFP','spikes'};
fs = cell2table(fs, 'VariableNames', col2);
%% remove_false_spikes_ce.m
% This may be unnecessary since we are only working with single clusters
% for each

window_size = 0.5;

spike_i = data1_SSDP(index,:);

numbins_rem = ceil(window_size*10^-3*fs.spikes);

%[4] minimum number of channels (or clusters) on which spikes have to occur
% simultaneously in order to be considered noise and removed:
numchan_rem = 4;

data1_spikes = data1_SSDP(8 + chan_size,:);
num_trials = length(trial_names);

rem_spikes = cell(num_trials,1);

data2_SSDP = data1_SSDP;

% Run artifact detection on a trial per trial basis:

for trial = 1:num_trials
    all_spikes = vertcat(data1_spikes{1,trial}{:});
    
    if isempty(all_spikes)
        %no spikes this trial, skip to next
        continue;
    end
    
    %define time vector spanning all relevant bins
    edges = min(all_spikes):max(all_spikes)+1;
    
    %count number of spikes per time bin:
    N = histcounts(all_spikes,edges);
    
    %sum number of spikes over specified time window
    N = movsum(N,numbins_rem);
    
 % find time indices 'numbins' around where spike sum is >= number of chans
    edges = edges(1:end-1); %last edge was added just for histcounts to work nicely
    invalid_bins = edges(logical(movsum(N>=numchan_rem,numbins_rem)));
    
    rem_spikes{trial} = all_spikes(ismember(all_spikes,invalid_bins));
    
    if ~isempty(invalid_bins)
        % remove invalid spikes from data, channel by channel:
        for chan = 1:num_chans
            spikes = data1_SSDP{data1_spikes(chan),trial}{:};
            data2_SSDP{data1_spikes(chan),trial}{:} = setdiff(spikes,invalid_bins);
        end
    end
end

%% Convert to binnedData and calculate VAF

for var = 1: size(data2_SSDP,2)
is_empty(var,1) = isempty(data2_SSDP{4,var}{:});
end

if any(is_empty ~= 0)
data2_SSDP = removevars(data2_SSDP,data2_SSDP.Properties.VariableNames{(is_empty)});
end

if max(data2_SSDP{4,1}{:}) < 20
    data2_SSDP = force_analog2grams(data2_SSDP,fs);
end

% Create binnedData
bin_size = 0.05; % seconds; set this bin size used during data acquisition
params.binsize = bin_size;

binnedData = extraire_bins_pour_force_et_clusters(data2_SSDP,fs,'tout','initiation',params);

force_pred_mototrak(binnedData);

%% Clear vars
clearvars -except data2_SSDP binnedData fs hold_time serial_session_time session_time mean_vaf stdev_vaf vaf_per_fold act_force pred_force



