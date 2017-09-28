function binnedData = convertMDF2binned(datatable,varargin)
%
% usage:  binnedData = convertMDF2binned(datastruct,[params])
%
% converts a "Michael Data Format" file to the binned format, according to parameters
% specified in the optional [params] argument.
%
%         datatable               : string of mdf.mat file path and name, or string of variable name in workspace, or actual data table
%
%         params fields:            [default values in brackets]
%                                   none, one or many of these fields can
%                                   be provided in the params argument
%                                   structure, any missing field will be
%                                   set to its default value.
%
%             binsize             : [0.05]   desired bin size in 
%             pre_capture         : [0.5]    duration of pre-trial recording
%             HP, Lp              : [50 10]  high pass and low pass cut off frequencies for EMG filtering
%             NormData            : [false]  specify whether the output data is to be normalized to unity
%             ArtRemEnable        : [false]  Whether or not to attempt detecting and deleting artifacts
%             NumChan             : [10]     Number of channels from which the artifact removal needs to detect simultaneous spikes to consider it an artifact
%             TimeWind            : [0.0005] time window, in seconds, over which the artifact remover will consider event to be "simultaneous"
%
%%%% Ethierlab 2017/09/14 -- CE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~istable(datatable)
    %Load the file or structure
    [datatable, fs] = load_data_table(datatable);
    if isempty(datatable)
        error('can''t load file');
    end
else
    %use the data table already in workspace, load fs too
    fs = evalin('base','fs');
end

%update missing params with default values
params = get_default_binning_params(varargin{:});

if isempty(params)
    disp('Invalid binning parameter(s)');
    return
end

%% Get general data table information
row_names    = datatable.Properties.RowNames;
num_rows     = length(row_names);
trial_labels = datatable.Properties.VariableNames;
num_trials   = length(trial_labels);

% get row indexes of every data type
EMG_i    = find(strncmpi('EMG',row_names,3));
LFP_i    = find(strncmpi('LFP',row_names,3));
force_i  = find(strcmpi('force',row_names));
spike_i  = find(strncmpi('clust',row_names,5));
jackpot_match = strfind(row_names,'jackpot');
jackpot_i = [];
for i = 1:num_rows
    if ~isempty(jackpot_match{i})
        jackpot_i = [jackpot_i; i];
    end
end

% get sampling frequency for each data type
fs_names = fs.Properties.VariableNames;
EMG_fs   = fs{1,strncmpi('EMG',fs_names,3)};
force_fs = fs{1,strncmpi('force',fs_names,5)};
LFP_fs   = fs{1,strncmpi('LFP',fs_names,3)};
spike_fs = fs{1,strncmpi('spike',fs_names,5)};

%% Initialize bin variables
    num_emgs    = length(EMG_i);
    num_force   = length(force_fs);
    num_spikes  = length(spike_i);
    num_LFP     = length(LFP_i);
    
    emg_bin     = cell(num_trials,num_emgs);
    force_bin   = cell(num_trials,num_force);
    spike_bin   = cell(num_trials,num_spikes);
    LFP_bin     = cell(num_trials,num_LFP);
    beta_bin    = cell(num_trials,num_LFP);
    lgamma_bin  = cell(num_trials,num_LFP);
    hgamma_bin  = cell(num_trials,num_LFP);
    vhgamma_bin = cell(num_trials,num_LFP);
    jackpot_bin = cell(num_trials,1);
    
    trial_dur   = cell(num_trials,1);
    timeframe   = cell(num_trials,1);

    LFPs = process_LFPs(datatable{LFP_i,:},LFP_fs);
    
    
%%  Trial per trial, extract and bin all data into new cells

for trial = 1:num_trials
    
    %use first force signal to infer duration and timeframe for this trial
    num_points       = numel(datatable{force_i(1),trial}{:});
    trial_dur{trial} = num_points/force_fs;
    num_bins         = floor(trial_dur{trial}/params.binsize);
    timeframe{trial} = params.binsize*(0:num_bins-1);
    
    %% 1-Bin EMG data
    if isempty(EMG_i)
        disp('No EMG data was found');
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMG data is hi-pass filtered at 50Hz, rectified and low pass filtered
        % at 10Hz, unless otherwise specified. It is downsampled to match the desired binsize.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % Filter for EMG data
        [bh,ah] = butter(4, params.EMG_hp*2/EMG_fs, 'high'); %highpass filter params
        [bl,al] = butter(4, params.EMG_lp*2/EMG_fs, 'low');  %lowpass filter params
        
        % original timeframe
        emgtimebins = (0:numel(datatable{EMG_i(1),trial}{:})-1)/EMG_fs;
        for E=1:num_emgs
            % Filter EMG data
            tempEMG = datatable{EMG_i(E),trial}{:};
            tempEMG = filtfilt(bh,ah,tempEMG); %highpass filter
            tempEMG = abs(tempEMG); %rectify
            tempEMG = filtfilt(bl,al,tempEMG); %lowpass filter
            %downsample EMG data to desired bin size
            emg_bin{trial,E} = interp1(emgtimebins, tempEMG, timeframe{trial},'linear','extrap')';
%             %Normalize EMGs
%             if params.NormData             
%                 %normalize to 99% percentile. Dont use the max because of possible large artefacts
%                 EMGNormRatio     = prctile(emg_bin{trial,E},99);
%                 emg_bin{trial,E} = emg_bin{trial,E}/EMGNormRatio;
%             end
        end
    end
    clear tempEMG emgtimebins E bh ah bl al EMGNormRatio;
    
    %% 2-Bin force data    
    if isempty(force_i)
        disp('No force data was found');
    else   
        forcetimebins = (0:numel(datatable{force_i(1),trial}{:})-1)/force_fs;
        for F = 1:num_force
            tempForce = datatable{force_i(F),trial}{:};
            force_bin{trial,F} = interp1(forcetimebins, tempForce, timeframe{trial},'linear','extrap')';
%             %Normalize Force
%             if params.NormData            
%                 %dont use the max because of possible outliars, use 99% percentile
%                 forceNormRatio     = prctile(force_bin{trial,F},99);
%                 force_bin{trial,F} = force_bin{trial,F}/forceNormRatio;
%             end
        end
    end
    clear tempForce forcetimebins F forceNormRatio;
        
    %% 3.1-Bin LFP data    
    if isempty(LFP_i)
        disp('No LFP data was found');
    else   
        LFPtimebins = (0:numel(datatable{LFP_i(1),trial}{:})-1)/LFP_fs;
        for L = 1:num_LFP
            LFP_bin{trial,L}     = interp1(LFPtimebins, datatable{LFP_i(L),trial}{:}, timeframe{trial},'linear','extrap')';
            beta_bin{trial,L}    = interp1(LFPtimebins, LFPs.beta{L,trial}          , timeframe{trial},'linear','extrap')';
            lgamma_bin{trial,L}  = interp1(LFPtimebins, LFPs.lgamma{L,trial}        , timeframe{trial},'linear','extrap')';
            hgamma_bin{trial,L}  = interp1(LFPtimebins, LFPs.hgamma{L,trial}        , timeframe{trial},'linear','extrap')';
            vhgamma_bin{trial,L} = interp1(LFPtimebins, LFPs.vhgamma{L,trial}       , timeframe{trial},'linear','extrap')';
%             %Normalize LFP
%             if params.NormData            
%                 %dont use the max because of possible outliars, use 99% percentile
%                 LFP_bin{trial,L} = LFP_bin{trial,L}/prctile(LFP_bin{trial,L},99);
%                 beta_bin{trial,L} = beta_bin{trial,L}/prctile(beta_bin{trial,L},99);
%                 LFP_bin{trial,L} = LFP_bin{trial,L}/prctile(LFP_bin{trial,L},99);
%                 LFP_bin{trial,L} = LFP_bin{trial,L}/prctile(LFP_bin{trial,L},99);
%                 LFP_bin{trial,L} = LFP_bin{trial,L}/prctile(LFP_bin{trial,L},99);
%             end
        end
    end
    clear tempLFP LFPtimebins LFPNormRatio L;

    %% 4-Bin spike data

    if isempty(spike_i)
        disp('No spike data was found');
    else
        for S = 1:num_spikes
            ts = datatable{spike_i(S),trial}{:}/spike_fs;
            spike_bin{trial,S} = train2bins(ts,timeframe{trial})';            
        end
    end
    
    %% 5-Bin jackpot
    if isempty(jackpot_i)
       disp('No jackpot info was found');
    else
        jackpot_bin{trial} = datatable{jackpot_i(1),trial}{:};
    end
    
end


%% Outputs
binnedData = table(jackpot_bin, force_bin, emg_bin, LFP_bin, beta_bin, lgamma_bin, hgamma_bin, vhgamma_bin, spike_bin);
        
end
