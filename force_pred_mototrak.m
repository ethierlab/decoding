function [vaf, act_force, pred_force] = force_pred_mototrak(binned_data,varargin)
% Created: Tuesday, June 12, 2019.
% Last edited: Wednesday, June 27, 2019. By Vajra.
%
%    This function trains neural firing rates on lever force and tests a decoder using k-fold
%    cross-validation.
%
%    USAGE: force_pred_mototrak(binned_data,'numfolds',int,'numlags',int,'zero_wind',int)
%
%    binned_data : processed from data2 into bins
%    vaf         : variance accounted for
%    act_force   : actual force
%    pred_force  : predicted force
%
%    numfolds    : number of k-folds to be used in cross-validation
%    numlags     : number of datasets which are stacked and shifted across
%    zero_wind   : bin size of the zero window for eliminating zero force
%
%
%%%% EthierLab 12/06/12 -- CE & VK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This sample file was used for reference to some comments below: jados-19-03-25-pm.mat

%% Argument handling

% Default parameters
params = struct(...
    'numfolds'         ,10, ... %sets 10 as default if params is not passed
    'numlags'          , 5,  ...
    'zero_wind'        ,10);

params = parse_input_params(params,varargin);


vaf        = nan(params.numfolds,1); %initializes vaf table
pred_force = []; %initializes pred_force array

num_trials = size(binned_data,1); %size of binned_data; 60 rows
all_trials = 1:num_trials;


%% pre-process data
% (remove data with zero force and duplicate and shift data)

processed_spikes = cell(num_trials,1); %initiates a cell array with length of num_trials (60)
processed_force  = cell(num_trials,1); %initiates a cell array with length of num_trials (60)

for t = 1:num_trials
    spikes   = cell2mat(binned_data.spike_bin(t,:)); %converts cell array into an ordinary array
    spikes = DuplicateAndShift(spikes,params.numlags); %duplicates the spikes and stacks them together, shifting by one bin each row up to 5 bins
    spikes = spikes(params.numlags:end,:); % subtracts the first 5 columns (corresponding to our blanks) from spike data i.e. the shift
    force    = binned_data.force_bin{t}(params.numlags:end); %subtracts the first 5 columns (corresponding to our shift) from force data
    
    % remove data where moving average of force is zero.
    valid_idx = movmean(force,params.zero_wind)>2; %sets up a moving mean. see >>help movmean.
    processed_force{t}     = force(valid_idx); %select forces only where valid_idx >2
    processed_spikes{t}    = spikes(valid_idx,:);%select spikes only where valid_idx >2
end

%% train and test decoder using multi-fold xval

%num_train_trials = round((1-1/params.numfolds)*num_trials); % 54; 10 folds should fit into 60

%num_test_trials = num_trials - num_train_trials;

num_test_trials  = floor(num_trials/params.numfolds);% 6; size of each training set

for f = 1:params.numfolds
    
    test_trials = (1:num_test_trials)+(f-1)*num_test_trials; %from 1,2,3,4,5,6 to 7,8,9,10,11,12,13...
    train_trials = all_trials(~ismember(all_trials,test_trials)); %this will shift forward through all_trials by one set (6) each loop according to test_trial
    
    inputs  = vertcat(processed_spikes{train_trials}); %vertically concatenates processed_spikes (where train_trials are located); 2778x105 double
    outputs = vertcat(processed_force {train_trials}); %vertically concatenates processed_force (where train_trials are located); 2778x1 double
    
    % train decoder
    W = filMIMO4(inputs,outputs,1,1,1); %inputs = columnwise input X, outputs = columnwise input Y. But why is params.numlags = 1??
    
    % test decoder
    test_spikes = vertcat(processed_spikes{test_trials}); %processed_spikes{1,2,3,4,5,6} then processed_spikes{7,8,9,10,11,12}...
    test_force  = vertcat(processed_force {test_trials});
    
    pred_force_temp = predMIMO4(test_spikes,W,test_force);
    
    % calc vaf
    vaf(f) = calc_vaf(pred_force_temp,test_force); %this will calculate the model's performance
    
    pred_force = [pred_force; pred_force_temp];
    
    %fprintf('fold %d :  vaf = %.2f\n',f, vaf(f));
    
    
end

assignin('base','vaf_per_fold', vaf)
assignin('base','mean_vaf', nanmean(vaf))
assignin('base','stdev_vaf', nanstd(vaf))
assignin('base','act_force',vertcat(processed_force{:}))
assignin('base','pred_force',pred_force)


% plot(act_force(1:200))
% hold on
% plot (pred_force(1:200))
% hold off

