function [vafs, act_force, pred_force] = force_pred_mototrak(binned_data)

% todo: pass as arguments
num_folds = 10;
numlags   = 5;
zero_wind = 10;
vafs      = nan(num_folds,1);
pred_force = [];

num_trials = size(binned_data,1);
all_trials = 1:num_trials;

%% pre-process data
% (remove data with zero force and duplicate and shift data)

processed_spikes = cell(num_trials,1);
processed_force  = cell(num_trials,1);

for t = 1:num_trials
    spikes   = cell2mat(binned_data.spike_bin(t,:));
    spikes = DuplicateAndShift(spikes,numlags);
    spikes = spikes(numlags:end,:);
    force    = binned_data.force_bin{t}(numlags:end);
    
    % remove data where moving average of force is zero.
    valid_idx = movmean(force,zero_wind)>2;
    processed_force{t}     = force(valid_idx);
    processed_spikes{t}    = spikes(valid_idx,:);
end

%% train and test decoder using multi-fold xval

num_train_trials = round((1-1/num_folds)*num_trials);
num_test_trials  = num_trials - num_train_trials;

for f = 1:num_folds
    
    test_trials = (1:num_test_trials)+(f-1)*num_test_trials;
    train_trials = all_trials(~ismember(all_trials,test_trials));
    
    inputs  = vertcat(processed_spikes{train_trials});
    outputs = vertcat(processed_force {train_trials});
     
    % train decoder
    W = filMIMO4(inputs,outputs,1,1,1);
    
    % test decoder
    test_spikes = vertcat(processed_spikes{test_trials});
    test_force  = vertcat(processed_force {test_trials});
    
    pred_force_temp = predMIMO4(test_spikes,W,test_force);
    
    % calc vaf
    vafs(f) = calc_vaf(pred_force_temp,test_force);
    
    pred_force = [pred_force; pred_force_temp];
    
    fprintf('fold %d :  vaf = %.2f\n',f, vafs(f));
    
end

act_force = vertcat(processed_force{:});
fprintf('\nmean vaf: %.2f\n',mean(vafs));
    
   

