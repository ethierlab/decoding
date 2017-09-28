function data_array = cat_data_from_bin_table(bd, var_names, varargin)
% usage: data_array = get_data_from_bin_table(tb, var_names, [trials])
%   this function extracts and concatenates data from a bin table (see convertMDF2binned.mat)
%   
%   data_array      :   returns an array of doubles containing all the data specified in var_names,
%                       concatenated from all the specified trials.
%   
%   bd              :   binned data table
%   var_names       :   string array with variable names to be extracted ( e.g. {'force', 'emg', 'LFP'} )
%                       (uses only first three letters, case insensitive)
%   trials          :   optional argument specifying specific trial numbers (unspecified: all trials)
%

data_array = [];
num_data_type = size(var_names,2);
num_trials = size(bd,1);

if nargin == 2
    trials =1:num_trials;
elseif nargin == 3
    trials = varargin{1};
    if isempty(trials)
        trials = 1:num_trials;
    end
else
    error('In:cat_data_from_bin_table, wrong number of arguments');
end 
    
for i=1:num_data_type
    data_type_i    = strncmpi(bd.Properties.VariableNames,var_names{i},3);
    num_sub_data_type = size(bd{1,data_type_i},2);
    for j=1:num_sub_data_type
        data_array = [data_array vertcat(bd{:,data_type_i}{trials,j})];
    end
end
