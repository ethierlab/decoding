function params = get_default_binning_params(varargin)

if nargin
    params = varargin{1};
    if nargin > 1
        warning('Wrong number of arguments');
        evalin('base','help convertBDF2binned');
        return;
    end
else
    params = [];
end

%% Default Parameters (all units are in seconds):
params_defaults = struct(...
    'binsize'       , 0.05,...
    'pre_capture'   , 0.5,...
    'EMG_hp'        , 50,...
    'EMG_lp'        , 10,...
    'NormData'      , false,...
    'ArtRemEnable'  , false,...
    'NumChan'       , 10,...
    'TimeWind'      , 0.0005);


%% Update missing values with defaults
all_param_names = fieldnames(params_defaults);
for i=1:numel(all_param_names)
    if ~isfield(params,all_param_names(i))
        params.(all_param_names{i}) = params_defaults.(all_param_names{i});
    end
end

%% Parameter Validation

% if mod(1,binsize)
%     disp('Please choose a binsize that is a factor of 1');
%     disp('data conversion aborted');
%     params = [];    
%     return
% end

