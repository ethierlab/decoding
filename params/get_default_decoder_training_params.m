function params = get_default_decoder_training_params(varargin)

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
        'numlags'         ,10,...
        'inputs'          ,{{'spike','beta','lgamma','hgamma','vhgamma'}},... %%{'force', 'emg', 'spike', 'LFP','beta','lgamma','hgamma','vhgamma'
        'trials'          ,[],...
        'PolynomialOrder' ,2,...
        'outputs'         ,{{'force','emg'}},...
        'lambda'          , 0 ...  %L2 regul
        );


%% Update missing values with defaults
all_param_names = fieldnames(params_defaults);
for i=1:numel(all_param_names)
    if ~isfield(params,all_param_names(i))
        params.(all_param_names{i}) = params_defaults.(all_param_names{i});
    end
end