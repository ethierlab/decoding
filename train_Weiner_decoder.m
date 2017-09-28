function [decoder, varargout]=train_Weiner_decoder(binnedData, varargin)
%    [decoder] = BuildModel(binnedData, params)
%
%       binnedData              : data table to build model from
%       params                  : structure with fields: [default values in brackets]
%           numlags             : [10] decoder length in number of bins (typically 10)
%           inputs              : ['spike','beta','lgamma','hgamma','vhgamma']
%                                   signals to use for decoding in a cell array containing any combination
%                                  of {'force', 'emg', 'spike', 'LFP','beta','lgamma','hgamma','vhgamma'} strings
%           trials              : []  array of trial numbers to use for training (empty to use all trials)
%           PolynomialOrder     : [2] order of the Weiner non-linearity (0=no Polynomial) Use 2 to predict rectified signals, 3 for non-rect
%           outputs             : ['force','emg'] signals to decode in a cell array containing any combination
%                                 of {'force', 'emg', 'spike', 'LFP'} strings
%       Note on params: not all the fields have to be present in the
%       'params' structure provided in arguments. Those that are not will
%       be filled with the values from 'get_default_decoder_traing_params.m'
%
%       decoder                 : structure of decoder data (H,P,numlags)
%       varargout               : {data_fit}, a struct with the fitted predictions, as well as reformatted inputs and outputs (with intial numlags bins removed)
%
%% Argument handling

if ~istable(binnedData)
    %Load the file or structure
    binnedData = load_data_table(binnedData);
    if isempty(datatable)
        error('can''t load file');
    end
end

%update missing params with default values
params = get_default_decoder_training_params(varargin{:});

if isempty(params)
    disp('Invalid decoder training parameter(s)');
    return
end

%% Inputs & Outputs
inputs  = double(cat_data_from_bin_table(binnedData,params.inputs,params.trials));

inputs = DuplicateAndShift(inputs,params.numlags);
mx = mean(inputs); 
inputs = detrend(inputs,'constant');

outputs = double(cat_data_from_bin_table(binnedData,params.outputs,params.trials));
my = mean(outputs);
outputs = detrend(outputs,'constant');
num_out = size(outputs,2);


%% Calculate Decoder

W = train_decoder(inputs,outputs,params);

% % W_temp = inputs\outputs;
% % % add back input and output means in first row:
% % W = [ nan(1,num_out); W_temp];
% % for out = 1:num_out
% %     W(1,out) = -sum(mx'.*W_temp(:,out)) + my(out);  
% % end

[data_fit, inputs_trim, outputs_trim] = predMIMO4(inputs,W,outputs);

%% Then, add non-linearity if applicable

if params.PolynomialOrder
    P = nan(params.PolynomialOrder+1,size(data_fit,2));
    
    %%%Find a Wiener Cascade Nonlinearity
    for out=1:num_out
        %Find and apply polynomial
        [P(:,out)] = WienerNonlinearity(data_fit(:,out), outputs_trim(:,out), params.PolynomialOrder);
        data_fit(:,out) = polyval(P(:,out),data_fit(:,out));
    end
end

%% Return decoder, prediction fit

decoder = struct(...
    'W', W,...
    'P', P,...
    'params',params...
    );

if nargout > 1
    data_fit = struct(...
        'inputs'        , inputs_trim,...
        'outputs'       , outputs_trim,...
        'preds'         , data_fit,...
        'vaf'           , calc_vaf(data_fit,outputs_trim)...
        );
    varargout = {data_fit};
end

