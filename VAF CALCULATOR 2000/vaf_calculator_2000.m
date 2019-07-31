function varargout = vaf_calculator_2000(varargin)
% VAF_CALCULATOR_2000 MATLAB code for vaf_calculator_2000.fig
%      VAF_CALCULATOR_2000, by itself, creates a new VAF_CALCULATOR_2000 or raises the existing
%      singleton*.
%
%      H = VAF_CALCULATOR_2000 returns the handle to a new VAF_CALCULATOR_2000 or the handle to
%      the existing singleton*.
%
%      VAF_CALCULATOR_2000('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VAF_CALCULATOR_2000.M with the given input arguments.
%
%      VAF_CALCULATOR_2000('Property','Value',...) creates a new VAF_CALCULATOR_2000 or raises the
%      existing singleton*.  Starting from the left, property Value pairs are
%      applied to the GUI before vaf_calculator_2000_OpeningFcn gets called.  An
%      unrecognized property name or invalid Value makes property application
%      stop.  All inputs are passed to vaf_calculator_2000_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vaf_calculator_2000

% Last Modified by GUIDE v2.5 25-Jul-2019 19:09:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @vaf_calculator_2000_OpeningFcn, ...
    'gui_OutputFcn',  @vaf_calculator_2000_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before vaf_calculator_2000 is made visible.
function vaf_calculator_2000_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vaf_calculator_2000 (see VARARGIN)

% Choose default command line output for vaf_calculator_2000
handles.output = hObject;

handles.loaded_data = [];
handles.binnedData = [];

set(handles.file_path, 'String', 'File path...');

set(handles.mu,'Value',1); 
set(handles.su, 'Value',1);

set(handles.numfolds_custom, 'String', 10);
set(handles.numlags_custom, 'String', 5);
set(handles.zero_custom, 'string', 10);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vaf_calculator_2000 wait for user response (see UIRESUME)
% uiwait(handles.vaf_table);


% --- Outputs from this function are returned to the command line.
function varargout = vaf_calculator_2000_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function file_path_Callback(hObject, eventdata, handles)
% hObject    handle to file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_path as text
%        str2double(get(hObject,'String')) returns contents of file_path as a double


%% Upload file
% --- Executes on button press in upload_file_button.
function upload_file_button_Callback(hObject, eventdata, handles)
% hObject    handle to upload_file_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if ~isempty(handles.file_path)
%     handles.loaded_data = load(get(handles.file_path,'string'))
% else
%     return
% end

[filename, pathname] = uigetfile(fullfile('*.mat'),'Upload matlab file that contains data2');

filelocation = fullfile(pathname,filename);

handles.loaded_data = load(filelocation);

% Fills file location
if ~isfield(handles.loaded_data,'data2')
    msgbox({'This file did not contain data2.'; 'Upload was unsuccessful.'},'Error');
    return
else
    msgbox('Upload successful!', 'Success')
end

set(handles.file_path, 'string', filelocation);



% Update handles structure
guidata(hObject, handles);



%% Process data

% --- Executes on button press in mu.
function mu_Callback(hObject, eventdata, handles)
% hObject    handle to mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mu

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in su.
function su_Callback(hObject, eventdata, handles)
% hObject    handle to su (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of su

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in process_data_button.
function process_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to process_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check and convert analog data to force if necessary
message = [];

if isempty(handles.loaded_data)
    message = msgbox('There was nothing to process!','Error');
    return
end

% Checkboxes 
mu_value = get(handles.mu, 'Value');
su_value = get(handles.su, 'Value');

if mu_value && su_value
    type = 'tout';
    
elseif mu_value && su_value == 0
    type = 'mu';
    
elseif su_value && mu_value == 0
    type = 'su';
    
else
    message = msgbox('Please select one or both options!', 'Error');
    return
end

data2 = handles.loaded_data.data2;

% Remove empty columns from data2
for var = 1:size(data2,2)
    is_empty(var,1) = isempty(data2{4,var}{:});
end

if any(is_empty ~= 0)
    data2 = removevars(data2,data2.Properties.VariableNames{(is_empty)});
end

% If data2 has not been converted...
if max(data2{4,1}{:}) < 5
    message = msgbox('Converting analogue forces to grams...','Loading');
    data2 = force_analog2grams(data2,handles.loaded_data.fs);
end

% Create binnedData
bin_size=.05; %s
params.binsize=bin_size;
handles.binnedData = extraire_bins_pour_force_et_clusters(data2,handles.loaded_data.fs,type,'initiation',params);

delete(message);

if ~isempty(handles.binnedData)
    msgbox('Processing complete!', 'Success')
end

% Update handles structure
guidata(hObject, handles);


function numfolds_custom_Callback(hObject, eventdata, handles)
% hObject    handle to numfolds_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numfolds_custom as text
%        str2double(get(hObject,'String')) returns contents of numfolds_custom as a double



% --- Executes during object creation, after setting all properties.
function numfolds_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numfolds_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function numlags_custom_Callback(hObject, eventdata, handles)
% hObject    handle to numlags_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numlags_custom as text
%        str2double(get(hObject,'String')) returns contents of numlags_custom as a double


% --- Executes during object creation, after setting all properties.
function numlags_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numlags_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zero_custom_Callback(hObject, eventdata, handles)
% hObject    handle to zero_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zero_custom as text
%        str2double(get(hObject,'String')) returns contents of zero_custom as a double


% --- Executes during object creation, after setting all properties.
function zero_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Calculate VAF
% --- Executes on button press in calculate_vaf_button.
function calculate_vaf_button_Callback(hObject, eventdata, handles)
% hObject    handle to calculate_vaf_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

message = [];

if isempty(handles.binnedData)
    message = msgbox('There was nothing to calculate!','Error');
    return
end

% Parameters to customize
numfolds = str2double(get(handles.numfolds_custom,'String'));
numlags = str2double(get(handles.numlags_custom,'String'));
zero_wind = str2double(get(handles.zero_custom,'String'));

% Predicted vs actual force
pred_force = [];
num_trials = size(handles.binnedData,1);
all_trials = 1:num_trials;

processed_spikes = cell(num_trials,1);
processed_force  = cell(num_trials,1);
vaf = zeros(numfolds,1);

for t = 1:num_trials
    spikes   = cell2mat(handles.binnedData.spike_bin(t,:));
    spikes = DuplicateAndShift(spikes,numlags);
    spikes = spikes(numlags:end,:);
    force    = handles.binnedData.force_bin{t}(numlags:end);
    
    % remove data where moving average of force is zero.
    valid_idx = movmean(force,zero_wind)>2;
    processed_force{t}     = force(valid_idx);
    processed_spikes{t}    = spikes(valid_idx,:);
end

num_test_trials  = (num_trials - mod(num_trials,numfolds))/numfolds;

for f = 1:numfolds
    
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
    vaf(f,1) = f;
    vaf(f,2) = calc_vaf(pred_force_temp,test_force);
    
    pred_force = [pred_force; pred_force_temp];
    
end

act_force = vertcat(processed_force{:});

% Fill table with fold # vs VAF
set(handles.mean_vaf_table, 'data', vaf);

% Mean VAF
set(handles.mean_vaf, 'string', mean(vaf(:,2)));

% SD VAF
set(handles.sd_vaf, 'string', std(vaf(:,2)));

axes(handles.force_figure);
% Plot actual force
plot(act_force(1:200),'-k');
xlabel('Bins')
ylabel('Force (g)')

hold on
% Plot predicted force
plot(pred_force(1:200),'-r');
legend('Actual force','Predicted force')

hold off

  message = msgbox('Calculations complete!', 'Success');
  
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function force_figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate force_figure


function mean_vaf_Callback(hObject, eventdata, handles)
% hObject    handle to mean_vaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_vaf as text
%        str2double(get(hObject,'String')) returns contents of mean_vaf as a double


% --- Executes during object creation, after setting all properties.
function mean_vaf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_vaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sd_vaf_Callback(hObject, eventdata, handles)
% hObject    handle to sd_vaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sd_vaf as text
%        str2double(get(hObject,'String')) returns contents of sd_vaf as a double


% --- Executes during object creation, after setting all properties.
function sd_vaf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sd_vaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in optimum_settings.
function optimum_settings_Callback(hObject, eventdata, handles)
% hObject    handle to optimum_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.binnedData)
    message = msgbox('There was nothing to calculate!','Error');
    return
end

numfolds = 10;
numlags = 10;
zero_wind = 20;
extra = 9;

loading = waitbar(0.1, 'Calculating number of lags...');

numlags_vaf = nan(numfolds,numlags + extra);

% Calculate max VAF for numlags
for i = numlags - extra : numlags + extra
    numlags_vaf(:,i) = force_pred_mototrak(handles.binnedData,'numlags',i);
end

numlags_mean = nanmean(numlags_vaf);

[~,max_numlags] = find(numlags_mean == max(numlags_mean));
set(handles.numlags_custom,'string',max_numlags)

waitbar(0.5, loading, 'Calculating size of zero window...');

% Calculate max VAF for zero_bin
 zero_wind_vaf = nan(numfolds,zero_wind + extra);

for i = zero_wind - extra:zero_wind + extra
    zero_wind_vaf(:,i) = force_pred_mototrak(handles.binnedData,'zero_wind',i);
end

zero_wind_mean = mean(zero_wind_vaf);
[~,max_zero_wind] = find(zero_wind_mean == max(zero_wind_mean));
set(handles.zero_custom,'string',max_zero_wind)

waitbar(1, loading, 'Calculations complete!');


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
