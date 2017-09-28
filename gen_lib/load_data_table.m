function [datatable,fs] = load_data_table(filename)
% loads a mat file
% usage:[datatable,fs] = load_data_table(filename)
%       filename   : path and filename of the mat file to be loaded
%                    or structure of table in workspace
%       datatable  : loaded trial table
%       fs         : loaded table of sampling frequencies
%
       
%default values:
datatable = struct([]);
if (nargin==0)
    disp('Please provide a string containing the file name');
    disp('of a data table in workspace or name and full path of a .mat file');
    return
end

WS = 'base'; %use structure in workspace

assignin(WS,'structname',filename);

if (evalin(WS,'exist(structname,''var'')'))
    fprintf('Using structure ''%s'' which was already in workspace\n',filename);
    datatable = evalin(WS,filename);
    fs = evalin(WS,'fs');
    evalin(WS,'clear structname');
    return
end
evalin(WS,'clear structname');

% if not, check that it's a file
if ~(exist(filename,'file')) 
    fprintf('%s ain''t no file or structure I ever heard of\n',filename);
    return
end

% assign content of structure from file to output
datastruct  = load(filename);
field_names = fieldnames(datastruct);
for i = 1:size(field_names,1)
    if strcmpi(field_names{i},'data')
        datatable = getfield(datatable, field_names{i});
    elseif strcmpi(field_names{i},'fs')
        fs = getfield(datatable, field_names{i});
    else
        warning('in:load_data_table; Unhandeled data field:%s',field_names{i});
    end
end




