%% Import data from text file



url = 'https://www.drugbank.ca/releases/5-1-4/downloads/all-drugbank-vocabulary';
disp(['Loading data from: ' url ]);

drugFile = unzip(url); %Get
filename = drugFile{1};

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Commonname", "Var4", "Var5", "Var6", "Var7"];
opts.SelectedVariableNames = "Commonname";
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
% [filename, pathname] = uigetfile('*.csv','Choose the drugbank vocabulary file');
% drugbankvocabulary = readtable(fullfile(pathname,filename), opts);

drugbankvocabulary = readtable(filename, opts);

%% Clear temporary variables
clear opts


%% Clean input words

doesNotContainNum = @(str) isempty(regexp(str,'[0-9]'));
%This seems a silly confusing way to remove things. But I can't be bothered
%to figure out the elegant way.
sel = table2array(rowfun(doesNotContainNum,drugbankvocabulary));
justLetterDrugs = drugbankvocabulary(sel,1);

%now convert the table to string array.
inputWordArray = table2array(justLetterDrugs);