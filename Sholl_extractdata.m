% Extract intersection information from Sholl Analysis '.csv' files

% Input: directory containing 1 folder per animal. In each folder is the
% .csv files for each microglia from the Sholl analysis
% Output: for each animal, saves the intersection matrix (i.e.
% intersections for every microglia) and the corresponding file names.
% Finds the per-animal average and the number of microglia analyzed.
% Compiles for  all animals into a matrix, and saves that matrix along with
% the animal names. 

% Get list of folders in input_dir
input_dir = uigetdir(); % Select input folder
files = dir(input_dir); % List items in input folder
% Extract only the sub-directories.
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
subFolders = files(dirFlags);
folder_list = struct2cell(subFolders);
folder_list = folder_list(1,:); % Generate list of sub-directories for input into preprocessing
num_folders = size(subFolders,1); % Number of conditions equals number sub-directories
folder_list = strtrim(folder_list); % Remove any lagging spaces

% Run through each folder to collect Sholl intersections data and then run
inters_avg_mat = []; % initialize empty matrix to gather data 
maxinters_avg_mat = []
for i=1:num_folders
    [inters_avg, maxinters_avg] = get_sholl_inters(folder_list{i},input_dir);
    inters_avg_mat(:,i) = inters_avg;
    maxinters_avg_mat(i) = maxinters_avg;
end

cd(input_dir);
results_filename = "Sholl_avgs_results.mat";
save(results_filename,"folder_list","inters_avg_mat","maxinters_avg_mat");

%{ 
% If you want a dialog box to input the animal_name manually
prompt = {'Enter animal name:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'name'};
animal_name = inputdlg(prompt,dlgtitle,dims,definput);
animal_name = animal_name{1,1};
%}

% Reads the Sholl results excel files (output from ImageJ plugin) in a
% folder and compiles the intersection data into one matrix ('inters_mat').
% Saves this to a .mat data file and outputs the average intersection data
% and the average max intersections of all files in folder.
function [inters_avg, maxinters_avg] = get_sholl_inters(animal_name,input_dir)

% Get file list
curr_dir = pwd;

cd(input_dir);
cd(animal_name);
input_names=list_directory(fullfile(input_dir,animal_name),'.csv');
input_names=strtrim(input_names);

% Initialize matrix for intersections and max intersections per microglia
n_files = size(input_names,1);
inters_mat = zeros(27,n_files); 
maxinters_mat = zeros(n_files,1);
% NOTE: Sholl anlaysis is set to do 0 to 50 um radius in 2 um increments,
% so 26 total intersection numbers. Add an extra row to inters_mat for the
% number of microglia analyzed

% Populate intersections matrix with each microglia 
for i = 1:n_files
    file_name = input_names{i};
    opts = detectImportOptions(file_name);
    opts.SelectedVariableNames = {'Inters_'};
    inters = readmatrix(file_name,opts);
    n_rows = size(inters,1);
    inters_mat(1:n_rows,i) = inters;
    maxinters_mat(i) = max(inters,[],'all'); 
end

% Calculate per animal average
inters_avg = mean(inters_mat,2);
inters_avg(27,1) = n_files;
maxinters_avg = mean(maxinters_mat,'all');
animal_file = strcat(animal_name,'.mat');
% Save file names, intersection matrix and average
save(animal_file,'input_names','inters_mat','inters_avg','maxinters_mat');

cd(curr_dir);

end


 
% Utility function that generates list of .tifs
function [ dirlist ] = list_directory( dir_path, file_ext )
% LIST specified folder and specified file name extension.
% It can be UNIX or Windows path to specify. Without specifying
% an second parameter it uses all files in that folder.
%
% Michael Tesar,
% 2016, Ceske Budejovice
% Version 1.0
%
%    INPUT:
%       dir_path - a string containing path to folder containing files
%       file_ext - a string which filter files
%
%    OUTPUT:
%       dirlist - string or char array of files meet criteria
%
% Example: [all_mp3] = list_directory ('/Users/usr/music/','.mp3');
%          [files] = list_directory ('/Users/usr/Documents/');
%
% Access data: listDir(char(index)) return string of indexed file meets
% specified criteria.
% Hint: It also can be used not for file extension but also for substring
% filtering of files.
%% Directory listing
% ==================
% Declare global manual iterator
ii = 1;
% Check if both parameters are string variable
if ischar (dir_path) && ischar(file_ext)
     listDir  = dir(dir_path);          % List directory
     listName = {listDir.name}';        % Extract only names
     allFiles = char(listName(3:end));  % Delete navigation strings . and ..
     N = length(listName(3:end));       % Get the number of total files
     % Loop for all files
     for i = 1:N
         % Check if meet file extension criteria
         if strfind(allFiles(i,:), file_ext);
             % If do then add to another variable and ++ manual iterator
             dirlist(ii,:) = {char(allFiles(i,:))};
             ii = ii + 1;
         else
             continue
         end
     end
else
    % Both parameters needs to be strings
    error('Enter strings!')
end
end
