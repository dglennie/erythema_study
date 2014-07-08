function [] = erymodel_fit()
% ERYMODEL_FIT
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

%% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, subfolders] = select_folder;

% For each subfolder (day), (subfolder format: Day# (no padding))

for i=1:length(subfolders)

    %% Step 2: Obtain signal spectra
    curdir = strcat(pathname,subfolders{1},'\');
    curcont = dir(curdir);
    cellstruc = struct2cell(curcont);
    cellstruc1 = cellstruc(1,:);
    curcont = cellstruc1(3:end);

    % Step 2.1: Pull 3 bkg files & average
    % Check that all 3 files are present (error and skip if not)
    exist_b(1) = cell2mat(strfind(curcont, 'b1.Master.Scope'));
    exist_b(2) = cell2mat(strfind(curcont, 'b2.Master.Scope'));
    exist_b(3) = cell2mat(strfind(curcont, 'b3.Master.Scope'));
    
    if sum(exist_b) == 3
        % a measurements exist
    else
        % don't exist, skip
    end
    % Pull 3 spec & average

    % Step 2.2: Pull 3 cal files & average

    % Step 2.3: Pull 3 arm files & average

    % Step 2.4: Pull 3 meas files & average

    %% Step 3: Extract parameters for arm & meas spectra
    
end

function [pathname, subfolders] = select_folder
%SELECT_FOLDER Get a list of files from selected folder

foldername = uigetdir;
pathname = strcat(foldername,'\');
contents = dir(foldername);

strucell = struct2cell(contents);
strucell1 = strucell(1,:);

subfolders = strucell1(3:end);

end