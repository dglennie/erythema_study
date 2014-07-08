function [] = erymodel_fit()
% ERYMODEL_FIT
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

%% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, subfolders] = select_folder;

% Step 1.2: Determine number of measurements/days/subfolders
% Subfolder format: Day# (no padding)

% For each subfolder (day),

%% Step 2: Obtain signal spectra

% Step 2.1: Pull 3 bkg files & average

% Step 2.2: Pull 3 cal files & average

% Step 2.3: Pull 3 arm files & average

% Step 2.4: Pull 3 meas files & average

%% Step 3: Extract parameters for arm & meas spectra

function [pathname, subfolders] = select_folder
%SELECT_FOLDER Get a list of files from selected folder

foldername = uigetdir;
pathname = strcat(foldername,'\');
contents = dir(foldername);

strucell = struct2cell(contents);
strucell1 = strucell(1,:);

subfolders = strucell1(3:end);

end