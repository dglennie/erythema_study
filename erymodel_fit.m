function [allrarm, allrmeas] = erymodel_fit()
% ERYMODEL_FIT
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

%% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, subfolders] = select_folder;

% For each subfolder (day), (subfolder format: Day# (no padding))

%jday = 1; Use to increment excel column
iday = 1; %for day=1:length(subfolders) %!!!! What if a day/subfolder is missing?

    %% Step 2: Obtain signal spectra
    curdir = strcat(pathname,subfolders{iday},'\');

    if exist(curdir) == 7 % if day's subfolder exists
        % Step 2.1-2.5: Get R_arm & R_meas
        [rarm, rmeas] = retrieve_spec(pathname, curdir, iday);
    else
        %???
    end

    allrarm(:,iday) = rarm;
    allrmeas(:,iday) = rmeas;

    %% Step 3: Extract parameters for arm & meas spectra
    
%end

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

function [rarm, rmeas] = retrieve_spec(pathname, curdir, day)

% Step 2.1-2.4: Pull rate spectra and average

pre = {'b', 'c', 'a', 'm'};

for jspec = 1:4
    kspec = 1;
    for ispec = 1:3
        curfilename = strcat(pre{jspec}, num2str(ispec), '.Master.Scope');
        if exist(strcat(curdir,curfilename)) == 2 % name exists
            % pull spec, convert to rate
            curfile = strcat(pathname, 'Day', num2str(day), '\', curfilename);
            data = dlmread(curfile,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
            inttime = dlmread(curfile,' ', [6,3,6,3]); % reads in the integration time, space delimited
            spec1 = (data(:,2)/(inttime/1000));
            spec2(:,kspec) = spec1(453:1069);
            kspec = kspec + 1;
        else
            % skip
        end
    end
    allspec(:, jspec) = mean(spec2,2); % average 3 measurements
end

% Step 2.5: Calculate arm and measurement reflectance

rarm = (allspec(:,3)-allspec(:,1))./(allspec(:,2)-allspec(:,1));

rmeas = (allspec(:,4)-allspec(:,1))./(allspec(:,2)-allspec(:,1));

end