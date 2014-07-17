function [xlsarm, xlshn] = erymodel_fit()
% ERYMODEL_FIT
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

%% Step 0: Load fitting files & set starting parameters
load('sbse_coeffs.mat')
load('model_params.mat')
param0 = [0.0076, 0.8421, -0.0017, 5.0798, 0.0321];
paramarmin = param0;
paramhnin = param0;
lb = [0, 0.5, -0.01, 0, -0.1];
ub = [0.08, 1, 0.01, 50, 0.4];
%options = optimset('Algorithm','levenberg-marquardt');
options = optimset('Algorithm','trust-region-reflective');
fracstdev = mua_param(1,:)/5000+0.005;

%% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, subfolders] = select_folder;

% For each subfolder (day), (subfolder format: Day# (no padding))

jday = 1; %Use to increment excel column
for iday=1:35
    iday
    %% Step 2: Obtain signal spectra
    %curdir = strcat(pathname,subfolders{iday},'\');
    curdir = strcat(pathname,'Day',num2str(iday),'\');

    if exist(curdir) == 7 % if day's subfolder exists
        [rarm, rhn] = retrieve_spec(pathname, curdir, iday); % Step 2.1-2.5: Get R_arm & R_meas
        
            %% Step 3: Correct rmeas for SBSE
        rarmsbse = corr_sbse(rarm, sbse_coeffs);
        stdevarm = fracstdev.*rarmsbse;
        [paramarmout,~,residual,~,~,~,jacobian] = lsqcurvefit(@(paramarmin,lambda_param) calc_rd(paramarmin,lambda_param,mua_param,musp),paramarmin,lambda_param,rarmsbse,lb,ub,options);
        eicorrarm = calc_ei(lambda_param, rarmsbse);

        xlsarm(jday,:) = [iday, paramarmout, eicorrarm];
        paramarmin = paramarmout;

        rhnsbse = corr_sbse(rhn, sbse_coeffs);
        stdevhn = fracstdev.*rhnsbse;
        [paramhnout,~,residual,~,~,~,jacobian] = lsqcurvefit(@(paramhnin,lambda_param) calc_rd(paramhnin,lambda_param,mua_param,musp),paramhnin,lambda_param,rhnsbse,lb,ub,options);
        eicorrhn = calc_ei(lambda_param, rhnsbse);
     
        xlshn(jday,:) = [iday, paramhnout, eicorrhn];
        paramhnin = paramhnout;

%     rhnslcf = calc_rd(paramhn,lambda_param,mua_param,musp);
%     figure
%     plot(lambda_param,rhnsbse,lambda_param,rhnslcf)
%     rarmslcf = calc_rd(paramarm,lambda_param,mua_param,musp);
%     figure
%     plot(lambda_param,rarmsbse,lambda_param,rarmslcf)
        
        jday = jday + 1;
    else
        %???
    end

%     allrarm(:,iday) = rarm;
%     allrmeas(:,iday) = rmeas;


    

    
end

print_to_excel(pathname, xlsarm, xlshn)

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

function [rarm, rhn] = retrieve_spec(pathname, curdir, day)

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

rhn = (allspec(:,4)-allspec(:,1))./(allspec(:,2)-allspec(:,1));

end

function rsbse = corr_sbse(rmeas, sbse_coeffs)

rsbse = ((sbse_coeffs(:,1).*rmeas + sbse_coeffs(:,2))./(rmeas + sbse_coeffs(:,3)))';

end

function rslcf = calc_rd(beta,lambda,extco,musp)
%RDSCAT Calculate diffuse reflectance with added scatter losses

% Step 1: Calculate diffuse reflectance using model (rcalc)

% beta(1) = total hemoglobin
% beta(2) = oxygen saturation
% beta(3) = melanin
% beta(4) = background
% beta(5) = background shift

mua = beta(1)*beta(2)*extco(1,:) + beta(1)*(1-beta(2))*extco(2,:) + beta(3)*extco(3,:) + beta(4)*extco(4,:) + beta(5);

mutp = mua + musp;
mueff = sqrt(3.*mua.*mutp);
D = 1./(3.*mutp);

%A = 2.348; % for an nrel = 1.33 (for water)
A = 2.745; % for an nrel = 1.4 (for tissue)

rcalc = musp./((mutp + mueff).*(1 + 2.*A.*D.*mueff));

% Step 2: Convolve spectrum with Gaussian
fwhm = 10;
mu = (lambda(end)-lambda(1))/2 + lambda(1);
goflambda = exp(-(lambda-mu).^2/(2*(fwhm/2.3548)^2));
result = conv(goflambda,rcalc);

cones = ones(1,length(lambda));
cal = conv(goflambda,cones);
norm = result./cal;

rconv = norm(302:918);

% Step 3: Calculate and add in scatter losses

for i=1:length(lambda)
    if mua(i) <= 0.05*musp(i)
        slcf(i) = 0.044*musp(i).^-0.6.*log(mua(i)) + 1.04;
         if slcf(i) >= 1, slcf(i) = 1; end
    else
        slcf(i) = 1;
    end
end

rslcf = slcf.*rconv;

end

%function [eicorr, del_eicorr] = calc_ei(lambda, rsbse, std)
function eicorr = calc_ei(lambda, rsbse)

% Calculate log inverse reflectance
lir = log10(1./rsbse);
%del_lir = 0.43*(std./rsbse);

% Determine location and values of EI input
lambdaei = [510 543 563 576 610];
for j=1:length(lambdaei)
    [~, eminloc(j)] = min(abs(lambda-lambdaei(j)));
    eival(j) = mean(lir(eminloc(j)-2:eminloc(j)+2));
    %del_eival(j) = mean(del_lir(eminloc(j)-2:eminloc(j)+2));
end

% Calculate EI
ei = 100*(eival(3) + 1.5*(eival(2) + eival(4)) - 2*(eival(1) + eival(5)));
%del_ei = 100*(del_eival(3) + 1.5*(del_eival(2) + del_eival(4)) + 2*(del_eival(1) + del_eival(5)));

% Determine location and values of MI input
lambdami = [650 700];
for k=1:length(lambdami)
    [~, mminloc(k)] = min(abs(lambda-lambdami(k)));
    mival(k) = mean(lir(mminloc(k)-15:mminloc(k)+15));
    %del_mival(k) = mean(del_lir(mminloc(k)-15:mminloc(k)+15));
    
end

% Calculate MI
mi = 100*(mival(1)-mival(2) + 0.015);
%del_mi = 100*(del_mival(1) + del_mival(2));

% Apply MI to get corrected EI
eicorr = ei*(1+0.04*mi);
%del_eicorr = eicorr*(del_ei/ei + 0.04*(del_ei/ei + del_mi/mi));

end

function [] = print_to_excel(pathname, xlsarm, xlshn)

pathsplit = strsplit(pathname,'\');
xlsfilename = pathsplit{end-1};
xlsname = strcat(pathname, pathsplit{end-1}, '(2).xlsx');
%xlsnamestr = xlsname{1};

headings = {'Day', 'tHb', 'SO2', 'mel', 'bkg', 'bkg_shift', 'EIcorr'};

xlswrite(xlsname, headings, 'Arm', 'A1')
xlswrite(xlsname, xlsarm, 'Arm', 'A2')

xlswrite(xlsname, headings, 'HN', 'A1')
xlswrite(xlsname, xlshn, 'HN', 'A2')


end