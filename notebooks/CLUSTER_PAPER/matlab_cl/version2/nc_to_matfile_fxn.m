%requires nc toolbox, for which it's all here:
%addpath('/Users/terezajana/Sync/phd_other/classes/EOSC_510/PROJEKT/nctoolbox/')
%setup_nctoolbox

function datamat = nc_to_matfile_fxn(var,year,noday)
close all
clc

addpath('/Users/terezajana/Sync/phd_other/classes/EOSC_510/PROJEKT/nctoolbox/')
setup_nctoolbox

addpath('/Users/terezajana/Sync/phd_other/at3/CLUSTER/verze2pt0/NC_HIND/')

fold = strcat(var,int2str(year))

dirnam = strcat('/Users/terezajana/Sync/phd_other/at3/CLUSTER/verze2pt0/NC_HIND/',fold,'/',var,'_TS/')

no_var = 1;
if strcmp(var,'BIO')
    no_var = 3;
    varnames =  { 'PHY', 'PHY2',  'MYRI'};
   
end

no_var
% %variables contained in netcdf files
% 
%preallocated matrices for wind variables with zeroes
%
datamat = zeros(580,noday,no_var);
% 
for stn = 0:579
    stn
        %stitch together name of nc file
        
        
        % fwi = stn_0_fwi4m_data_sp10_threshold50.nc
        % ved = stn_0avg_ved_sp10.nc
        % wind = stn_0_wind_data_sp10.nc
        
        if strcmp(var,'FWI')
            nf = strcat(dirnam,'stn_',int2str(stn),'_fwi4m_data_sp10_threshold50.nc');
            varnames = {'freshwater_index'};
        end
        if strcmp(var,'VED')
            nf = strcat(dirnam,'stn_',int2str(stn),'avg_ved_sp10.nc');
            varnames = {'daily_ved'};
        end
        if strcmp(var,'HALO')
            nf = strcat(dirnam,'stn_',int2str(stn),'halo_depth_sp10.nc');
            varnames = {'halocline_depth'};
        end
        if strcmp(var,'WIND')
            nf = strcat(dirnam,'stn_',int2str(stn),'_wind_data_sp10.nc');
            varnames = { 'wind_energy', 'wind_stresses', 'wind_mags'};   
        end
        
        if strcmp(var,'NIT')
            nf = strcat(dirnam,'stn_',int2str(stn),'surf_nit_sp10.nc');
            varnames = {'surface_nitrogen'};
        end
            
        if strcmp(var,'BIO')
            nf = strcat(dirnam,'stn_',int2str(stn),'_sp10.nc');
            varnames =  { 'PHY', 'PHY2',  'MYRI'};   
        end
        
        
        nc1 = ncdataset(nf);
        for i = 1:no_var
            pl = nc1.data(varnames{i});
            size(pl)
            if strcmp(var,'BIO')
                datamat(stn+1,:,i) = pl(1:noday);
            else
                datamat(stn+1,:) = pl;
            end
        end

end

if strcmp(var,'BIO')
    datamat_linear = zeros(580,noday*3);
    datamat_linear(:,1:noday) = datamat(:,:,1);
    datamat_linear(:,noday+1:noday*2) = datamat(:,:,2);
    datamat_linear(:,noday*2+1:end) = datamat(:,:,3);
end
return 


% %for biology only:
% datamat_linear = zeros(580,365*3);
% datamat_linear(:,1:365) = datamat(:,:,1);
% datamat_linear(:,366:365*2) = datamat(:,:,2);
% datamat_linear(:,365*2+1:end) = datamat(:,:,3);



