%% ============================================================
%  Paths & colormaps
%% ============================================================
cd("/home/donglaiyang/Documents/Georgia-Tech/Research/ASE-perturbed")
CESM_path = '/data-archive/ISMIP6-climate-forcing/CESM2-WACCM_ssp585/';
CCSM_path = '/data-archive/ISMIP6-climate-forcing/CCSM4-RCP85/';

colormap_path = '/home/donglaiyang/Documents/Georgia-Tech/Research/common-data-set/ScientificColourMaps8/';
load([colormap_path 'nuuk/nuuk.mat'])
load([colormap_path 'vik/vik.mat'])

secinyear = 365*24*3600;
rho_ice   = 917;
dt_target = 10;     % years -> subsampling interval
refYear   = 2015;   % new anomaly baseline

pp = projcrs(3031); % EPSG:3031, Antarctic Polar Stereographic

%% ============================================================
%  One struct per model describing its files
%% ============================================================
models(1).name     = 'CESM2-WACCM';
models(1).f_to2100 = [CESM_path 'CESM2-WACCM_16km_anomaly_ssp585_1995-2100.nc'];
models(1).f_to2300 = [CESM_path 'CESM2-WACCM_16km_anomaly_ssp585_2101-2299.nc'];
models(1).outfile  = 'data/projection-climate-forcing/CESM2-WACCM_16km_anomaly2015_ssp585_1995-2300.mat';

models(2).name     = 'CCSM4';
models(2).f_to2100 = [CCSM_path 'CCSM4_16km_anomaly_1995-2100.nc'];
models(2).f_to2300 = [CCSM_path 'CCSM4_16km_anomaly_2101-2300.nc'];
models(2).outfile  = 'data/projection-climate-forcing/CCSM4_16km_anomaly2015_RCP85_1995-2300.mat';

%% ============================================================
%  Main loop: process each model
%% ============================================================
for m = 1:numel(models)
    modelName = models(m).name;
    fprintf('Processing %s ...\n', modelName);

    % ---- Period 1 (...to 2100) ----
    f1   = models(m).f_to2100;
    lat1 = ncread(f1,'lat');
    lon1 = ncread(f1,'lon');
    smb1 = ncread(f1,'smb_anomaly');
    ts1  = ncread(f1,'ts_anomaly');
    yr1  = local_time2year(f1);

    % ---- Period 2 (...to 2300) ----
    f2   = models(m).f_to2300;
    lat2 = ncread(f2,'lat');
    lon2 = ncread(f2,'lon');
    smb2 = ncread(f2,'smb_anomaly');
    ts2  = ncread(f2,'ts_anomaly');
    yr2  = local_time2year(f2);

    % --- sanity check: same spatial grid across periods ---
    if ~isequal(lat1,lat2) || ~isequal(lon1,lon2)
        error('Grid mismatch between the two periods for %s', modelName);
    end
    lat = lat1; lon = lon1;

    % ---- Concatenate along time ----
    smb_all = cat(3, smb1, smb2);
    ts_all  = cat(3, ts1,  ts2);
    yr_all  = [yr1(:); yr2(:)];

    % guard against a duplicated boundary year
    [yr_all, uidx] = unique(yr_all, 'stable');
    smb_all = smb_all(:,:,uidx);
    ts_all  = ts_all(:,:,uidx);

    % ---- Re-baseline anomaly to 2015 ----
    idx2015 = find(yr_all == refYear, 1);
    if isempty(idx2015)
        error('Year %d not found in combined series for %s', refYear, modelName);
    end
    smb_base = smb_all - smb_all(:,:,idx2015);
    ts_base  = ts_all  - ts_all(:,:,idx2015);

    % ---- Unit conversion for SMB: kg m^-2 s^-1  ->  m ice yr^-1 ----
    smb_base = smb_base * secinyear / rho_ice;

    % ---- Subsample to every 10 years, anchored on 2015 ----
    keep    = mod(yr_all - refYear, dt_target) == 0;
    yr_dec  = yr_all(keep);
    smb_dec = smb_base(:,:,keep);
    ts_dec  = ts_base(:,:,keep);

    % ---- Project coordinates once ----
    [xg, yg] = projfwd(pp, lat(:), lon(:));

    % only keep year 2015 and forward
    idx2015 = find(yr_dec == refYear,1);
    yr_dec  = yr_dec(idx2015:end);
    smb_dec = smb_dec(:,:,idx2015:end);
    ts_dec  = ts_dec(:,:,idx2015:end);

    % ---- Package & save ----
    S.model        = modelName;
    S.lat          = lat;
    S.lon          = lon;
    S.x            = xg;
    S.y            = yg;
    S.year         = yr_dec;
    S.smb_anom2015 = smb_dec;   % m ice eq. / year, relative to 2015
    S.ts_anom2015  = ts_dec;    % K, relative to 2015

    save(models(m).outfile, '-struct', 'S');
    models(m).S = S;   % keep in memory for plotting below
end

%% ============================================================
%  Plot 1: SMB anomaly, 2015 vs 2300, rows = model
%% ============================================================
figure;
for m = 1:numel(models)
    S = models(m).S;
    idxStart = find(S.year == refYear, 1);
    idxEnd   = numel(S.year);

    subplot(2,2,(m-1)*2+1)
    slice = S.smb_anom2015(:,:,idxStart);
    scatter(S.x, S.y, 8, slice(:), 'filled');
    colormap(vik); colorbar; axis equal tight
    title(sprintf('%s SMB anomaly, %d', S.model, S.year(idxStart)))

    subplot(2,2,(m-1)*2+2)
    slice = S.smb_anom2015(:,:,idxEnd);
    scatter(S.x, S.y, 8, slice(:), 'filled');
    clim([-2,2])
    colormap(vik); colorbar; axis equal tight
    title(sprintf('%s SMB anomaly, %d', S.model, S.year(idxEnd)))
end
sgtitle('SMB anomaly relative to 2015 (m ice yr^{-1})')

%% ============================================================
%  Plot 2: Surface temperature anomaly, 2015 vs 2300, rows = model
%% ============================================================
figure;
for m = 1:numel(models)
    S = models(m).S;
    idxStart = find(S.year == refYear, 1);
    idxEnd   = numel(S.year);

    subplot(2,2,(m-1)*2+1)
    slice = S.ts_anom2015(:,:,idxStart);
    scatter(S.x, S.y, 8, slice(:), 'filled');
    colormap(vik); colorbar; axis equal tight
    title(sprintf('%s Ts anomaly, %d', S.model, S.year(idxStart)))

    subplot(2,2,(m-1)*2+2)
    slice = S.ts_anom2015(:,:,idxEnd);
    scatter(S.x, S.y, 8, slice(:), 'filled');
    colormap(vik); colorbar; axis equal tight
    title(sprintf('%s Ts anomaly, %d', S.model, S.year(idxEnd)))
end
sgtitle('Surface temperature anomaly relative to 2015 (K)')

%% ============================================================
%  Local function: netCDF time -> calendar year
%% ============================================================
function yrs = local_time2year(fname)
    t      = ncread(fname, 'time');
    units  = ncreadatt(fname, 'time', 'units');   % e.g. 'days since 1995-1-1 0:0:0'
    tok    = regexp(units, 'since\s+([\d\-]+)', 'tokens');
    refdate = datetime(tok{1}{1}, 'InputFormat', 'yyyy-M-d');
    dates   = refdate + days(t);
    yrs     = year(dates);
end
