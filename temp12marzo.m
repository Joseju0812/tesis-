clc; clear; close all;

%% ---------------- 1) Paths y carga ----------------
rutaDatos = fullfile('C:\Users\josue\OneDrive\Documentos\codigos de modlar\copia de datos\01_Pre\Subj001_Sess001.snirf');

addpath('C:\Users\josue\OneDrive\Documentos\modular');
addpath('C:\Users\josue\OneDrive\Documentos\codigos\funciones');
addpath('C:\Users\josue\OneDrive\Documentos\codigos de modlar')

datos = cargarDatos(rutaDatos,'/');

t = datos.nirs.data1.time;
fs = 1/mean(diff(t));

I_raw = datos.nirs.data1.dataTimeSeries;

[nTime_total, nMeas] = size(I_raw);

%% ---------------- 2) Preprocesado inicial ----------------

I_tddr = TDDR_(I_raw, fs);

medianWin_s = 0.4;
medianWinSamples = max(1, round(medianWin_s * fs));

I_med = movmedian(I_tddr, medianWinSamples,1,'Endpoints','shrink');

%% ---------------- 3) Baseline ----------------

baseline_seconds = 5;
n0 = max(1, round(baseline_seconds * fs));

I0 = median(I_med(1:n0,:),1);

I0(I0<=0)=eps;
I_med(I_med<=0)=eps;

%% ---------------- 4) Construir pares y distancias ----------------

wavelengthsVec = datos.nirs.probe.wavelengths(:)';

shortThresh_mm = 10;
maxPairDist_mm = 40;

[pairs, rho_chan_mm, channels, measArr_fileorder] = ...
    parespordistanciayrho(datos, shortThresh_mm, maxPairDist_mm,'Units','mm','Plot',false);

rho_chan_mm = rho_chan_mm(:)';
rho_chan_cm = rho_chan_mm ./ 10;

shortIdx = find(rho_chan_mm <= shortThresh_mm);
longIdx  = find(rho_chan_mm > shortThresh_mm);

nChan = numel(channels);
nWl   = numel(wavelengthsVec);
nTime = size(I_med,1);

fprintf('Canales totales: %d. Shorts: %d ; Longs: %d\n',numel(rho_chan_mm),numel(shortIdx),numel(longIdx));

%% ---------------- 5) Convertir a OD ----------------

OD = -log(I_med ./ I0);

%% ---------------- 6) Filtro Bandpass ----------------

bpHP = 0.01;
bpLP = 0.2;

filterOrder = 4;
nyq = fs/2;

Wn = [bpHP bpLP] / nyq;

[b,a] = butter(filterOrder,Wn,'bandpass');

OD_filt = filtfilt(b,a,OD);

%% ---------------- 7) Detrend ----------------

OD_detr = detrend(OD_filt);

%% ================= SCR MULTI-SHORT EN OD =================

disp('Aplicando Short Channel Regression (multi-short) en OD...')

OD_corr = OD_detr;

N = size(OD_detr,1);
onesN = ones(N,1);

short_ml = [];

for s = 1:length(shortIdx)
    short_ml = [short_ml channels(shortIdx(s)).ml_indices];
end

short_ml = unique(short_ml);

Xshort = OD_detr(:,short_ml);

for ch = 1:nChan

    if ismember(ch,shortIdx)
        continue
    end

    rows_long = channels(ch).ml_indices;

    for w = 1:nWl

        rowL = rows_long(w);

        yL = OD_detr(:,rowL);

        X = [Xshort onesN];

        beta = X \ yL;

        beta_short = beta(1:size(Xshort,2));

        y_corr = yL - Xshort * beta_short;

        OD_corr(:,rowL) = y_corr;

    end
end

disp('SCR finalizado.')

%% ---------------- 8) Reconstruir intensidad ----------------

I_scr = I0 .* exp(-OD_corr);

I_proc = I_scr;

%% ---------------- 9) Reordenar intensidad a 3D ----------------

I_3D = zeros(nChan,nTime,nWl);
I0_3D = zeros(nChan,nWl);

for ch = 1:nChan

    rows = channels(ch).ml_indices;

    if isempty(rows)
        warning('Canal %d sin indices',ch)
        continue
    end

    wlIdxs = arrayfun(@(r) measArr_fileorder(r).wavelengthIndex, rows);

    [~,ord] = sort(wlIdxs,'ascend');

    rows_sorted = rows(ord);

    if numel(rows_sorted) ~= nWl
        if numel(rows_sorted) < nWl
            rows_sorted = [rows_sorted repmat(rows_sorted(end),1,nWl-numel(rows_sorted))];
        end
    end

    for w = 1:nWl
        I_3D(ch,:,w) = I_proc(:,rows_sorted(w));
        I0_3D(ch,w) = I0(rows_sorted(w));
    end

end

%% ---------------- 10) Matriz de extinciones ----------------

eps_mat = GetExtinctions(wavelengthsVec);

%% ---------------- 11) Calibración CW-fNIRS ----------------

[DeltaC_cal, S_est, G_est, OD_cal] = cwfnirs_calibrate( ...
    I_3D, I0_3D, eps_mat, rho_chan_cm, ...
    'Detrend', true, ...
    'RegLambda',1e-4, ...
    'BaselineSamples',n0);

%% ---------------- 12) Ganancia ----------------

gain_linear = exp(-mean(G_est,2));

figure
plot(gain_linear,'LineWidth',1.2)
title('Ganancia efectiva por canal')
xlabel('Canal')
ylabel('Ganancia')
grid on

hold on

for k = 1:numel(shortIdx)
    xline(shortIdx(k),'y','LineWidth',2)
end

%% ---------------- 13) Reconstruir OD calibrada ----------------

OD_recon = zeros(nTime,nMeas);

for ch = 1:nChan

    rows = channels(ch).ml_indices;

    wlIdxs = arrayfun(@(r) measArr_fileorder(r).wavelengthIndex, rows);

    [~,ord] = sort(wlIdxs,'ascend');

    rows_sorted = rows(ord);

    for w = 1:min(nWl,numel(rows_sorted))

        OD_recon(:,rows_sorted(w)) = reshape(OD_cal(ch,:,w),[nTime 1]);

    end
end

%% ---------------- 14) Obtener concentraciones ----------------

conc_HbO_cal = reshape(DeltaC_cal(1,:,:),[nTime nChan]) * 1e6;
conc_HbR_cal = reshape(DeltaC_cal(2,:,:),[nTime nChan]) * 1e6;

conc_HbT_cal = conc_HbO_cal + conc_HbR_cal;

%% ---------------- 15) Corrección Sato Wavelet ----------------

disp('Detectando artefactos con wavelets...')

conc_HbO_orig = conc_HbO_cal;
conc_HbR_orig = conc_HbR_cal;
conc_HbT_orig = conc_HbT_cal;

escala_sato = 9;
umbral_sato = 43;
block_dur = 1;

saltos_guardados = cell(nChan,1);

for ch = 1:nChan

    senal_hbt = conc_HbT_orig(:,ch);

    res = detect_wavelet_artifacts_fast(senal_hbt,fs,block_dur,...
        'scale',escala_sato,'threshold',umbral_sato,...
        'scales',1:20,'padSec',1);

    blocksWithArtifact = find(res.artifactMask(:,1));

    if ~isempty(blocksWithArtifact)

        jumps_idx = floor(res.maxPos(blocksWithArtifact,1))+1;

        jumps_idx(isnan(jumps_idx))=[];
        jumps_idx = jumps_idx(jumps_idx < length(senal_hbt));

        if ~isempty(jumps_idx)

            saltos_guardados{ch}=jumps_idx;

            conc_HbO_cal(:,ch) = remove_step_offsets(conc_HbO_orig(:,ch),jumps_idx);
            conc_HbR_cal(:,ch) = remove_step_offsets(conc_HbR_orig(:,ch),jumps_idx);

            conc_HbT_cal(:,ch) = conc_HbO_cal(:,ch) + conc_HbR_cal(:,ch);

        end
    end

end

disp('Corrección Sato finalizada.')

%% ---------------- 16) GLM ----------------

nsample = nTime;

if isfield(datos.nirs,'stim1') && isfield(datos.nirs,'stim2')

    stim_onset = datos.nirs.stim1.data(:,1);
    stim_duration = datos.nirs.stim2.data(:,1) - datos.nirs.stim1.data(:,1);

    onsets = round(stim_onset*fs);
    durations = max(1,round(stim_duration*fs));

    st = zeros(nsample,1);

    for k=1:numel(onsets)

        on = min(max(1,onsets(k)),nsample);
        off = min(nsample,on+durations(k));

        st(on:off)=1;

    end

else

    st=zeros(nsample,1);

end

HRF = respuestahemodinamica(30,fs);

reg = conv(st,HRF,'same');

dreg = [0; diff(reg)];

X = [reg dreg ones(nsample,1)];

Y = [conc_HbO_cal conc_HbR_cal];

B = pinv(X)*Y;

nChan = size(conc_HbO_cal,2);

beta_HbO = B(1,1:nChan);
beta_HbR = B(1,nChan+1:end);

%% ---------------- 17) Visualización ejemplo ----------------

ch_ej = 13;

figure

subplot(2,1,1)
plot(t,conc_HbO_cal (:,ch_ej) + 0.1,'r')
hold on
plot(t,conc_HbR_cal(:,ch_ej),'b')
title('Concentraciones HbO HbR')
% 
% subplot(2,1,2)
% bar([beta_HbO(ch_ej) beta_HbR(ch_ej)])
% title('Betas GLM')


% %% 11.5) t-test sobre los betas del GLM
% 
% N = size(X,1);     % muestras
% p = size(X,2);     % regresores
% 
% % Varianza residual por canal
% sigma2 = sum(Res.^2) ./ (N - p);     % [1 x 2*nChan]
% 
% % (X'X)^-1
% XtX_inv = inv(X' * X);
% 
% % error estándar del beta HRF (primer regresor)
% se_beta = sqrt( sigma2 .* XtX_inv(1,1) );   % [1 x 2*nChan]
% 
% % t-statistic
% t_stat = B(1,:) ./ se_beta;
% 
% % p-value bilateral
% p_val = 2 * (1 - tcdf(abs(t_stat), N - p));
% 
% % separar HbO y HbR
% t_HbO = t_stat(1:nChan);
% t_HbR = t_stat(nChan+1:end);
% 
% p_HbO = p_val(1:nChan);
% p_HbR = p_val(nChan+1:end);
% 
% % canales significativos
% alpha = 0.05;
% 
% sig_HbO = p_HbO < alpha;
% sig_HbR = p_HbR < alpha;
% 
% fprintf('Canales significativos HbO: %d de %d\n', sum(sig_HbO), nChan);
% fprintf('Canales significativos HbR: %d de %d\n', sum(sig_HbR), nChan);
% figure
% subplot(2,1,1)
% stem(t_HbO,'filled')
% hold on
% yline(2,'r--')
% yline(-2,'r--')
% title('t-stat HbO')
% xlabel('Canal')
% ylabel('t')
% grid on
% 
% subplot(2,1,2)
% stem(t_HbR,'filled')
% hold on
% yline(2,'r--')
% yline(-2,'r--')
% title('t-stat HbR')
% xlabel('Canal')
% ylabel('t')
% grid on

%%
%% ---------------- 18) Estadística GLM ----------------

disp('Calculando estadística GLM...')

Y_hat = X * B;
Res   = Y - Y_hat;

[nSample,nReg] = size(X);
nSignals = size(Y,2);

t_stats = zeros(size(B));
p_values = zeros(size(B));

XtX_inv = pinv(X'*X);

for s = 1:nSignals
    
    sigma2 = sum(Res(:,s).^2) / (nSample - nReg);
    
    CovB = sigma2 * XtX_inv;
    
    se = sqrt(diag(CovB));
    
    t_stats(:,s) = B(:,s) ./ se;
    
end

df = nSample - nReg;

p_values = 2 * (1 - tcdf(abs(t_stats),df));
p_HbO = p_values(1,1:nChan);
p_HbR = p_values(1,nChan+1:end);

t_HbO = t_stats(1,1:nChan);
t_HbR = t_stats(1,nChan+1:end);
figure

subplot(2,1,1)
bar(p_HbO)
hold on
yline(0.05,'r--','LineWidth',2)
title('p-values HbO')
xlabel('Canal')
ylabel('p')

subplot(2,1,2)
bar(p_HbR)
hold on
yline(0.05,'r--','LineWidth',2)
title('p-values HbR')
xlabel('Canal')
ylabel('p')

%%

%% ---------------- 20) Coordenadas de canales ----------------
%% ---------------- 20) Coordenadas de canales ----------------
% 
% srcPos = datos.nirs.probe.sourcePos3D;
% detPos = datos.nirs.probe.detectorPos3D;
% 
% chanPos = zeros(nChan,3);
% 
% for ch = 1:nChan
% 
%     rows = channels(ch).ml_indices;
% 
%     if isempty(rows)
%         continue
%     end
% 
%     % tomar la primera longitud de onda del canal
%     r = rows(1);
% 
%     src = measArr_fileorder(r).sourceIndex;
%     det = measArr_fileorder(r).detectorIndex;
% 
%     chanPos(ch,:) = (srcPos(src,:) + detPos(det,:)) / 2;
% 
% end
% 
% x = chanPos(:,1);
% y = chanPos(:,2);
% z = chanPos(:,3);
% 
% activation = -log10(p_HbO);
% activation(p_HbO > 0.05) = 0;
% figure
% 
% scatter(x,y,140,activation,'filled')
% 
% colormap(jet)
% colorbar
% 
% title('fNIRS Brain Activation Map')
% xlabel('X')
% ylabel('Y')
% 
% axis equal
% grid on
% xi = linspace(min(x),max(x),100);
% yi = linspace(min(y),max(y),100);
% 
% [XI,YI] = meshgrid(xi,yi);
% 
% ZI = griddata(x,y,activation,XI,YI,'cubic');
% 
% figure
% 
% imagesc(xi,yi,ZI)
% set(gca,'YDir','normal')
% 
% hold on
% scatter(x,y,60,'k','filled')
% 
% colormap(jet)
% colorbar
% 
%  title('Interpolated fNIRS Activation Map')
%  plot_fnirs_brainmap(datos,channels,measArr_fileorder,p_HbO,shortIdx)
% %%

%% ---------------- Brain activation map (t-statistics) ----------------

srcPos = datos.nirs.probe.sourcePos3D;
detPos = datos.nirs.probe.detectorPos3D;

chanPos = zeros(nChan,3);

for ch = 1:nChan
    
    rows = channels(ch).ml_indices;
    
    if isempty(rows)
        continue
    end
    
    r = rows(1);
    
    src = measArr_fileorder(r).sourceIndex;
    det = measArr_fileorder(r).detectorIndex;
    
    chanPos(ch,:) = (srcPos(src,:) + detPos(det,:))/2;
    
end

x = chanPos(:,1);
y = chanPos(:,2);


% usar t-stats como activación

activation = t_HbR;

% interpolación espacial

xi = linspace(min(x),max(x),120);
yi = linspace(min(y),max(y),120);

[XI,YI] = meshgrid(xi,yi);

ZI = griddata(x,y,activation,XI,YI,'cubic');

% graficar mapa

figure

imagesc(xi,yi,ZI)
set(gca,'YDir','normal')

hold on

colormap(jet)
colorbar

% dibujar canales

scatter(x,y,80,'k','filled')

% canales significativos

sigIdx = find(abs(t_HbO) > 2);

scatter(x(sigIdx),y(sigIdx),150,'w','LineWidth',2)

% fuentes

scatter(srcPos(:,1),srcPos(:,2),200,'r','filled')

% detectores

scatter(detPos(:,1),detPos(:,2),200,'b','filled')

title('fNIRS Brain Activation Map (t-statistics)')

xlabel('X')
ylabel('Y')

axis equal
grid on

%%
% filtrado espacial, se calcula el promedio a lo largo de todos loa canales
% lo que es la rescpuesta global, se substrae del canal el promedio,
% aplicado de la reconstruccion y el SCR, pero antes de cualquier filtrado
% temporal 
%
% un umbral de los coeficientes de variacion
%
%filtros de savitzki golay, despues de interpolacion por splines cubicos 
% algoritmos de correccion de saltos
% 
% tomar parametros para cada sujeto no crea un sesgo 


OD = -log(I_med ./ I0);
OD_filt = filtfilt(b,a,OD);
OD_detr = detrend(OD_filt);

% ---------- Preparar measArr si no existe (recuperarlo de datos) ----------
if ~exist('measArr','var') || isempty(measArr)
    fn = fieldnames(datos.nirs.data1);
    mlF = fn(startsWith(fn,'measurementList'));
    measArr = struct([]);
    for ii=1:numel(mlF)
        e = datos.nirs.data1.(mlF{ii});
        if isempty(e), continue; end
        if isstruct(e) && numel(e)>1
            measArr = [measArr; e(:)]; %#ok<AGROW>
        else
            measArr = [measArr; e]; %#ok<AGROW>
        end
    end
    fprintf('Reconstruido measArr desde datos -> %d entradas\n', numel(measArr));
end

% ---------- Asegurar shortIdx (si no lo tienes) ----------
if ~exist('shortIdx','var') || isempty(shortIdx)
    try
        shortThresh_mm = 10;
        rho_chan_mm = rho_chan_mm; % asume ya calculado en tu pipeline
        shortIdx = find(rho_chan_mm <= shortThresh_mm);
    catch
        shortIdx = [];
    end
end

% ---------- Lista de procedimientos que queremos graficar ----------
procList = {'i_raw','i_med','od_cal','od_recon','conc_cal','gain_linear'};

outFolder = fullfile(pwd,'plots_all'); 
if ~exist(outFolder,'dir'), mkdir(outFolder); end

for p = 1:numel(procList)
    procName = procList{p};
    opts = struct();
    opts.datos = datos;
    opts.channels = channels;
    opts.measArr = measArr;
    opts.t = t;
    opts.proc = procName;
    opts.channelList = 1:numel(channels);
    opts.insetSize = 0.090;
    opts.showPairs = true;
    opts.pairsValidos = pairs;    % si tienes 'pairs'
    opts.shortIdx = shortIdx;
    opts.markShortChannels = true;
    opts.targetP2P = 18;
    opts.scalingMode = 'local';
    opts.annotateScale = true;
    opts.maxGain = 500;

    % pasar datos según proc
    switch procName
        case 'i_raw'
            opts.I_raw = I_raw;
        case 'i_med'
            opts.I_med = I_med;
        case 'od'
            opts.OD = OD;
        case 'od_filt'
            opts.OD_filt = OD_filt;
        case 'od_detr'
            opts.OD_detr = OD_detr;
        case 'od_cal'
            % si tienes OD_cal en formato ch x t x wl pásalo
            if exist('OD_cal','var'), opts.OD_cal = OD_cal; end
        case 'od_recon'
            if exist('OD_recon','var'), opts.OD_recon = OD_recon; end
        case 'conc'
            if exist('conc_HbO','var'), opts.conc_HbO = conc_HbO; end
            if exist('conc_HbR','var'), opts.conc_HbR = conc_HbR; end
        case 'conc_cal'
            if exist('conc_HbO_cal','var'), opts.conc_HbO_cal = conc_HbO_cal; end
            if exist('conc_HbR_cal','var'), opts.conc_HbR_cal = conc_HbR_cal; end
        case 'gain_linear'
            if exist('gain_linear','var'), opts.gain_linear = gain_linear; end
    end

        % fig = plotMontageWithInserts_v8(opts);

end
%%

fig = plotMontage2D_only(opts);
%% ---------------- Función local ----------------

function datos = cargarDatos(filepath,nodo)

dataInfo = h5info(filepath);

datos = leernodo(dataInfo,filepath,nodo);

end