% function [DeltaC, S_est, G_est, OD_cal] = cwfnirs_calibrate(I_3D, I0_3D, eps_mat, rho_chan_cm, varargin)
% % CWFNIRS_CALIBRATE Calibración CW-fNIRS basada en OD = E*S*dc + G.
% %   I_3D      : [nChan x nTime x nWl] intensidades
% %   I0_3D     : [nChan x nWl] baseline por canal y lambda
% %   eps_mat   : [nWl x nChrom] coeficientes de extinción (típico nChrom=2)
% %   rho_chan_cm : [1 x nChan] distancia fuente-detector en cm
% %
% % Opciones:
% %   'Detrend'         : true/false (default true)
% %   'RegLambda'       : regularización ridge para inversión (default 1e-4)
% %   'BaselineSamples' : nº muestras para estimar baseline en OD (default [] => todo)
% %
% % Salidas:
% %   DeltaC : [nChrom x nTime x nChan]
% %   S_est  : [nChan x 1] factor de trayectoria efectivo (escala tipo DPF)
% %   G_est  : [nChan x nWl] offset por canal/lambda
% %   OD_cal : [nChan x nTime x nWl] OD calibrada (offset removido)
% 
% p = inputParser;
% p.addParameter('Detrend', true, @(x) islogical(x) || isnumeric(x));
% p.addParameter('RegLambda', 1e-4, @(x) isnumeric(x) && isscalar(x) && x>=0);
% p.addParameter('BaselineSamples', [], @(x) isempty(x) || (isscalar(x) && x>=1));
% p.parse(varargin{:});
% opts = p.Results;
% 
% [nChan, nTime, nWl] = size(I_3D);
% nChrom = size(eps_mat, 2);
% 
% rho_chan_cm = rho_chan_cm(:);
% if isempty(opts.BaselineSamples)
%     baseIdx = 1:nTime;
% else
%     baseIdx = 1:min(nTime, round(opts.BaselineSamples));
% end
% 
% DeltaC = zeros(nChrom, nTime, nChan);
% S_est = zeros(nChan,1);
% G_est = zeros(nChan,nWl);
% OD_cal = zeros(nChan,nTime,nWl);
% 
% for ch = 1:nChan
%     I_ch = squeeze(I_3D(ch,:,:)); % [nTime x nWl]
%     I0_ch = I0_3D(ch,:);
% 
%     I_ch(I_ch <= 0) = eps;
%     I0_ch(I0_ch <= 0) = eps;
% 
%     OD = -log(bsxfun(@rdivide, I_ch, I0_ch)); % [nTime x nWl]
% 
%     if opts.Detrend
%         OD = detrend(OD);
%     end
% 
%     G = mean(OD(baseIdx,:), 1);
%     OD0 = bsxfun(@minus, OD, G);
% 
%     A = eps_mat * rho_chan_cm(ch);   % [nWl x nChrom], falta factor S
%     dc_base = pinv(A) * OD0(baseIdx,:).';
%     OD_base_pred = (A * dc_base).';
% 
%     num = sum(sum(OD0(baseIdx,:) .* OD_base_pred));
%     den = sum(sum(OD_base_pred.^2)) + eps;
%     S = num / den;
% 
%     if ~isfinite(S) || S <= 0
%         S = 1;
%     end
% 
%     A_cal = A * S;
%     AtA = A_cal.' * A_cal;
%     reg = opts.RegLambda * mean(diag(AtA));
%     W = (AtA + reg * eye(nChrom)) \ A_cal.';
% 
%     dc_t = (W * OD0.').';
% 
%     DeltaC(:,:,ch) = dc_t.';
%     S_est(ch) = S;
%     G_est(ch,:) = G;
%     OD_cal(ch,:,:) = reshape(OD0, [1, nTime, nWl]);
% end
% end
% % 
function [DeltaC, S_est, G_est, OD_cal] = cwfnirs_calibrate(I_3D, I0_3D, eps_mat, rho_chan_cm, varargin)
% CWFNIRS_CALIBRATE_FAST Versión optimizada y vectorizada del Algoritmo 2.


    detrend_flag = true; 
    reg_lambda   = 1e-4;
    base_samps   = 5; 
    % Procesar varargin manualmente si es necesario (simplificado)
    if ~isempty(varargin)
        for k = 1:2:length(varargin)
            switch lower(varargin{k})
                case 'detrend', detrend_flag = varargin{k+1};
                case 'reglambda', reg_lambda = varargin{k+1};
                case 'baselinesamples', base_samps = varargin{k+1};
            end
        end
    end

    % --- 2. Preparación de Datos  ---
    [nChan, nTime, nWl] = size(I_3D);
    nChrom = size(eps_mat, 2);
    rho_chan_cm = rho_chan_cm(:); 

    % Definir índices de baseline
    if isempty(base_samps), baseIdx = 1:nTime; 
    else, baseIdx = 1:min(nTime, round(base_samps)); end

    I_perm  = permute(I_3D, [2, 3, 1]);    % [Time, Wl, Chan]
    I0_perm = permute(I0_3D, [3, 2, 1]);   % [1, Wl, Chan] (para broadcasting)

    % Evitar log(0)
    I_perm(I_perm <= 0) = eps;
    I0_perm(I0_perm <= 0) = eps;

    % --- 3. Cálculo de OD y Pre-procesamiento (Todo en bloque) ---
   
    OD_all = -log(I_perm ./ I0_perm);      % [Time, Wl, Chan]

    if detrend_flag
        OD_all = detrend(OD_all); 
    end

    % Baseline G (promedio en dim 1) y centrado
    G_all = mean(OD_all(baseIdx, :, :), 1); % [1, Wl, Chan]
    OD0_all = OD_all - G_all;               % [Time, Wl, Chan]

    % --- 4. Inversión del Modelo
    DeltaC = zeros(nChrom, nTime, nChan);
    S_est  = ones(nChan, 1);

    % Matriz identidad pre-calculada para regularización
    EyeChrom = eye(nChrom); 

    for ch = 1:nChan
    
        % reshape fuerza a 2D [Time x Wl] eliminando dimensiones singleton
        OD0_ch = reshape(OD0_all(:, :, ch), nTime, nWl); 

        DPF = 6.25;

        % Modelo físico completo
        A = eps_mat * (DPF * rho_chan_cm(ch));
        
        AtA = A' * A;
        reg = reg_lambda * mean(diag(AtA));
        W = (AtA + reg * EyeChrom) \ A';
        
        DeltaC(:, :, ch) = W * OD0_ch';
        S_est(ch) = DPF;   % solo guardas el valor fijo

        % B. Inversión final usando S estimado
        A_cal = A * DPF;                         % Ajustar modelo
        AtA_cal = A_cal' * A_cal;
        reg2 = reg_lambda * mean(diag(AtA_cal));

        % Solución mínimos cuadrados para todo el tiempo
        % (AtA + lambda*I)^-1 * A' * OD'
        W_final = (AtA_cal + reg2 * EyeChrom) \ A_cal'; 

        % Calcular concentraciones finales [nChrom x nTime]
        DeltaC(:, :, ch) = W_final * OD0_ch'; 
    end

    % --- 5. Formatear Salidas Restantes ---+
    % G_est: [nChan x nWl]
    G_est = permute(G_all, [3, 2, 1]); 
    % OD_cal: regresar a [nChan x nTime x nWl]
    OD_cal = permute(OD0_all, [3, 1, 2]);

end