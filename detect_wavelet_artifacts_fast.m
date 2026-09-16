function result = detect_wavelet_artifacts_fast(signals, fs, blocks, varargin)

% -------- Parse inputs
p = inputParser;
addRequired(p,'signals');
addRequired(p,'fs',@(x)isnumeric(x)&&isscalar(x)&&x>0);
addRequired(p,'blocks');
addParameter(p,'scale',9,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'threshold',43,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'scales',1:30,@isvector);
addParameter(p,'padSec',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'dim','channelsFirst',@(s) ismember(s,{'channelsFirst','channelsLast'}));
addParameter(p,'saveNormCWT',false,@islogical);
addParameter(p,'plotBlock',false,@islogical);
addParameter(p,'plotChannel',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
parse(p,signals,fs,blocks,varargin{:});

sc = p.Results.scale;
TH = p.Results.threshold;
scales = p.Results.scales(:)';
padSec = p.Results.padSec;
dim = p.Results.dim;
saveNormCWT = p.Results.saveNormCWT;
plotBlock = p.Results.plotBlock;
plotChannel = p.Results.plotChannel;

% -------- Prepara señales como channels x samples
sig = p.Results.signals;
if isvector(sig)
    sig = reshape(sig,1,[]);
else
    if strcmp(dim,'channelsLast') && size(sig,1) ~= 1
        sig = sig'; % samples x channels -> channels x samples
    end
end
[nch, N] = size(sig);

% -------- Bloques (conversión a índices)
if isscalar(blocks) && blocks>0
    blockLenS = blocks;
    blockLenSamples = round(blockLenS * fs);
    nBlocks = floor(N / blockLenSamples);
    blocksIdx = zeros(nBlocks,2);
    for k=1:nBlocks
        si = (k-1)*blockLenSamples + 1;
        ei = si + blockLenSamples - 1;
        blocksIdx(k,:) = [si, ei];
    end
else
    blocksIdx = double(blocks);
    if size(blocksIdx,2)~=2, error('blocks must be Px2 indices or scalar seconds'); end
    nBlocks = size(blocksIdx,1);
end

% -------- Precompute kernels (Haar discrete-like) por escala
nScales = length(scales);
kernels = cell(1,nScales);
kernelLens = zeros(1,nScales);
for si = 1:nScales
    s = scales(si);
    kLen = max(2, round(s));          % kernel length en muestras (al menos 2)
    if mod(kLen,2)~=0, kLen = kLen+1; end
    hlen = kLen/2;
    kernel = [ones(1, hlen), -ones(1, hlen)];
    kernel = kernel / sqrt(s);        % normalización por sqrt(scale)
    kernels{si} = kernel;
    kernelLens(si) = kLen;
end

% -------- Result containers
artifactMask = false(nBlocks, nch);
maxNormValue = nan(nBlocks, nch);
if saveNormCWT
    normCWT = cell(nBlocks, nch);
else
    normCWT = [];
end

padSamps = round(padSec * fs);

% -------- Computar conv por canal y escala (una vez) y luego responder bloques

% máscara por muestra: channels x samples
sampleMask = false(nch, N);

% posición del máximo normalizado en muestras dentro de cada bloque y canal
maxPos = nan(nBlocks, nch);

for ch = 1:nch
    x = sig(ch,:);
    % Prealocar coefMat (scales x N) - usar 'single' si N grande para ahorrar memoria
    coefMat = zeros(nScales, N, 'like', x);
    for si = 1:nScales
        k = kernels{si};
        % conv 'same' para mantener longitud N
        coefMat(si,:) = conv(x, k, 'same');
    end
    b2full = abs(coefMat).^2; % scales x N
    
    % Para cada bloque extraer la porción interior (sin pad) y calcular mediana por escala
    for b = 1:nBlocks
        si = max(1, blocksIdx(b,1) - padSamps);
        ei = min(N, blocksIdx(b,2) + padSamps);
        innerStart = blocksIdx(b,1);
        innerEnd   = blocksIdx(b,2);
        % índices relativos en b2full (misma matriz)
        relInnerStart = innerStart;
        relInnerEnd = innerEnd;
        % medianas por escala sobre la parte interior (tal como Sato)
        med = median( b2full(:, relInnerStart:relInnerEnd), 2 );
        med(med==0) = eps;
        % Normalización (solo para el interior si queremos; aquí normalizamos con med del interior)
        wnorm_block = bsxfun(@rdivide, b2full(:, relInnerStart:relInnerEnd), med);
        % Índice de la escala seleccionada
        [~, sc_idx] = min(abs(scales - sc));
        % máximo en la escala elegida dentro del bloque interior
        maxval = max( wnorm_block(sc_idx, :), [], 2 );
        maxNormValue(b,ch) = maxval;
        artifactMask(b,ch) = (maxval > TH);
        
        %%
        % localizar índice del máximo dentro del bloque (en la escala sc_idx)
        [~, idxMax] = max(wnorm_block(sc_idx,:));         % índice relativo dentro del bloque interior
        tmax_global = (relInnerStart + idxMax - 1);       % índice de muestra global (1..N)
        maxPos(b,ch) = tmax_global;
        
        % si ese bloque fue marcado como artefacto, marcar todas las muestras del bloque
        if artifactMask(b,ch)
            sIdx = relInnerStart:relInnerEnd;            % rango de muestras del bloque interior
            sampleMask(ch, sIdx) = true;                 % marcar en la máscara por muestra
        end


        if saveNormCWT
            normCWT{b,ch} = wnorm_block; % scales x time_interior
        end
        % Plot opcional (por bloque y por canal)
        if plotBlock && ch == plotChannel
            t_block = (relInnerStart:relInnerEnd)/fs;
            figure;
            imagesc(t_block, scales, wnorm_block);
            set(gca,'YDir','normal');
            xlabel('Tiempo [s]'); ylabel('Escala');
            title(sprintf('Escalograma - canal %d - bloque %d (max=%.2f)', ch, b, maxval));
            colorbar;
            hold on;
            yline(scales(sc_idx),'w--','LineWidth',1.2); % escala usada
            % marcar tiempo del maximo dentro del bloque
            [~, idxMax] = max(wnorm_block(sc_idx,:));
            tmax_global = (relInnerStart + idxMax -1)/fs;
            plot(tmax_global, scales(sc_idx),'wo','MarkerFaceColor','w');
        end
    end
    % Liberar memoria de coefMat / b2full antes del siguiente canal (opcional)
    clear coefMat b2full
end

% -------- Package result
result.artifactMask = artifactMask;
result.maxNormValue = maxNormValue;
result.normCWT = normCWT;
result.params.scale = sc;
result.params.threshold = TH;
result.params.scales = scales;
result.params.fs = fs;
result.params.padSec = padSec;
result.params.blocks = blocksIdx;
result.sampleMask = sampleMask;            % nch x N logical (true = artefacto)
result.sampleMaskAny = any(sampleMask,1);  % 1 x N logical (artefacto en cualquier canal)
result.maxPos = maxPos;   

end
