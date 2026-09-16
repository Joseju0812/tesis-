function fig = plotMontage2D_only(opts)
% Plot limpio: fuentes, detectores, canales y pares (sin señales)

arguments
    opts struct
end

% ----------- Required -----------
req = {'datos','channels'};
for k=1:numel(req)
    if ~isfield(opts, req{k})
        error('Falta opts.%s', req{k});
    end
end

datos = opts.datos;
channels = opts.channels;

% ----------- Opcionales -----------
pairsValidos = optget(opts,'pairsValidos',[]);
shortIdx     = optget(opts,'shortIdx',[]);
showPairs    = optget(opts,'showPairs',true);
markShort    = optget(opts,'markShortChannels',true);

% ----------- Posiciones -----------
sp = datos.nirs.probe.sourcePos3D;
dp = datos.nirs.probe.detectorPos3D;

if size(sp,1)==3, sp = sp'; end
if size(dp,1)==3, dp = dp'; end

% Midpoints (canales)
nCh = numel(channels);
mids = nan(nCh,3);

for i = 1:nCh
    s = channels(i).sourceIdx;
    d = channels(i).detectorIdx;
    if s>0 && d>0
        mids(i,:) = (sp(s,:) + dp(d,:))/2;
    end
end

xy = mids(:,1:2);

% ----------- Figura -----------
fig = figure('Color','w');
ax = axes(fig); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

title(ax,'Montage 2D (Sources–Detectors–Channels)');
xlabel(ax,'X'); ylabel(ax,'Y');

% ----------- Fuentes -----------
scatter(ax, sp(:,1), sp(:,2), 60, 'r','filled');
text(sp(:,1)+1, sp(:,2), "S"+(1:size(sp,1))', 'Color','r');

% ----------- Detectores -----------
scatter(ax, dp(:,1), dp(:,2), 60, 'b','filled');
text(dp(:,1)+1, dp(:,2), "D"+(1:size(dp,1))', 'Color','b');

% ----------- Canales (midpoints) -----------
valid = ~any(isnan(xy),2);
scatter(ax, xy(valid,1), xy(valid,2), 40, 'k','filled');

% ----------- Pares short-long -----------
% if showPairs && ~isempty(pairsValidos)
%     for p = 1:numel(pairsValidos)
%         L = pairsValidos(p).long_idx;
%         S = pairsValidos(p).short_idx;
% 
%         if all([L S] > 0) && all([L S] <= nCh)
%             pL = xy(L,:);
%             pS = xy(S,:);
% 
%             if ~any(isnan([pL pS]))
%                 plot(ax,[pL(1) pS(1)], [pL(2) pS(2)], 'k-','LineWidth',1.2);
%             end
%         end
%     end
% end

for ch = 1:nCh
    s = channels(ch).sourceIdx;
    d = channels(ch).detectorIdx;

    if s>0 && d>0
        pS = sp(s,1:2); % fuente
        pD = dp(d,1:2); % detector

        plot(ax, [pS(1) pD(1)], [pS(2) pD(2)], ...
            'b-', 'LineWidth', 1.5); % azul
    end
end
% ----------- Marcar canales cortos -----------
if markShort && ~isempty(shortIdx)
    idx = shortIdx(shortIdx <= nCh);
    idx = idx(~any(isnan(xy(idx,:)),2));

    plot(ax, xy(idx,1), xy(idx,2), 'yo', ...
        'MarkerSize',10,'LineWidth',1.5);
end

end

% ----------- helper -----------
function out = optget(s, field, def)
    if isfield(s, field) && ~isempty(s.(field))
        out = s.(field);
    else
        out = def;
    end
end