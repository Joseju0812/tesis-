

function [pairsValidos, rho_chan_mm, channels,measArr] = parespordistanciayrho(datos, shortThresh_mm, maxPairDist_mm, varargin)
% ----- parse inputs -----
p = inputParser;
addRequired(p,'datos');
addRequired(p,'shortThresh_mm',@isnumeric);
addRequired(p,'maxPairDist_mm',@isnumeric);
addParameter(p,'Units','mm',@(x) ismember(x,{'mm','cm','m'}));
addParameter(p,'Plot',false,@islogical);
parse(p, datos, shortThresh_mm, maxPairDist_mm, varargin{:});
units = datos.nirs.metaDataTags.LengthUnit;
doPlot = p.Results.Plot;


% ----- collect measurementList entries -----
fn = fieldnames(datos.nirs.data1);
mlFields = fn(startsWith(fn,'measurementList'));
measArr = struct([]);
for i=1:numel(mlFields)
    entry = datos.nirs.data1.(mlFields{i});
    if isempty(entry), continue; end
    if isstruct(entry) && numel(entry)>1
        measArr = [measArr; entry(:)];
    else
        measArr = [measArr; entry];
    end
end
M_meas = numel(measArr);
if M_meas==0
    error('No measurementList found in datos.nirs.data1');
end

% ----- probe positions normalized to [N x 3] -----
sp = datos.nirs.probe.sourcePos3D;
dp = datos.nirs.probe.detectorPos3D;
if size(sp,1)==3 && size(sp,2)~=3, srcPos = sp'; else srcPos = sp; end
if size(dp,1)==3 && size(dp,2)~=3, detPos = dp'; else detPos = dp; end

% % ----- compute rho per entry (in original units) -----
% srcIdx = nan(M_meas,1); detIdx = nan(M_meas,1);
% for i=1:M_meas
%     srcIdx(i) = measArr(i).sourceIndex;
%     detIdx(i) = measArr(i).detectorIndex;
% end
% rho_entry = nan(M_meas,1);
% for i=1:M_meas
%     s = srcIdx(i); d = detIdx(i);
%     if any([s,d] < 1) || s>size(srcPos,1) || d>size(detPos,1)
%         rho_entry(i) = NaN;
%     else
%         rho_entry(i) = norm(srcPos(s,:) - detPos(d,:)); % same units as srcPos/detPos
%     end
% end
% 
% % ----- convert rho_entry -> mm according to 'Units' param -----
% switch units
%     case 'mm'
%         rho_entry_mm = rho_entry;
%     case 'cm'
%         rho_entry_mm = rho_entry * 10;
%     case 'm'
%         rho_entry_mm = rho_entry * 1000;
%     otherwise
%         rho_entry_mm = rho_entry; % fallback
% end

% ----- compute rho per entry (in original units) ----- 
srcIdx = nan(M_meas,1); detIdx = nan(M_meas,1);
for i=1:M_meas
    srcIdx(i) = measArr(i).sourceIndex;
    detIdx(i) = measArr(i).detectorIndex;
end

rho_entry = nan(M_meas,1);
validEntry = false(M_meas,1);
for i=1:M_meas
    s = srcIdx(i); d = detIdx(i);
    if any([s,d] < 1) || s>size(srcPos,1) || d>size(detPos,1)
        rho_entry(i) = NaN;
        validEntry(i) = false;
    else
        rho_entry(i) = norm(srcPos(s,:) - detPos(d,:)); % still in srcPos units
        validEntry(i) = true;
    end
end

% Remove invalid measurementList entries (they do not map to valid src/det)
if any(~validEntry)
    warning('parespordistanciayrho: %d measurementList entries invalid (src/det out of range). These will be ignored.', sum(~validEntry));
    measArr = measArr(validEntry);
    srcIdx = srcIdx(validEntry);
    detIdx = detIdx(validEntry);
    rho_entry = rho_entry(validEntry);
    M_meas = numel(measArr);
end

% ----- At this point decide scale of srcPos/detPos and convert to mm if needed -----
% estimate typical src-det distance (after removing invalid)
try
    sampleN = min(6, numel(srcIdx));
    samplePairs = 1:sampleN;
    sample_d = mean(arrayfun(@(i) norm(srcPos(srcIdx(i),:)-detPos(detIdx(i),:)), samplePairs),'omitnan');
catch
    sample_d = mean(vecnorm(srcPos - detPos(1:min(end,size(detPos,1)),:),2,2),'omitnan');
end
% heuristics: if average < 1.5 -> meters; if between 1.5 and 15 -> cm; else assume mm
if sample_d < 1.5
    scaleFactor = 1000;   % m -> mm
elseif sample_d < 15
    scaleFactor = 10;     % cm -> mm
else
    scaleFactor = 1;      % mm
end
if scaleFactor ~= 1
    warning('parespordistanciayrho: scaling src/det positions by %g to convert to mm (heuristic).', scaleFactor);
    srcPos = srcPos * scaleFactor;
    detPos = detPos * scaleFactor;
    rho_entry = rho_entry * scaleFactor;
end

% ----- convert rho_entry -> mm according to 'Units' param (defensive) -----
switch units
    case 'mm'
        rho_entry_mm = rho_entry;
    case 'cm'
        rho_entry_mm = rho_entry * 10;
    case 'm'
        rho_entry_mm = rho_entry * 1000;
    otherwise
        rho_entry_mm = rho_entry;
end


% ----- build unique channels per src-det pair and compute rho per channel -----
pairs = [srcIdx detIdx];
[uniqPairs, ~, icPairs] = unique(pairs,'rows','stable');
nChan = size(uniqPairs,1);
channels = repmat(struct('sourceIdx',[],'detectorIdx',[],'canalID',[],'ml_indices',[],'rho_entries_mm',[],'rho_chan_mm',[],'posNom_mm',[]),1,nChan);
for ch=1:nChan
    s = uniqPairs(ch,1); d = uniqPairs(ch,2);
    idxs = find(icPairs==ch);
    channels(ch).sourceIdx = s;
    channels(ch).detectorIdx = d;
    channels(ch).canalID = sprintf('S%02d-D%02d',s,d);
    channels(ch).ml_indices = idxs;
    channels(ch).rho_entries_mm = rho_entry_mm(idxs);
    channels(ch).rho_chan_mm = median(channels(ch).rho_entries_mm,'omitnan');
    % midpoint: note srcPos/detPos are in original units -> convert if needed
    mid = (srcPos(s,:) + detPos(d,:)) / 2;
    switch units
        case 'mm', channels(ch).posNom_mm = mid;
        case 'cm', channels(ch).posNom_mm = mid * 10;
        case 'm',  channels(ch).posNom_mm = mid * 1000;
    end
end
rho_chan_mm = [channels.rho_chan_mm];

% ----- separate short / long by rho_chan_mm -----
shortIdx = find(rho_chan_mm <= shortThresh_mm);
longIdx  = find(rho_chan_mm >  shortThresh_mm);

% ----- many-to-one matching: each long chooses nearest short (shorts can repeat) -----
pairsValidos = struct('long_chan',{}, 'short_chan',{}, 'dist_mm',{}, 'long_idx',{}, 'short_idx',{});
if isempty(longIdx) || isempty(shortIdx)
    % nothing to match
    return;
end
% filter out channels with NaN midpoints
validLong = ~any(isnan(vertcat(channels(longIdx).posNom_mm)),2);
validShort = ~any(isnan(vertcat(channels(shortIdx).posNom_mm)),2);

% map indices back to original indices arrays
longIdx_valid = longIdx(validLong);
shortIdx_valid = shortIdx(validShort);

if isempty(longIdx_valid) || isempty(shortIdx_valid)
    pairsValidos = struct([]);
    return;
end

posLong = vertcat(channels(longIdx_valid).posNom_mm);
posShort = vertcat(channels(shortIdx_valid).posNom_mm);
D = pdist2(posLong, posShort);


k = 1;
for iL = 1:size(D,1)
    [dmin, jShort] = min(D(iL,:));
    if dmin <= maxPairDist_mm
        pairsValidos(k).long_chan  = channels(longIdx(iL)).canalID;
        pairsValidos(k).short_chan = channels(shortIdx(jShort)).canalID;
        pairsValidos(k).dist_mm    = dmin;
        pairsValidos(k).long_idx   = longIdx(iL);
        pairsValidos(k).short_idx  = shortIdx(jShort);
        k = k + 1;
    end
end

% ----- optional plot -----
if doPlot
    figure('Name','Pairs Short-Long (many-to-one)','NumberTitle','off'); hold on; grid on; axis equal;
    scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 60, 'g','filled'); % sources
    scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 60, 'b','filled'); % detectors
    mid = vertcat(channels.posNom_mm);
    scatter3(mid(:,1), mid(:,2), mid(:,3), 40, 'k','filled','MarkerFaceAlpha',0.5);
    colors = lines(max(1,numel(pairsValidos)));
    for p=1:numel(pairsValidos)
        Lidx = pairsValidos(p).long_idx;
        Sidx = pairsValidos(p).short_idx;
        p1 = channels(Lidx).posNom_mm; p2 = channels(Sidx).posNom_mm;
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], '-','LineWidth',1.5,'Color',colors(mod(p-1,size(colors,1))+1,:));
    end
    xlabel(sprintf('X (%s)', units)); ylabel(sprintf('Y (%s)', units)); zlabel(sprintf('Z (%s)', units));
    title('Pairs (many-to-one)');
    view(3);
end

end
