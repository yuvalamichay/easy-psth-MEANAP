function grammObject = RasterPSTH(spikeTimes, events, varargin)
%   easy.RasterPSTH(spikeTimes,events) plots a raster and psth for the
%   spike times aligned to an event. spikeTimes should be a column vector of
%   S x 1 spike time events. events should be a S x 2 cell array, where S 
%   is the number of alignments required. First column should contain a string label
%   for that alignment. Second column should contain a column vector of
%   times for that event.
%   The function returns the grammObject (graphical object for the
%   figure). This function requires gramm and spikes toolboxes to be
%   installed.

default_window = [-1 +1];
default_psthBinWidth = 1/1000;
default_psthSmoothWidth = 50/1000;

%Validate inputs
p = inputParser;
addRequired(p,'spikeTimes',@(x) iscolumn(x) & ~isempty(x));
addRequired(p,'events',@iscell);
addParameter(p,'splitBy',[],@iscell);
addParameter(p,'sortBy',[],@iscell);
addParameter(p,'psthBinWidth',default_psthBinWidth,@isscalar);
addParameter(p,'psthSmoothWidth',default_psthSmoothWidth,@isscalar);
addParameter(p,'window',default_window,@isrow);
addParameter(p,'ylim',[],@(x) length(x)==2 & diff(x)>0);
addParameter(p,'titleText','');
addParameter(p,'splitByColours',[],@iscell);
addParameter(p,'saveFigure',[],@ischar);
parse(p,spikeTimes,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('gramm')),'Please install the gramm toolbox https://github.com/piermorel/gramm');
assert(~isempty(which('psthAndBA')),'Please install the spikes toolbox https://github.com/cortex-lab/spikes/');

spikeTimes = p.Results.spikeTimes;
events = p.Results.events;
num_alignments = size(p.Results.events,1);

if isempty(p.Results.splitBy)
    splitBy = cell(num_alignments,2);
else
    splitBy = p.Results.splitBy;
end

if isempty(p.Results.sortBy)
    sortBy = cell(num_alignments,2);
else
    sortBy = p.Results.sortBy;
end

%For each event type, compute the aligned spike times and binned psths
bins                    = cell(num_alignments,1);
binnedArray_smoothed    = cell(num_alignments,1);
spikeTimes_byEvent      = cell(num_alignments,1);
for i = 1:num_alignments
    [bins{i}, binnedArray_smoothed{i}, spikeTimes_byEvent{i}] = alignEventTimes(spikeTimes, events{i,2}, p.Results.window, p.Results.psthBinWidth, p.Results.psthSmoothWidth);
    numBinsToTrim = round((p.Results.psthSmoothWidth*2)/p.Results.psthBinWidth);
    binnedArray_smoothed{i}(:,1:numBinsToTrim) = NaN;

    if ~isempty(sortBy{i,2})
        sortBy{i,2} = sortBy{i,2} - events{i,2}; %Calculate sort time relative to event time
        [sortBy{i,2},sortIdx] = sort(sortBy{i,2},'ascend');
        spikeTimes_byEvent{i} = spikeTimes_byEvent{i}(sortIdx);
        binnedArray_smoothed{i} = binnedArray_smoothed{i}(sortIdx,:);
        if ~isempty(splitBy{i,2})
            splitBy{i,2} = splitBy{i,2}(sortIdx);
        end
    end
end

clear grammObject;

%For each event type, create gramm visual object
for i = 1:num_alignments
    title = events(i,1);
    if ~isempty(splitBy{i,2})
        title = [title,[' split by ' splitBy{i,1}]];
    end
    if ~isempty(sortBy{i,2})
        title = [title,[' sorted by ' sortBy{i,1}]];
    end

    %create gramm raster object
    grammObject(1,i) = gramm('x',spikeTimes_byEvent{i},'color',splitBy{i,2});
    grammObject(1,i).geom_raster('geom','point');
    grammObject(1,i).geom_vline('xintercept', 0);
    grammObject(1,i).set_names('y','Event number','x',events{i,1},'color','');
    grammObject(1,i).set_title(title);
    grammObject(1,i).set_point_options('base_size',2);
    grammObject(1,i).axe_property('XLim',p.Results.window,'Ydir','reverse');
    grammObject(1,i).axe_property('Ytick','','Ycolor','none');

    %create gramm psth object
    grammObject(2,i) = gramm('x',bins{i},'y',binnedArray_smoothed{i},'color',splitBy{i,2});
    grammObject(2,i).stat_summary('setylim', true,'type', 'sem');
    grammObject(2,i).geom_vline('xintercept', 0);
    grammObject(2,i).set_names('y','Spikes/sec','x',events{i,1},'color','');
    grammObject(2,i).axe_property('XLim',p.Results.window);

    if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{i})
        numSplits = length(unique(splitBy{i,2}));
        numColours = size(p.Results.splitByColours{i},1);
        assert(numSplits == numColours, sprintf('%s: %d colours specified but %d splitting conditions',splitBy{i,1},numColours,numSplits));
        grammObject(:,i).set_color_options('map',p.Results.splitByColours{i},'n_color',numColours,'n_lightness',1);
    end
end

%Render figure;
grammObject.set_title(p.Results.titleText);
grammObject.axe_property('TickDir','out');
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
grammObject.draw();

% Add sorting icons (if any)
for i = 1:num_alignments
    if ~isempty(sortBy{i,2})
        grammObject(1,i).update('x',num2cell(sortBy{i,2}));
        grammObject(1,i).geom_raster('geom','point');
        grammObject(1,i).no_legend();
    end
end

if any(~cellfun(@isempty,sortBy(:,2)))
    grammObject.draw();
end

%Change icons to black triangles
grammObjectRasters = grammObject(1,:);
hasSort = ~cellfun(@isempty,sortBy(:,2));
if any(hasSort)
    hx = [grammObjectRasters(1,hasSort).results];
    for i = 1:length(hx)
        set([hx(i).geom_raster_handle],'Marker','v','MarkerSize',4,'MarkerFaceColor',[0 0 0]);
    end
end

%Match all PSTH firing rate axes
if num_alignments > 1
    linkaxes([grammObject(2,:).facet_axes_handles],'y');
end

%If ylim provided, set all YLims
if ~isempty(p.Results.ylim)
    set([grammObject(2,:).facet_axes_handles],'ylim',p.Results.ylim);
end

%add small text indicating psth parameters
text = sprintf('%d ms binning\n%d ms smoothing',1000*p.Results.psthBinWidth, 1000*p.Results.psthSmoothWidth);
axPos = grammObject(2,1).facet_axes_handles.Position; %first psth plot handle position
annotation('textbox', axPos, 'string',text,'HorizontalAlignment','left','EdgeColor','none');

%Save figure if desired
if ~isempty(p.Results.saveFigure)
    [fpath,fname,fext] = fileparts(p.Results.saveFigure);
    grammObject.export('file_name',fname,'export_path',fpath,'file_type',fext(2:end));
end
end

function [bins, binnedArray_smoothed, spikeTimes_byEvent] = alignEventTimes(spikeTimes, eventTimes, window, psthBinWidth, psthSmoothWidth)
%This function computes two things: 1) spike times in a window surrounding
%an event time. 2) binned spike counts around the window (useful for PSTH
%later).
[~, bins, ~, ~, ~, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinWidth);
binnedArray(isnan(eventTimes),:) = NaN;
smoothFilt = myGaussWin(psthSmoothWidth, 1/psthBinWidth);
smoothFilt(1:round(numel(smoothFilt)/2)-1) = 0;
smoothFilt = smoothFilt./sum(smoothFilt);
binnedArray_smoothed = conv2(smoothFilt,1,binnedArray', 'same')';
binnedArray_smoothed = binnedArray_smoothed./psthBinWidth;
out = WithinRanges(spikeTimes, eventTimes + window,(1:length(eventTimes))','matrix');
spikeTimes_byEvent = arrayfun( @(n) spikeTimes(logical(out(:,n))) - eventTimes(n) , 1:length(eventTimes), 'uni', 0)';
end
