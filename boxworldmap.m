function h = boxworldmap(varargin) %ltlim, lnlim, Opt.latgrid, Opt.longrid)
%BOXWORLDMAP Worldmap axis expanded to box
% 
% h = boxworldmap;
% h = boxworldmap(ltlim, lnlim)
% h = boxworldmap(region)
% h = boxworldmap(Z, R)
% h = boxworldmap(..., Name, Value)
%
% This function is a wrapper around worldmap (Mapping Toolbox).  It
% modifies worldmap results by turning on visibility of the underlying axis
% and extending grid lines and grid labels to fill that entire box rather
% than just the map axis frame.  This allows one to more easily plot beyond
% the boundaries of the map axis frame (although data must be manually
% projected first).
%
% Note: This function is still a WIP... haven't found a robust way to label
% the grid lines that works for all regions without label overlap.
%
% Input variables:
%
%   See worldmap (Mapping Toolbox)
%
% Optional input variables (passed as parameter/value pairs
%
%   latgrid:    vector of latitude values where gridlines should be placed.
%               Use [] for none, NaN to inherit from the default map axis.
%               [NaN]
%
%   longrid:    vector of longitude values where gridlines should be placed.
%               Use [] for none, NaN to inherit from the default map axis.
%               [NaN]
%
% Output variables:
%
%   h:          structure of graphics objects:
%
%               ax: axis handle
%
%               gpar:   handles to latitude (parallels) grid line objects
%
%               gmer:   handles to longitude (meridian) grid line objects
%
%               lblpar: handles to parallel label text objects
%
%               lblmer: handles to meridian label text objects

% Copyright 2023-2024 Kelly Kearney


% Parse input

p = inputParser;

p.addOptional('wm1', NaN, @(x) validateattributes(x, {'string', 'char', 'numeric'}, {}));
p.addOptional('wm2', NaN, @(x) validateattributes(x, {'numeric'}, {}));
p.addParameter('latgrid', NaN, @(x) validateattributes(x, {'numeric'}, {'vector', '>=', -90, '<=', 90}));
p.addParameter('longrid', NaN, @(x) validateattributes(x, {'numeric'}, {'vector', '>=', -180, '<=', 360}));
p.parse(varargin{:});
Opt = p.Results;

% Call worldmap

h.ax = gca;

if isnan(Opt.wm2)
    if isnan(Opt.wm1)
        worldmap;
    else
        worldmap(Opt.wm1);
    end
else
    worldmap(Opt.wm1, Opt.wm2);
end

mstruct = getm(h.ax);

% Get grid info if needed

if isnan(Opt.latgrid)
    if isscalar(mstruct.plinelocation)
        Opt.latgrid = sort([0:mstruct.plinelocation:90 -(0:mstruct.plinelocation:90)]);
    else
        Opt.latgrid = mstruct.plinelocation;
    end
end
if isnan(Opt.longrid)
    if isscalar(mstruct.mlinelocation)
        Opt.longrid = sort([0:mstruct.mlinelocation:180 -(0:mstruct.plinelocation:180)]);
    else
        Opt.longrid = mstruct.plinelocation;
    end
end

set(h.ax, 'visible', 'on');
setm(h.ax, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off');

% Grid lines

xbox = [h.ax.XLim mean(h.ax.XLim)];
ybox = [h.ax.YLim mean(h.ax.YLim)];
[xbox,ybox] = ndgrid(xbox, ybox);
[ltbox, lnbox] = projinv(mstruct, xbox(:), ybox(:));

ltlim2 = [min(ltbox) max(ltbox)];
lnlim2 = [min(wrapTo360(lnbox)) max(wrapTo360(lnbox))];

[xpar, ypar] = deal(cell(size(Opt.latgrid)));
for ip = 1:length(Opt.latgrid)
    [ltpar, lnpar] = interpm(Opt.latgrid(ip)*ones(1,2), lnlim2, 0.5);
    [xpar{ip}, ypar{ip}] = projfwd(mstruct, ltpar, lnpar);
end
[xmer, ymer] = deal(cell(size(Opt.longrid)));
for im = 1:length(Opt.longrid)
    [ltmer, lnmer] = interpm(ltlim2, Opt.longrid(im)*ones(1,2), 0.5);
    [xmer{im}, ymer{im}] = projfwd(mstruct, ltmer, lnmer);
end

% h.gpar = cellfun(@(x,y) patch([x;NaN],[y;NaN], 'r'), xpar, ypar);
% h.gmer = cellfun(@(x,y) patch([x;NaN],[y;NaN], 'r'), xmer, ymer);

h.gpar = cellfun(@(x,y) plot(x,y, 'r'), xpar, ypar);
h.gmer = cellfun(@(x,y) plot(x,y, 'r'), xmer, ymer);

set([h.gpar h.gmer], 'color', [0.84706 0.86275 0.83922], 'linewidth', mstruct.glinewidth, 'linestyle', mstruct.glinestyle);

% Grid line labels

pbox = polyshape(h.ax.XLim([1 1 2 2 1]), h.ax.YLim([1 2 2 1 1]));

[xplbl, yplbl] = deal(nan(size(Opt.latgrid)));
for ii = 1:length(Opt.latgrid)
    [in,out] = intersect(pbox, [xpar{ii} ypar{ii}]);
    if  ~isempty(in)
        xplbl(ii) = in(1,1); % first in is left
        yplbl(ii) = in(1,2);
    end
end
[xmlbl, ymlbl] = deal(nan(size(Opt.longrid)));
for ii = 1:length(Opt.longrid)
    [in,out] = intersect(pbox, [xmer{ii} ymer{ii}]);
    if ~isempty(in)
        xmlbl(ii) = in(1,1); % last in is top
        ymlbl(ii) = in(1,2); % last in is top
    end
end


h.lblpar = text(xplbl+diff(h.ax.XLim)*0.005, yplbl, strrep(cellstr(angl2str(Opt.latgrid,'ns')),'^',''), 'horiz', 'left');
h.lblmer = text(xmlbl, ymlbl, strrep(cellstr(angl2str(wrapTo180(Opt.longrid),'ew')),'^',''), 'horiz', 'center', 'vert', 'top');



