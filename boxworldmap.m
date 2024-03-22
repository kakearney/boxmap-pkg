function h = boxworldmap(ltlim, lnlim, glt, gln)
%BOXWORLDMAP Worldmap axis expanded to box
% 
% h = boxworldmap(ltlim, lnlim, glt, gln)
%
% 


% ltlim = [46 70];
% lnlim = [160 205];
% glt = 50:10:70;
% gln = 145:15:220;


% h = plotgrid('size', [1 1], 'mar', 0.01, 'ml', 0.03, 'mb', 0.03);
% setpos(h.fig, '# # 24cm 20cm');
% set(h.fig, 'color', 'w');

h.ax = gca;

worldmap(ltlim, lnlim);
set(h.ax, 'visible', 'on');
setm(h.ax, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off');

% Borders

% [bx,by] = mfwdtran(bp.Vertices(:,2), bp.Vertices(:,1));
% h.cst = plot(bx,by, 'color', rgb('gray'));

% tmp = plotromsrho(Grd, Grd.h, false);

% [px, py] = mfwdtran(Grd.lat_psi, Grd.lon_psi);
% tmp = pcolorpad(px, py, Grd.h(2:end-1,2:end-1));
% shading flat;
% uistack(tmp, 'bottom');

% Grid lines

xbox = [h.ax.XLim mean(h.ax.XLim)];
ybox = [h.ax.YLim mean(h.ax.YLim)];
[xbox,ybox] = ndgrid(xbox, ybox);
[ltbox, lnbox] = minvtran(xbox(:), ybox(:));

ltlim2 = [min(ltbox) max(ltbox)];
lnlim2 = [min(wrapTo360(lnbox)) max(wrapTo360(lnbox))];

[xpar, ypar] = deal(cell(size(glt)));
for ip = 1:length(glt)
    [ltpar, lnpar] = interpm(glt(ip)*ones(1,2), lnlim2, 0.5);
    [xpar{ip}, ypar{ip}] = mfwdtran(ltpar, lnpar);
end
[xmer, ymer] = deal(cell(size(gln)));
for im = 1:length(gln)
    [ltmer, lnmer] = interpm(ltlim2, gln(im)*ones(1,2), 0.5);
    [xmer{im}, ymer{im}] = mfwdtran(ltmer, lnmer);
end

% h.gpar = cellfun(@(x,y) patch([x;NaN],[y;NaN], 'r'), xpar, ypar);
% h.gmer = cellfun(@(x,y) patch([x;NaN],[y;NaN], 'r'), xmer, ymer);

h.gpar = cellfun(@(x,y) plot(x,y, 'r'), xpar, ypar);
h.gmer = cellfun(@(x,y) plot(x,y, 'r'), xmer, ymer);

set([h.gpar h.gmer], 'color', [0.84706 0.86275 0.83922], 'linewidth', 0.5, 'linestyle', '-');

% Grid line labels

pbox = polyshape(h.ax.XLim([1 1 2 2 1]), h.ax.YLim([1 2 2 1 1]));

[xplbl, yplbl] = deal(nan(size(glt)));
for ii = 1:length(glt)
    [in,out] = intersect(pbox, [xpar{ii} ypar{ii}]);
    if  ~isempty(in)
        xplbl(ii) = in(1,1); % first in is left
        yplbl(ii) = in(1,2);
    end
end
[xmlbl, ymlbl] = deal(nan(size(gln)));
for ii = 1:length(gln)
    [in,out] = intersect(pbox, [xmer{ii} ymer{ii}]);
    if ~isempty(in)
%         xmlbl(ii) = in(end,1); % last in is top
%         ymlbl(ii) = in(end,2); % last in is top
        xmlbl(ii) = in(1,1); % last in is top
        ymlbl(ii) = in(1,2); % last in is top
    end
%     lblpt = [lblpt; intersect(in, out, 'rows')];
end



h.lblpar = text(xplbl+diff(h.ax.XLim)*0.005, yplbl, strrep(cellstr(angl2str(glt,'ns')),'^',''), 'horiz', 'left');
h.lblmer = text(xmlbl, ymlbl, strrep(cellstr(angl2str(wrapTo180(gln),'ew')),'^',''), 'horiz', 'center', 'vert', 'top');



