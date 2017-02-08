function h = globezoommap(ax, origin, zoomlev, roundlev)
%GLOBEZOOMMAP Plot zoomed-in orthographic projection map
%
% h = globezoommap(ax, origin, zoomlev, roundlev)
%
% This function provides a cosmetic overlay to a map axis, constraining the
% plotted region to the largest rectangle available within the map axis
% frame.  The frame and meridian/parallel labels are removed from the map
% axis, and two new cartesian axes are added that include appropriate x/y
% ticks where the map grid lines intersect the new axis rectangle.
%
% Input variables:
%
%   ax: handle to axis
%
% Output variables:
%
%   h:  structure with handles to new graphics objects
%
%       mask:   a patch object, covering up all objects in ax that fall
%               outside the maximum enclosed rectange of the map canvas
%
%       ax:     1 x 2 array of axis objects.  ax(1) holds the left and
%               bottom axis, and ax(2) holds the right and top axis, ax(3)
%               holds the axis-masking axis

% Copyright 2016 Kelly Kearney

axes(ax);
axesm('ortho', 'grid', 'on', 'origin', origin);
setm(ax, 'MeridianLabel', 'off', 'ParallelLabel', 'off', 'frame', 'off', 'mlinelocation', 15, 'gcolor', rgb('light gray'));

frm = findall(ax, 'tag', 'Frame');
mer = findall(ax, 'tag', 'Meridian');
par = findall(ax, 'tag', 'Parallel');
lbl = findall(ax, 'tag', 'PLabel', '-or', 'tag', 'MLabel');
fig = ancestor(ax, 'figure');

zoom(zoomlev);

% Where is the axes box in coordinate space?

% xedge = linspace(ax.XLim(1), ax.XLim(2), 21);
% yedge = linspace(ax.YLim(1), ax.YLim(2), 21);
% xbox = [xedge ones(size(yedge))*xedge(end) xedge(end:-1:1) ones(size(yedge))*xedge(1)];
% ybox = [yedge(end)*ones(size(xedge)) yedge(end:-1:1) yedge(1)*ones(size(xedge)) yedge];
% [ltbox, lnbox] = minvtran(xbox, ybox);

bx = ax.XLim([1 2 2 1 1]);
by = ax.YLim([2 2 1 1 2]);

% h.axnew = axes('Position', ax.Position, 'xlim', 

mstruct = getm(ax);

% Masking polygon

pos = plotboxpos(ax);

xaxlim = [pos(1) pos(1)+pos(3)];
yaxlim = [pos(2) pos(2)+pos(4)];
[xpoly, ypoly] = polybool('-', [0 0 1 1 0], [0 1 1 0 0], xaxlim([1 1 2 2 1]), yaxlim([1 2 2 1 1]));
[f,v] = poly2fv(xpoly, ypoly);

% Where will new axis need ticks

[xmer, ymer] = polyxpoly(bx, by, mer.XData, mer.YData);
[xpar, ypar] = polyxpoly(bx, by, par.XData, par.YData);

tol = 1e-6;
ist = abs(ax.YLim(2) - ymer) < tol;
isb = abs(ax.YLim(1) - ymer) < tol;
isl = abs(ax.XLim(1) - xpar) < tol;
isr = abs(ax.XLim(2) - xpar) < tol;

xtk1 = xmer(ist);
[~, xtklbl1] = minvtran(mstruct, xmer(ist), ymer(ist));
xtk2 = xmer(isb);
[~, xtklbl2] = minvtran(mstruct, xmer(isb), ymer(isb));

ytk1 = ypar(isl);
[ytklbl1, ~] = minvtran(mstruct, xpar(isl), ypar(isl));
ytk2 = ypar(isr);
[ytklbl2, ~] = minvtran(mstruct, xpar(isr), ypar(isr));

xtkt = sortrows(unique([xtk1 round(xtklbl1,roundlev)], 'rows'));
xtkb = sortrows(unique([xtk2 round(xtklbl2,roundlev)], 'rows'));

ytkl = sortrows(unique([ytk1 round(ytklbl1,roundlev)], 'rows'));
ytkr = sortrows(unique([ytk2 round(ytklbl2,roundlev)], 'rows'));

% text(xmer(ist), ymer(ist), strtrim(cellstr(num2str(xtklbl1))));
% text(xmer(isb), ymer(isb), strtrim(cellstr(num2str(xtklbl2))));
% text(xpar(isl), ypar(isl), strtrim(cellstr(num2str(ytklbl1))));
% text(xpar(isr), ypar(isr), strtrim(cellstr(num2str(ytklbl2))));
% plot(xmer, ymer, 'r.');
% plot(xpar, ypar, 'b.');

% unit = fig.Units;
% fig.Units = 'normalized';
% [bxfig, byfig] = axescoord2figurecoord(bxlim, bylim, ax);
% fig.Units = unit;

h.ax(1) = axes('position', pos, 'xlim', ax.XLim, 'ylim', ax.YLim, ...
               'box', 'off', 'yaxisloc', 'left', 'xaxisloc', 'bottom', ...
               'xtick', xtkb(:,1), 'ytick', ytkl(:,1), ...
               'xticklabel', cellstr(num2str(xtkb(:,2))), ...
               'yticklabel', cellstr(num2str(ytkl(:,2))), ...
               'color', 'none');
           
h.ax(2) = axes('position', pos, 'xlim', ax.XLim, 'ylim', ax.YLim, ...
               'box', 'off', 'yaxisloc', 'right', 'xaxisloc', 'top', ...
               'xtick', xtkt(:,1), 'ytick', ytkr(:,1), ...
               'xticklabel', cellstr(num2str(xtkt(:,2))), ...
               'yticklabel', cellstr(num2str(ytkr(:,2))), ...
               'color', 'none');
h.ax(3) = axes('position', [0 0 1 1], 'color', 'none', 'visible', 'off', 'xlim', [0 1], 'ylim', [0 1]);
h.mask = patch('faces', f, 'vertices', v, 'facecolor', fig.Color, ...
               'edgecolor', 'none', 'parent', h.ax(3));
set(h.mask, 'FaceAlpha', 0.9);

uistack(h.ax(1:2), 'top');


           
           
%            blah
           
% h.ax(1) = axes('position', [bxfig(1) byfig(1) diff(bxfig) diff(byfig)], ...
%                   'box', 'off', 'yaxisloc', 'left', 'xaxisloc', 'bottom', ...
%                   'xtick', xtkb(:,1), 'ytick', ytkl(:,1), ...
%                   'xticklabel', cellstr(num2str(xtkb(:,2))), ...
%                   'yticklabel', cellstr(num2str(ytkl(:,2))), ...
%                   'xlim', bxlim, 'ylim', bylim, 'color', 'none');
% h.ax(2) = axes('position', [bxfig(1) byfig(1) diff(bxfig) diff(byfig)], ...
%                   'box', 'off', 'yaxisloc', 'right', 'xaxisloc', 'top', ...
%                   'xtick', xtkt(:,1), 'ytick', ytkr(:,1), ...
%                   'xticklabel', cellstr(num2str(xtkt(:,2))), ...
%                   'yticklabel', cellstr(num2str(ytkr(:,2))), ...
%                   'xlim', bxlim, 'ylim', bylim, 'color', 'none');
% 
% delete(frm);
% delete(lbl);