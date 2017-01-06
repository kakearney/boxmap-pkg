function h = maskmap(ax)
%MASKMAP Mask a map axis to make it look like a normal cartesian axis
%
% h = maskmap(ax)
%
% This function provides a cosmetic overlay to a map axis, constraining the
% plotted region to the largest rectangle available within the map axis
% frame.  The frame and meridian/parallel labels are removed from the map
% axis, and two new cartesian axes are added that include appropriate x/y
% ticks where the map grid lines intersect the new axis rectangle.
%
% Input variables:
%
%   ax: handle to a map axis
%
% Output variables:
%
%   h:  structure with handles to new graphics objects
%
%       mask:   a patch object, covering up all objects in ax that fall
%               outside the maximum enclosed rectange of the map canvas
%
%       ax:     1 x 2 array of axis objects.  ax(1) holds the left and
%               bottom axis, and ax(2) holds the right and top axis.

% Copyright 2016 Kelly Kearney


frm = findall(ax, 'tag', 'Frame');
mer = findall(ax, 'tag', 'Meridian');
par = findall(ax, 'tag', 'Parallel');
lbl = findall(ax, 'tag', 'PLabel', '-or', 'tag', 'MLabel');
fig = ancestor(ax, 'figure');

mstruct = getm(ax);

xmask = linspace(ax.XLim(1), ax.XLim(2), 100);
ymask = linspace(ax.YLim(1), ax.YLim(2), 100);
[xg, yg] = meshgrid(xmask, ymask);
isax = inpolygon(xg, yg, frm.XData, frm.YData);

[cc,hh,ww,mm] = FindLargestRectangles(isax);
[~, pos] = max(cc(:)); % max by area

b = bwboundaries(mm);
bx = xmask(b{1}(:,2));
by = ymask(b{1}(:,1));

bxlim = minmax(bx);
bylim = minmax(by);

bx = bxlim([1 1 2 2 1]);
by = bylim([1 2 2 1 1]);

xax = ax.XLim([1 1 2 2 1]);
yax = ax.YLim([1 2 2 1 1]);
[xpoly, ypoly] = polybool('-', xax, yax, bx, by);
[f,v] = poly2fv(xpoly, ypoly);

h.mask = patch('faces', f, 'vertices', v, 'facecolor', fig.Color, ...
               'edgecolor', 'none', 'parent', ax);
% set(h.mask, 'FaceAlpha', 0.9);

[xmer, ymer] = polyxpoly(bx, by, mer.XData, mer.YData);
[xpar, ypar] = polyxpoly(bx, by, par.XData, par.YData);

tol = 1e-6;
ist = abs(bylim(2) - ymer) < tol;
isb = abs(bylim(1) - ymer) < tol;
isl = abs(bxlim(1) - xpar) < tol;
isr = abs(bxlim(2) - xpar) < tol;

xtk1 = xmer(ist);
[~, xtklbl1] = minvtran(mstruct, xmer(ist), ymer(ist));
xtk2 = xmer(isb);
[~, xtklbl2] = minvtran(mstruct, xmer(isb), ymer(isb));

ytk1 = ypar(isl);
[ytklbl1, ~] = minvtran(mstruct, xpar(isl), ypar(isl));
ytk2 = ypar(isr);
[ytklbl2, ~] = minvtran(mstruct, xpar(isr), ypar(isr));

xtkt = sortrows([xtk1 round(xtklbl1,2)]);
xtkb = sortrows([xtk2 round(xtklbl2,2)]);

ytkl = sortrows([ytk1 round(ytklbl1,2)]);
ytkr = sortrows([ytk2 round(ytklbl2,2)]);

unit = fig.Units;
fig.Units = 'normalized';
[bxfig, byfig] = axescoord2figurecoord(bxlim, bylim, ax);
fig.Units = unit;

h.ax(1) = axes('position', [bxfig(1) byfig(1) diff(bxfig) diff(byfig)], ...
                  'box', 'off', 'yaxisloc', 'left', 'xaxisloc', 'bottom', ...
                  'xtick', xtkb(:,1), 'ytick', ytkl(:,1), ...
                  'xticklabel', cellstr(num2str(xtkb(:,2))), ...
                  'yticklabel', cellstr(num2str(ytkl(:,2))), ...
                  'xlim', bxlim, 'ylim', bylim, 'color', 'none');
h.ax(2) = axes('position', [bxfig(1) byfig(1) diff(bxfig) diff(byfig)], ...
                  'box', 'off', 'yaxisloc', 'right', 'xaxisloc', 'top', ...
                  'xtick', xtkt(:,1), 'ytick', ytkr(:,1), ...
                  'xticklabel', cellstr(num2str(xtkt(:,2))), ...
                  'yticklabel', cellstr(num2str(ytkr(:,2))), ...
                  'xlim', bxlim, 'ylim', bylim, 'color', 'none');

delete(frm);
delete(lbl);