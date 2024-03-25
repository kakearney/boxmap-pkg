function [m, Lim, hax] = boxmap(latlim, lonlim, varargin)
%BOXMAP Plot map to a regtangular axis
%
% [m, Lim, hax] = boxmap(latlim, lonlim, p1, v1, ...)
%
% This function plots sets up an axis to plot projected map data, but to
% maintain a rectangular plotting region.  It is primarily intended for
% conical projections.
%
% Input variables:
%
%   latlim:     desired latitude limits
%
%   lonlim:     desired longitude limits
%
% Optional input variables (passed as parameter/value pairs):
%
%   proj:       projection ['eqdconic']
%
%   ax:         handle of axis to modify [gca]
%
%   lontick:    tick interval for longitude ticks [1]
%
%   lattick:    tick interval for latitude ticks [1]
%
%   format:     function to format tick labels, takes as input a scalar
%               number and returns a string [@(x) num2str(x, '%.2f')]
%
%   plot:       logical scalar, true to plot and false to simply return m
%               and Lim outputs
%
%   dx:         minimum distance between points, in km, for drawing outline
%               around region of interest [100]
%
%   npt:        number of points used to draw box (increase if corners are
%               being cut off in your x/y limits [200]
%
% Output variables:
%
%   hax:        1 x 2 array of axis handles.  The first axis is the
%               same as the 'ax' input parameter, and controls the left and
%               bottom axes for tick purposes; the second is layered on top
%               and controls the right and top axes. 
%
%   m:          map projection struct used to project lat/lon data onto the
%               resulting axis coordinate systems
%
%   Lim:        1 x 1 structure
%
%               lat:    latitude limits for the resulting plotted space
%
%               lon:    longitude limits for the resulting plotted space
%
%               x:      x limits of the plot
%
%               y:      y limits of the plot

% Copyright 2015 Kelly Kearney

% Parse input

p = inputParser;
p.addParameter('proj', 'eqdconic');
p.addParameter('plot', true);
p.addParameter('ax', []);
p.addParameter('lontick', 1);
p.addParameter('lattick', 1);
p.addParameter('format', @(x) num2str(x, '%.2f'));
p.addParameter('npt', 200);
p.addParameter('dx', 100);
p.parse(varargin{:});

Opt = p.Results;

if Opt.plot && isempty(Opt.ax)
    Opt.ax = gca;
end

% Map projection structure

m = defaultm(Opt.proj);
m.maplatlimit = latlim;
m.maplonlimit = lonlim;
m = defaultm(m);

% Calculate box limits

minmax = @(x) [min(x(:)) max(x(:))];

[tklat, tklon] = interpm(latlim([1 2 2 1 1]), lonlim([1 1 2 2 1]), km2deg(Opt.dx), 'rh');
[xlim, ylim] = mfwdtran(m, tklat, tklon);
xlim = minmax(xlim);
ylim = minmax(ylim);
xybox = interparc(linspace(0,1,Opt.npt), xlim([1 1 2 2 1]), ylim([1 2 2 1 1]), 'linear');

[ltbox, lnbox] = minvtran(m, xybox(:,1), xybox(:,2));
lnbox = wrapTo360(lnbox);

latlim2 = minmax(ltbox);
lonlim2 = minmax(lnbox);

Lim.lat = latlim2;
Lim.lon = lonlim2;
Lim.x = xlim;
Lim.y = ylim;

if Opt.plot

    % Tick locations

    rlo = @(x,tol) floor(x./tol).*tol;
    rhi = @(x,tol) ceil(x./tol).*tol;

    ist = diff(xybox(:,1)) > 0;
    isb = diff(xybox(:,1)) < 0;
    isl = diff(xybox(:,2)) > 0;
    isr = diff(xybox(:,2)) < 0;

    lntk = rlo(lonlim2(1),Opt.lontick):Opt.lontick:rhi(lonlim2(2),Opt.lontick);
    lttk = rlo(latlim2(1),Opt.lattick):Opt.lattick:rhi(latlim2(2),Opt.lattick);

    tickt = interp1(lnbox(ist), xybox(ist,1), lntk);
    tickb = interp1(lnbox(isb), xybox(isb,1), lntk); 
    tickl = interp1(ltbox(isl), xybox(isl,2), lttk);
    tickr = interp1(ltbox(isr), xybox(isr,2), lttk); 

    lnlbl = arrayfun(Opt.format, lntk, 'uni', 0);
    ltlbl = arrayfun(Opt.format, lttk, 'uni', 0);

    tick = {tickl tickt tickr tickb};
    lbl = {ltlbl lnlbl ltlbl lnlbl};

    for ii = 1:length(tick)
        isn = isnan(tick{ii});
        tick{ii} = tick{ii}(~isn);
        lbl{ii}   = lbl{ii}(~isn);
    end

    % Adjust plot

    axprop = {'xlim', xlim, ...
              'ylim', ylim, ...
              'dataaspectratio', [1 1 1], ...
              'box', 'off'};
    lbax = Opt.ax;
    set(lbax, axprop{:}, ...
        'yaxislocation', 'left', ...
        'xaxislocation', 'bottom', ...
        'xtick', tick{4}, ...
        'xticklabel', lbl{4}, ...
        'ytick', tick{1}, ...
        'yticklabel', lbl{1});
    rtax = axes('position', get(lbax, 'position'), axprop{:}, ...
        'yaxislocation', 'right', ...
        'xaxislocation', 'top', ...
        'xtick', tick{2}, ...
        'xticklabel', lbl{2}, ...
        'ytick', tick{3}, ...
        'yticklabel', lbl{3}, ...
        'color', 'none');

    hax = [lbax rtax];
else
    hax = [];
end








