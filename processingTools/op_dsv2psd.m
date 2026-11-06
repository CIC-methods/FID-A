% op_dsv2psd.m 
% Jamie Near, Sunnybrook Research Institute, 2026
% 
% USAGE:
% PSD_outputs = op_dsv2psd(dsvFolder, startTime, endTime)
%
% DESCRIPTION:
% Read Siemens IDEA DSV files and make a publication-ready pulse-sequence 
% diagram.
%
% INPUTS:
% dsvFolder     = folder containing one ore more of the Siemens DSP-Simulator 
%                 DSV files:
%                   *_RFD.dsv  RF amplitude
%                   *_NC1.dsv  RF phase (degrees)
%                   *_GRX.dsv  X gradient
%                   *_GRY.dsv  Y gradient
%                   *_GRZ.dsv  Z gradient
%                   *_ADC.dsv  ADC activity waveform
%
%
% startTime     = start time of the plot in [seconds] (Optional, Default 0).
% endTime       = end time of the plot in [seconds] (Optional, Default end).
% invertGrads   = Boolean variable to invert the gradient polarities for 
%                 better visualization (Optional.  Default = False);
%
% OUTPUTS:
% PSD_outputs   = struct of structs variable with fields:
%                   .RFD  RF amplitude      
%                   .RFP  RF phase
%                   .NC1  NCO phase
%                   .GX   X gradient        
%                   .GY   Y gradient      
%                   .GZ  Z gradient
%                   .ADC  ADC activity waveform
%                   
%                   ->Each of the above has fields:
%                       .t      time vector [s] (column, absolute time; sample 1 is at t = 0)
%                       .amp    waveform in the DSV physical units (column, same length as .t)
%                       .units  vertical units from the DSV header
%                       .dt     sample interval [s]
%
% DSV files are read with io_readdsv.m and unpacked with op_dsv2timecourse.m.
% ------------------------------------------------------------------------

function PSD_outputs = op_dsv2psd(dsvFolder, startTime, endTime, invertGrads)

if nargin < 1 || isempty(dsvFolder),    dsvFolder   = pwd; end
if nargin < 2 || isempty(startTime),    startTime   = 0;   end
if nargin < 3 || isempty(endTime),      endTime     = Inf; end
if nargin < 4 || isempty(invertGrads),  invertGrads = False; end

% field name, file pattern, label for the plot
spec = {
    'RFD', '*RFD.dsv', '|RF|'
    'RFP', '*RFP.dsv', 'RF<'
    'NC1', '*NC1.dsv', 'NCO'
    'GX',  '*GRX.dsv', 'Gx'
    'GY',  '*GRY.dsv', 'Gy'
    'GZ',  '*GRZ.dsv', 'Gz'
    'ADC', '*ADC.dsv', 'ADC'
    };

% ------------------------------ Load ------------------------------------
PSD_outputs = struct();
for k = 1:size(spec,1)
    PSD_outputs.(spec{k,1}) = readChannel(dsvFolder, spec{k,2}, ...
        startTime, endTime, spec{k,3});
end

% Complex RF waveform, if possible
rfd = PSD_outputs.RFD;
rfp = PSD_outputs.RFP;
if ~isempty(rfd) && ~isempty(rfp) && numel(rfd.t) == numel(rfp.t) ...
        && all(abs(rfd.t - rfp.t) < 1e-12)
    PSD_outputs.RFD.complex = rfd.amp .* exp(1i*deg2rad(rfp.amp));
end

% ------------------------------ Plot ------------------------------------
plotPulseSeq(PSD_outputs, spec, invertGrads);

end


% ========================================================================
function out = readChannel(dsvFolder, pattern, startTime, endTime, label)
%READCHANNEL Load one DSV channel and trim it to [startTime, endTime].
% Returns [] if the file is absent or the window contains no samples.

out = [];

f = dir(fullfile(dsvFolder, pattern));
if isempty(f)
    fprintf('%s file (%s) not found - skipping.\n', label, pattern);
    return
elseif numel(f) > 1
    error('makePulseSeqFig:duplicateFiles', ...
        'More than one file matches %s in %s.', pattern, dsvFolder);
end

dsv  = io_readdsv(fullfile(f.folder, f.name));
wave = op_dsv2timecourse(dsv);
wave = wave(:);
N    = numel(wave);
dt   = dsv.HORIDELTA * timeUnitScale(dsv.HORIUNITNAME);   % seconds

% Index range computed ONCE; t and amp both come from the same idx, so
% their lengths always agree. Sample n is at time (n-1)*dt.
i1 = max(1, round(startTime/dt) + 1);
i2 = min(N, round(endTime/dt)   + 1);
if i2 < i1
    warning('makePulseSeqFig:emptyWindow', ...
        '%s: no samples between %g s and %g s (data span 0 to %g s).', ...
        label, startTime, endTime, (N-1)*dt);
    return
end
idx = (i1:i2).';

out.t     = (idx - 1) * dt;
out.amp   = wave(idx);
out.units = dsv.VERTUNITNAME;
out.dt    = dt;
end


% ========================================================================
function s = timeUnitScale(unitName)
%TIMEUNITSCALE Multiplier converting a DSV horizontal unit to seconds.
u = strtrim(lower(char(unitName)));
if contains(u, 'ms')
    s = 1e-3;
elseif contains(u, 'ns')
    s = 1e-9;
elseif strcmp(u, 's')
    s = 1;
elseif endsWith(u, 's')
    s = 1e-6;     % micro sign in any encoding (µs, μs, or garbled)
else
    error('makePulseSeqFig:unknownUnit', ...
        'Unrecognised DSV time unit "%s".', u);
end
end


% ========================================================================
function plotPulseSeq(PSD, spec, invertGrads)
%PLOTPULSESEQ Single-axes cartoon pulse-sequence diagram (unitless traces).

fields  = spec(:,1);
labels  = spec(:,3);
present = find(cellfun(@(f) ~isempty(PSD.(f)), fields));
if isempty(present)
    warning('makePulseSeqFig:nothingToPlot', 'No channels to plot.');
    return
end
nPlot = numel(present);

% Gx/Gy/Gz share one scale so their relative amplitudes are preserved.
gradFields = {'GX','GY','GZ'};
adcField = {'ADC'};
gradScale = 0;
for g = 1:numel(gradFields)
    if ~isempty(PSD.(gradFields{g}))
        gradScale = max(gradScale, max(abs(PSD.(gradFields{g}).amp)));
    end
end

tmin = min(cellfun(@(f) PSD.(f).t(1),   fields(present)));
tmax = max(cellfun(@(f) PSD.(f).t(end), fields(present)));
if tmax <= tmin, tmax = tmin + 1; end   % single-sample guard

fig = figure('Color','w','Name','Pulse Sequence','NumberTitle','off');
ax = axes(fig);
hold(ax,'on');

halfHeight = 0.65;                        % trace half-height (row = 1 unit)
traceColor = [0.05 0.30 0.62];

for r = 1:nPlot
    k  = present(r);
    ch = PSD.(fields{k});
    yc = nPlot - r + 1;                  % first channel at the top

    if ismember(fields{k}, gradFields)
        if invertGrads
            scale = -2*gradScale;
        else
            scale = 2*gradScale;
        end
        yc = (nPlot/2) - 0.5*r + 0.5;    
    elseif ismember(fields{k}, adcField)
        scale = 2*max(abs(ch.amp));
        yc = (nPlot/2) - 0.5*r + 0.5;
    else
        if r==1
            scale = max(abs(ch.amp));
            yc = (nPlot/2) - r + 1;
        else
            scale = 4*max(abs(ch.amp));
            yc = (nPlot/2) - 0.45*r + 0.45;
        end
    end
    if ~isfinite(scale) || scale == 0, scale = 1; end

    plot(ax, [tmin-tmin tmax-tmin], [yc yc], ':', 'Color', [0.7 0.7 0.7],'LineWidth',1.4);   % baseline
    plot(ax, ch.t-tmin, yc + halfHeight*ch.amp/scale, ...
        'LineWidth', 2, 'Color', traceColor);

    text(ax, tmin-tmin - 0.01*(tmax-tmin), yc, labels{k}, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontName','Arial', 'FontSize',20, 'Interpreter','none', ...
        'Clipping','off');
end

xlim(ax, [tmin-tmin tmax-tmin]);
ylim(ax, [0.25 (nPlot/2)+0.7]);
ax.YAxis.Visible = 'off';
box(ax,'off');
set(ax, 'FontName','Arial', 'FontSize',14, 'Layer','top', ...
    'Position',[0.10 0.14 0.87 0.80],'LineWidth',1);
xlabel(ax, 'Time (s)', 'FontSize',16);
title(ax, 'Pulse Sequence', 'FontWeight','bold','FontSize',20);
end
