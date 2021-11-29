% op_CSIplot.m
% Brenden Kadota, McGill University 2019.
% 
% USAGE:
% sim_plotCSI(in)
% sim_plotCSI(in, 'min_range', 0, 'max_range', 4, 'coilNum', 3, 'averageNum', 1, 
%             'x_indecies', [10,20], 'y_range', [-100, 100])
% 
% DESCRIPTION:
% This function takes a processed MRSI twix file and plots the data using 
% the matlab plot function. Horizontal coordinates are the x coordinates
% and vertical are y coordinates. The plots at specific (x,y) positions are
% the spectral coordintaes
%
% INPUTS:
% in          = input cell array of simulated spectra from a spatially resolved simulation
%   NAME VALUE ARGUMENTS (ALL OPTIONAL)
%   Name: 'min_range' Value: double, minimum ppm range or time range to plot
%   Name: 'max_range' Value: double, maximum ppm range or time range to plot
%   Name: 'coilNum' Value: integer, index of coil to plot
%   Name: 'averageNum' Value: integer, index of average to plot
%   Name: 'x_indecies' Value: (2,1) integer, bounds (inclusive) of x indecies to plot.
%                              First index must be smaller than the second.
%   Name: 'y_indecies' Value: (2,1) integer, bounds (inclusive) of y indecies to plot.
%                              First index must be smaller than the second.
%   Name: 'x_range' Value: (2,1) double, bounds (inclusive) of x values to plot.
%                              First index must be smaller than the second.
%   Name: 'y_range' Value: (2,1) double, bounds (inclusive) of y indecies to plot.
%                              First index must be smaller than the second.
%                                           
%                
%
% OUTPUTS: 
% fig   = figure handle.

function fig = op_CSIPlot(in, plot_type, ppm_arguments, indecies)
arguments
    in (1,1) struct
    plot_type.plane_type char {mustBeMemberi(plot_type.plane_type, {'real', 'imag', 'abs'})} = 'real'
    plot_type.plot_dim char {mustBeMemberi(plot_type.plot_dim, {'ppm', 't', 'time'})} = 't'
    ppm_arguments.min_range (1,1) double = get_default_range(in, 'min')
    ppm_arguments.max_range (1,1) double = get_default_range(in, 'max')
    indecies.coilNum (1,1) double = 1
    indecies.averageNum (1,1) double = 1
    indecies.x_indecies (2,1) double
    indecies.y_indecies (2,1) double
    indecies.x_range (2,1) double
    indecies.y_range (2,1) double
end

%Argument checks
if(~isfield(in, 'specs'))
    error('please fourier transform along the spatial dimension before plotting')
end
if(in.dims.coils == 0)
    if(indecies.coilNum > 1)
        error('op_CSIPlot:argumentError', 'No coils to plot')
    end
elseif(in.sz(in.dims.coils) < indecies.coilNum)
    error('op_CSIPlot:argumentError', 'Coil index is larger than number of coils')
end

if(in.dims.averages == 0) 
    if(indecies.averageNum > 1)
        error('op_CSIPlot:argumentError', 'No averages to plot')
    end
elseif(in.sz(in.dims.averages) < indecies.averageNum)
    error('op_CSIPlot:argumentError', 'Average index is larger than number of averages') 
end

if(exist('indecies.x_indecies', 'var') && exist('indecies.x_range', 'var'))
    error('op_CSIPlot:argumentError', 'Only x_indecies or x_range should be used, not both');
end
if(exist('indecies.y_indecies', 'var') && exist('indecies.y_range', 'var'))
    error('op_CSIPlot:argumentError', 'Only y_indecies or y_range should be used, not both');
end


if(~isfield(indecies, 'x_indecies') && ~isfield(indecies, 'x_range'))
    x = [1, length(in.xCoordinates)];
elseif(isfield(indecies, 'x_range'))
    x_val = in.xCoordinates > indecies.x_range(1) & in.xCoordinates < indecies.x_range(2);
    x = [find(x_val, 1 ), find(x_val, 1, 'last' )];
else
    x = indecies.x_indecies;
end

if(~isfield(indecies, 'y_indecies') && ~isfield(indecies, 'y_range'))
    y = [1, length(in.yCoordinates)];
elseif(isfield(indecies, 'y_range'))
    y_val = in.yCoordinates > indecies.y_range(1) & in.yCoordinates < indecies.y_range(2);
    y = [find(y_val, 1 ), find(y_val, 1, 'last' )];
else
    y = indecies.y_indecies;
end



ppmmin = ppm_arguments.min_range;
ppmmax = ppm_arguments.max_range;
%check if ppm exists to plot
if ~isfield(in, 'ppm')
    range_bool = in.t >= ppmmin & in.t <= ppmmax;
    ppm = in.t(range_bool);
    
else
    %if ppm doesn't exist plot time domain
    range_bool = in.ppm >= ppmmin & in.ppm <= ppmmax;
    ppm = in.ppm(range_bool);
    plot_type.plot_dim = 'ppm';
end


%TODO: figure out permutations when there is zero dim. (ie. no coil or
%average dimensions)
[~, non_z_idx, vals] = find([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, ...
                           in.dims.averages]);
specs = permute(in.specs, vals);

idx = {range_bool, x(1):x(2), y(1):y(2), indecies.coilNum, indecies.averageNum};
idx = idx(non_z_idx);
specs = specs(idx{:});


%get lowercase
complex_plot = lower(plot_type.plane_type);

%convert to plot type
if(strcmp(complex_plot, 'real'))
    specs = real(specs);
elseif (strcmp(complex_plot, 'imag'))
    specs = imag(specs);
else
    specs = abs(specs);
end

%max difference in a voxel
yrange = max(max(specs, [], 1) - min(specs, [], 1), [], 'all');


%scale factors to fit at each (x,y) coordinates
scalefactorX=(0.8*in.deltaX)/(ppm(end)-ppm(1));
scalefactorY=(0.8*in.deltaY)/yrange;
specs = specs.*scalefactorY;
%scale the intensity of the specs to the scalefactorY


%create figure and hold 
figure
ax = axes;
hold(ax, 'on');

if(strcmp(plot_type.plot_dim, 'ppm'))
    specs = flip(specs, 1);
end

%set the x and y limits
xlim(ax, [(in.xCoordinates(x(1))-in.deltaX), (in.xCoordinates(x(2))+in.deltaX)]);
xticks(ax, in.xCoordinates(x(1):x(2)));
ylim(ax, [in.yCoordinates(y(1))-in.deltaX, in.yCoordinates(y(2))+in.deltaY]);
yticks(ax, in.yCoordinates(y(1):y(2)));

for x_idx = 1:size(specs,2)
    for y_idx = 1:size(specs,3)
        %first scale the ppm scale so that the range is correct;
        ppm_plot=(ppm-ppm(1))*scalefactorX;
        %now shift it to the correct x-position;
        ppm_plot = ppm_plot + (x_idx-1)*in.deltaX - (0.8*in.deltaX)/2 + in.xCoordinates(x(1));
        y_coords = specs(:,x_idx,y_idx) + (size(specs,3)-y_idx)*in.deltaY + in.yCoordinates(y(1)) - specs(end, x_idx,y_idx);
        %Now start plotting
        plot(ax, ppm_plot, y_coords, 'PickableParts', 'none');        
    end
end
hold(ax, 'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%FUNCTION END%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%






function bound = get_default_range(in, str)
if(isfield(in, 'ppm'))
    range = in.ppm;
else
    range = in.t;
end

if(strcmp(str, 'min'))
    bound = min(range);
else
    bound = max(range);
end
end

function mustBeMemberi(a, members)
    if ~(strcmpi(a, members))
        eidType = 'mustBeMemberi:notAMember';
        msgType = strcat("Value must be a member of this set:  ", strjoin(string(members), ', '));
        throwAsCaller(MException(eidType,msgType))
    end
end
