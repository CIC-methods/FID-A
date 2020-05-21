% sim_plotCSI.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% sim_plotCSI(in,coilNum)
% sim_plotCSI(in)
% 
% DESCRIPTION:
% This function takes a processed MRSI twix file and plots the data using 
% the matlab plot function. Horizontal coordinates are the x coordinates
% and vertical are y coordinates. The plots at specific (x,y) positions are
% the spectral coordintaes
%
% INPUTS:
% in          = input cell array of simulated spectra from a spatially resolved simulation
% coilNum     = coil to plot if coils have not been combined
%
% OUTPUTS:
% displays graph of processed MRSI data in twix struct. 


function sim_plotCSI(in, coilNum)

%check to see if coilNum has been defined
if ~exist('coilNum','var')
    coilNum = 1;
end
 

if(in.dims.coils == 0)
    yrange=max(real(in.specs),[],'all') - min(real(in.specs),[],'all');
else
    yrange=max(real(in.specs(:,coilNum,:,:)),[],'all') - min(real(in.specs(:,coilNum,:,:)),[],'all');
end

%scale factors to fit at each (x,y) coordinates
scalefactorX=(0.8*in.deltaX)/max(in.t);
scalefactorY=(0.8*in.deltaY)/yrange;

%scale the intensity of the specs to the scalefactorY
tempSpec=op_ampScale(in,scalefactorY);

%flip the specs along the y axis 
tempSpec.specs = flip(tempSpec.specs, in.dims.y);

%permute the specs dimensions so t, x, y are the first 3 dimensionss
dimsToPlot = [in.dims.t, in.dims.x, in.dims.y];
extraDims = setdiff(numel(size(in.sz)), dimsToPlot);
tempSpec = permute(tempSpec.specs, [dimsToPlot, extraDims]);
tempSpec = reshape(tempSpec, [in.sz(dimsToPlot), prod(in.sz(extraDims))]);

%create figure and hold 
figure;
hold on;

%set the x and y limits
xlim([(in.xCoordinates(1)-in.deltaX), (in.xCoordinates(end)+in.deltaX)]);
xticks(in.xCoordinates);
ylim([in.yCoordinates(1)-in.deltaX, in.yCoordinates(end)+in.deltaY]);
yticks(in.yCoordinates);

for x = 1:size(in.specs,in.dims.x)
    for y = 1:size(in.specs,in.dims.y)
        %first scale the ppm scale so that the range is correct;
        time=in.t*scalefactorX;
        %now shift it to the correct x-position;
        time = time + (x-1)*in.deltaX - (0.8*in.deltaX)/2 + in.xCoordinates(1);
        %Now start plotting
        plot(time,(tempSpec(:,x,y,coilNum) + ((y-1)*in.deltaY + in.yCoordinates(1))));
    end
end
hold off;
end
