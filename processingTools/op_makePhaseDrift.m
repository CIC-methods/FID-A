% op_makePhaseDrift.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,phDrift]=op_makePhaseDrift(in,totalDrift,noise);
% 
% DESCRIPTION:
% Add phase drift to a dataset containing multiple averages.  This is generally 
% used to generate simulated datasets with phase drift. 
% 
% INPUTS:
% in         = input data in matlab structure format.
% totalDrift = total amount of phase drift (in degrees) to add over the whole scan.
%               If totalDrift is a scalar, then a constant slope of drift
%               will be added.  If totalDrift is a vector or matrix with dimensions
%               equal to the dimensions of the input data, then this
%               vector specifies the drift applied to each average.
% noise      = the standard deviation of noise to add to the phase drift
%             function.
%
% OUTPUTS:
% out        = Output dataset with phase drift added.
% phDrift     = Vector of phase drift values that were added (in degrees).


function [out,phDrift]=op_makePhaseDrift(in,totalDrift,noise);
%out=op_makedrift(in,totalDrift);

%First make the matrices needed for multiplication
if any(totalDrift)
    if isequal(size(totalDrift),[in.averages/in.subspecs,in.subspecs])
        ph=totalDrift;
    else
        ph=linspace(0,totalDrift,in.averages*in.subspecs);
        ph=reshape(ph,in.subspecs,in.averages);
        ph=ph';
    end
else
    ph=zeros(in.averages,in.subspecs);
end
phnoise=noise*randn(size(ph));
phDrift=ph+phnoise;
%PH=repmat(phDrift',in.sz(1),1);
phDrift_shft=shiftdim(phDrift,-1);
PH=repmat(phDrift_shft,in.sz(1),1,1);

%Now apply the drift to the fids;
fids=in.fids.*exp(-1i*PH*pi/180);

%Now re-calculate specs using ifft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
