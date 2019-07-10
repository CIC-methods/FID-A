% op_makeFreqDrift.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fDrift]=op_makeFreqDrift(in,totalDrift,noise);
% 
% DESCRIPTION:
% Add frequency drift to a dataset containing multiple averages.  This is 
% generally used to generate simulated datasets with phase drift. 
% 
% INPUTS:
% in         = input data in matlab structure format.
% totalDrift = total amount of frequency drift (in Hz) to add over the whole scan.
%               If totalDrift is a scalar, then a constant slope of drift
%               will be added.  If totalDrift is a vector or matrix with dimensions 
%               equal to the dimensions of the input data, then this
%               vector specifies the drift applied to each average.
% noise      = the standard deviation of noise to add to the frequency drift
%             function.
%
% OUTPUTS:
% out        = Output dataset with frequency drift added.
% fDrift     = Vector of frequency drift values that were added (in Hz).

function [out,fDrift]=op_makeFreqDrift(in,totalDrift,noise);
%

%First make the matrices needed for multiplication
T=repmat(in.t',[1 in.sz(2:end)]);
if any(totalDrift)
    if isequal(size(totalDrift),[in.averages/in.subspecs,in.subspecs])
        f=totalDrift;
    else
        f=linspace(0,totalDrift,in.averages*in.subspecs);
        f=reshape(f,in.subspecs,in.averages);
        f=f';
    end
else
    f=zeros(in.averages,in.subspecs);
end
fnoise=noise*randn(size(f));
fDrift=f+fnoise;
%F=repmat(fDrift',in.sz(1),1);
fDrift_shft=shiftdim(fDrift,-1);
F=repmat(fDrift_shft,in.sz(1),1,1);

%Now apply the drift to the fids;
fids=in.fids.*exp(1i*T.*F*2*pi);

%Now re-calculate specs using ifft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;




