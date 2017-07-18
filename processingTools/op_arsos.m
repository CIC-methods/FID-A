% op_arsos.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_arsos(in,domain);
% 
% DESCRIPTION:
% Perform all rank statistic order filter (see Slotboom J et al, Meas Sci
% Technol.  20 (2009)).  This effectively rank-orders the data along the
% averages dimension.  This can be done in either the time domain (default)
% or the frequency domain.  
% 
% INPUTS:
% in	 = input data in matlab structure format.
% domain = time domain ('t', (default)) or frequency domain ('f').
%
% OUTPUTS:
% out    = Output following arsos filtering.

function out=op_arsos(in,domain);

if in.flags.averaged || in.dims.averages==0 || in.averages<2
    error('ERROR:  Averaging has already been performed!  Aborting!');
end

if nargin<2
    domain='t';
end

%sort the data along the averages dimension;
if strcmp(domain,'t');
    fids=sort(real(in.fids),in.dims.averages)+(1i*sort(imag(in.fids),in.dims.averages));
    
    %Do a DC shift correction using the last 25% of the time domain data:
    tminIndex=ceil(0.75*length(in.t));
    
    medianReal=median(real(fids),in.dims.averages);
    medianImag=median(imag(fids),in.dims.averages);
    
    for n=1:in.sz(in.dims.averages)
        if in.dims.averages==2
            DC0real=mean(real(fids(tminIndex:end,n)))-mean(medianReal(tminIndex:end));
            DC0imag=mean(imag(fids(tminIndex:end,n)))-mean(medianImag(tminIndex:end));
            fids(:,n)=fids(:,n)-((DC0real+(1i*DC0imag)));
        end
    end
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
elseif strcmp(domain,'f');
    specs=sort(real(in.specs),in.dims.averages)+(1i*sort(imag(in.specs),in.dims.averages));
    
    %Do a DC shift correction using everything below 0 ppm and everything 
    %above 10 ppm:
    medianReal=median(real(specs),in.dims.averages);
    medianImag=median(imag(specs),in.dims.averages);
    
    for n=1:in.sz(in.dims.averages)
        if in.dims.averages==2
            DC0real=mean(real(specs(in.ppm<0 & in.ppm>10,n)))-mean(medianReal(in.ppm<0 & in.ppm>10));
            DC0imag=mean(imag(specs(in.ppm<0 & in.ppm>10,n)))-mean(medianImag(in.ppm<0 & in.ppm>10));
            specs(:,n)=specs(:,n)-((DC0real+(1i*DC0imag)));
        end
    end
    
    
    %convert back to time domain
    %if the length of Fids is odd, then you have to do a circshift of one to
    %make sure that you don't introduce a small frequency shift into the fids
    %vector.
    if mod(size(specs,in.dims.t),2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids=fft(fftshift(specs,in.dims.t),[],in.dims.t);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids=fft(circshift(fftshift(specs,in.dims.t),1),[],in.dims.t);
    end
end

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;


