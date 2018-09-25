% op_addScans.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_addScans(in1,in2,subtract);
% 
% DESCRIPTION:
% Add or subtract two scans together.
% 
% INPUTS:
% in1        = First spectrum to add (in matlab structure format)
% in2        = Second spectrum to add (also in matlab structure format).
% subtract	 = (optional).  Add or subtract?  (0 = add, 1=subtract).
%             Default=0;
%
% OUTPUTS:
% out        = Result of adding inputs in1 and in2.  


function out=op_addScans(in1,in2,subtract);

if nargin<3
    subtract=0;
end

%If in1 is an empty structure (in1=struct()), then the output is just in2.  
%This functionality is useful for looping, when you are starting from zero 
%and adding many spectra together.  JN 17/05/2013.
if isempty(fieldnames(in1))
    if subtract==1
        fids=-in2.fids;
    else
        fids=in2.fids;
    end
           
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in2.dims.t),in2.dims.t);
    
    %FILLING IN DATA STRUCTURES
    out=in2;
    out.fids=fids;
    out.specs=specs;

    %FILLING IN THE FLAGS
    out.flags=in2.flags;


else
    if subtract
        fids=in1.fids-in2.fids;
    else
        fids=in1.fids+in2.fids;
    end
    
    if in1.sz ~= in2.sz
        error('ERROR:  Spectra must be the same number of points');
    end
    if in1.dwelltime ~= in2.dwelltime
        error('ERROR:  Spectra must have the same dwelltime');
    end
    if in1.spectralwidth ~= in2.spectralwidth
        error('ERROR:  Spectra must have the same spectral width');
    end
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in1.dims.t),in1.dims.t);
    
    %FILLING IN DATA STRUCTURES
    out=in1;
    out.fids=fids;
    out.specs=specs;

    %FILLING IN THE FLAGS
    out.flags=in1.flags;
end