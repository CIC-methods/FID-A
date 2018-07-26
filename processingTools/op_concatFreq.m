% op_concatFreq.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_concatFreq(in1,in2);
% 
% DESCRIPTION:
% Concatenate two scans along the averages dimension.  Two scans with 50
% averages each will now look like a single scan with 100 averages.
% 
% INPUTS:
% in1    = first input in matlab structure format.
% in2    = second input in matlab structure format.
% shift  = desired shift between stitched spectra in ppm (optional 
%           default=in1.spectralwidth).  This shift must be at least as 
%           large as the spectral witdth. 
%
% OUTPUTS:
% out    = Output following concatenation of inputs along the 
%           averages dimension.  
 
function out=op_concatFreq(in1,in2,shift);
 
%RULE:  Averages dimenstion must be the same for both inputs UNLESS one of
%the inputs has an averages dimension of zero:
 
if in1.dims.t~=0 && in2.dims.t~=0 && in1.dims.t ~= in2.dims.t
    error('time/frequency dimensions must be the same for both inputs');
end
 
%if in1.dims.averages==0 && in2.dims.averages==0
%     dim=2;
%elseif in1.dims.averages==0
%     dim=in2.dims.averages;
%elseif in2.dims.averages==0;
%     dim=in1.dims.averages;
%else
%     dim=in1.dims.averages;
%end
 
if nargin<3 || shift==in1.spectralwidth/(in1.txfrq/1e6)
 
    %concatenante spectra
    specs=cat(1,in1.specs,in2.specs);
    sz=size(specs);
    
    %make new ppm axis
    amountToShift=((in1.ppm(1)-in1.ppm(end))+(in1.ppm(1)-in1.ppm(2)));
    ppmstitch=in1.ppm+amountToShift;
    ppmnew=[ppmstitch in1.ppm];
elseif shift<in1.spectralwidth/(in1.txfrq/1e6)
    error('ERROR:  shift must be at least as large as the spectral width in ppm! ABORTING!!!');
elseif shift>in1.spectralwidth/(in1.txfrq/1e6)
    %make dummy spectrum that is the same spectral width as the desired gap
    dummy=in1;
    if shift>2*(dummy.spectralwidth/(dummy.txfrq/1e6))
        %will do this later
    else 
        gap=shift-in1.spectralwidth/(in1.txfrq/1e6);
        dummy=op_freqrange(dummy,dummy.ppm(1)-gap,dummy.ppm(1)+0.00000001); 
        dummy=op_ampScale(dummy,0); 
    end
    
    %concatenate spectra and gap
    specs=cat(1,in1.specs,dummy.specs);
    specs=cat(1,specs,in2.specs); 
    sz=size(specs);
 
    %make new ppm axis
    amountToShift=((dummy.ppm(1)-dummy.ppm(end))+(dummy.ppm(1)-dummy.ppm(2)));
    ppmstitch=dummy.ppm+amountToShift; 
    ppmnew=[ppmstitch in1.ppm];
    amountToShift=((ppmnew(1)-ppmnew(end))+(ppmnew(1)-ppmnew(2)));
    ppmstitch=in2.ppm+amountToShift;
    ppmnew=[ppmstitch ppmnew]; 
end 
    
 
%convert back to time domain
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(size(specs,in1.dims.t),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids=fft(fftshift(specs,in1.dims.t),[],in1.dims.t);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids=fft(circshift(fftshift(specs,in1.dims.t),1),[],in1.dims.t);
end
 
%Calculate the new spectral width and dwelltime; 
dppm=abs(ppmnew(2)-ppmnew(1));
ppmrange=abs((ppmnew(end)-ppmnew(1)))+dppm; 
spectralwidth=ppmrange*in1.Bo*42.577;
dwelltime=1/spectralwidth;
 
%Calculate the time scale
t=[0:dwelltime:(sz(1)-1)*dwelltime];
 
%FILLING IN DATA STRUCTURE
out=in1;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.t=t;
out.ppm=ppmnew;
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
 
%FILLING IN THE FLAGS
out.flags=in1.flags;
out.flags.writtentostruct=1;
out.flags.averaged=0;