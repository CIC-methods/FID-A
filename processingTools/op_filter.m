% op_filter.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,lor]=op_filter(in,lb);
% 
% DESCRIPTION:
% Perform line broadening by multiplying the time domain signal by an
% exponential decay function.  
% 
% INPUTS:
% in     = input data in matlab structure format.
% lb     = line broadening factor in Hz.
%
% OUTPUTS:
% out    = Output following alignment of averages.  
% lor    = Exponential (time domain) filter envelope that was applied.

function [out,lor]=op_filter(in,lb);

if lb==0
    out=in;
else
    
    if in.flags.filtered
        %cont=input('WARNING:  Line Broadening has already been performed!  Continue anyway?  (y or n)','s');
        cont='y';
        if cont=='y'
            %continue;
        else
            error('STOPPING');
        end
    end
    
    fids=in.fids;
    
    t2=1/(pi*lb);
    
    %Create an exponential decay (lorentzian filter):
    lor=exp(-in.t/t2);
    %plot(in.t,lor);
    
    %first make a bunch of vectors of ones that are the same lengths as each of
    %the dimensions of the data.  Store them in a cell array for ease of use.
    for n=1:length(in.sz)
        p{n}=ones(in.sz(n),1);
    end
    
    %Now, now take the lorentzian filter vector that we made earlier (lor) and use it
    %to populate a filter array that has the same dimensions as the data.  To
    %do this, we have to use the ndgrid function, which is essentially the same
    %as the meshgrid function, except in multiple dimensions.  b, c, and d are
    %dummy variables and are not used.
    if length(in.sz)==1
        fil=lor;
    end
    if length(in.sz)==2
        [fil,b]=ndgrid(lor,p{2});
    end
    if length(in.sz)==3
        [fil,b,c]=ndgrid(lor,p{2},p{3});
    end
    if length(in.sz)==4
        [fil,b,c,d]=ndgrid(lor,p{2},p{3},p{4});
    end
    if length(in.sz)==5
        [fil,b,c,d,e]=ndgrid(lor,p{2},p{3},p{4},p{5});
    end
    
    %Now multiply the data by the filter array.
    fids=fids.*fil;
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.filtered=1;
    
end
