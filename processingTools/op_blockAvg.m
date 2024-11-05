% op_blockAvg.m
% Jamie Near, Sunnybrook Research Institute 2024.
% 
% USAGE:
% out=op_blockAvg(in,N);
% 
% DESCRIPTION:
% Compress a dataset along the averages dimension by combining the averages 
% into blocks of N averages. The N averages in each block are then averaged, 
% leaving a new dataset with Navg/N transients.  The use cases for this 
% function are 1) for improving the SNR of transients in low SNR data, so
% that retrospective drift correction will work better and 2) for grouping
% averages in a functional MRS acquisition. 
% 
% INPUTS:
% in	= input data in matlab structure format.
% N     = Number of averages to include per block.  
%
% OUTPUTS:
% out   = Output following block averaging.  

function out=op_blockAvg(in,N);

    %Check if there are avea
    if in.flags.averaged || in.averages<2 || in.dims.averages==0
        %DO NOTHING
        disp('WARNING: No averages found. Returning input without modification!');
        out=in;
        return;
    end
    
    %check if the number of averages is an even multiple of N:
    if mod(in.sz(in.dims.averages)/N,1)
        error('ERROR:  The number of averages does not divide evenly into N.');
    end
    
    %Make new fids array:
    if in.dims.averages==2;
        fids=zeros(size(in.fids(:,in.sz(in.dims.averages)/N,:,:,:)));
    elseif in.dims.averages==3;
        fids=zeros(size(in.fids(:,:,in.sz(in.dims.averages)/N,:,:)));
    elseif in.dims.averages==4;
        fids=zeros(size(in.fids(:,:,:,in.sz(in.dims.averages)/N,:)));
    elseif in.dims.averages==5;
        fids=zeros(size(in.fids(:,:,:,:,in.sz(in.dims.averages)/N)));
    end

    %Sum in groups of N:
    for n=1:in.sz(in.dims.averages)/N
        if in.dims.averages==2
            fids(:,n,:,:,:)=sum(in.fids(:,N*(n-1)+1:N*(n-1)+N,:,:,:),2);
        elseif in.dims.averages==3
            fids(:,:,n,:,:)=sum(in.fids(:,:,N*(n-1)+1:N*(n-1)+N,:,:),3);
        elseif in.dims.aveages==4
            fids(:,:,:,n,:)=sum(in.fids(:,:,:,N*(n-1)+1:N*(n-1)+N,:),4);
        elseif in.dims.averages==5
            fids(:,:,:,:,n)=sum(in.fids(:,:,:,:,N*(n-1)+1:N*(n-1)+N),5);
        end
    end

    fids=squeeze(fids);
    fids=fids/N; %divide by number of averages per block;
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
      
    %re-calculate the sz variable
    sz=size(fids);

    %re-caulculate the number of averages:
    averages=in.averages/N;
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.dims=in.dims;
    out.averages=averages;
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;

end



