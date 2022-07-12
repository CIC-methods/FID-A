%op_CSISpecReg.m
%
% USAGE:
% [out, freqs, phases]=op_CSISpecReg(in,averagingMode,pointsPerLoop);
%
% DESCRIPTION:
% Takes a CSI FID-A structure and applies frequency and phase drift correction.
% Uses spectral registration to align the FIDS acquired from a common
% k-space point. Applies the measured frequency and phase drift correction
% to the entire k-space. This requires an MRSI sequence that has a common
% k-space point at every TR (eg. This function will work for a rosette
% sequence but not EPSI/Concentric rings/standard CSI).
%
% INPUTS:
% in        = input data in FID-A CSI structure format.
%
% OUTPUTS:
% out      = corrected output data in FID-A CSI structure format
% freqs    = measured frequency drift
% phases    = measured phase drift


function [out,freqs,phases,fid0,fid0_fc]=op_CSISpecReg(in,averagingMode,pointsPerLoop)
% arguments
%     in (1,1) struct
% end

if in.dims.extras
    error('ERROR:  this function cannot be run on a structure that contains an ''extras'' dimension.');
end

new=in;

if(new.dims.averages==1||new.dims.averages==0)
    fid0=op_CSImakeK0fid(new,1,pointsPerLoop);
    out=in;
    freqs=zeros(fid0.sz(fid0.dims.coils),1);
    phases=zeros(fid0.sz(fid0.dims.coils),1);

end



if(new.dims.averages>1)
    
    if (strcmp(averagingMode,'longTerm'))
        Data4D=new.data;
        %Reshapes the 4D Data to 3D data. Keeps the first averages shots
        %first and keep others followed by. 
        Data3D=zeros(size(Data4D,1),size(Data4D,2),size(Data4D,3)*size(Data4D,4));
        for i=1:size(Data4D,3)
            Data3D(:,:,1+(size(Data4D,4)*(i-1)):size(Data4D,4)*i)=Data4D(:,:,i,:);
        end
        new.data=Data3D;
        clear Data3D;
        new.sz=size(new.data);
        new.dims.ky=0;
        new.dims.averages=3;
        new.averages=new.sz(3);
        new.rawAverages=new.sz(3);
        
    end

     if (strcmp(averagingMode,'shortTerm'))
         new.data=permute(new.data,[1,2,4,3]);
         new.dims.ky=3;
         new.dims.averages=4;
         new.data=reshape(new.data, [new.sz(1),new.sz(2),new.sz(3)*new.sz(4)]);
         new.sz=size(new.data);

        %rename the dimensions.  Call all shots averages:
        new.dims.ky=0;
        new.dims.averages=3;
        new.averages=new.sz(3);
        new.rawAverages=new.sz(3);
   
     end

    %Extract the k0 fids:
      fid0=op_CSImakeK0fid(new,1,pointsPerLoop);

    %Making measures across all coils and all averages/Nsh:
      freqs=zeros(fid0.sz(fid0.dims.coils),fid0.sz(fid0.dims.averages));
      phases=zeros(fid0.sz(fid0.dims.coils),fid0.sz(fid0.dims.averages));

     %Now loop through all coils and calculate the frequency/phase drift
     %for each channel and shot. 
        for n=1:fid0.sz(fid0.dims.coils)
            fid0_coiln=op_takecoils(fid0,n);
            [temp{n},f,p]=op_CSIalignAverages(fid0_coiln,0.5, 'y');
            freqs(n,:)=f;
            phases(n,:)=p;
        end
        figure;
        plot([1:fid0.sz(fid0.dims.averages)],freqs);
        figure;
        plot([1:fid0.sz(fid0.dims.averages)],phases);


        %Now that we have estimated the frequency and phase offsets for
        %each FID, we want to apply the frequency and phase drift
        %corrections.

       
        %First make a t matrix for multiplication:
        t_mat=repmat(new.adcTime',[1,new.sz(new.dims.coils),fid0.sz(fid0.dims.averages)]);

        %Now make frequency and phase matrices for multiplcation:
        f_mat=repmat(freqs,[1,1,new.sz(new.dims.t)]);
        f_mat=permute(f_mat,[3,1,2]);
        p_mat=repmat(phases,[1,1,new.sz(new.dims.t)]);
        p_mat=permute(p_mat,[3,1,2]);


        %Now do the frequency and phase drift correction:
        mrsi_fc=new;
        mrsi_fc.data=mrsi_fc.data.*exp(1i*t_mat.*f_mat*2*pi).*exp(1i*p_mat*pi/180);

        %extract the corrected k0 fids;
        fid0_sr=op_CSImakeK0fid(mrsi_fc,1,pointsPerLoop);

   if (strcmp(averagingMode,'longTerm'))

       % reshaping 3D data to 4D for output
    
       out=mrsi_fc;
       Out_Data3D=out.data;
    
       New4D=zeros(size(in.data));
    
       for i=1:size(New4D,3)
          New4D(:,:,i,:)  = Out_Data3D(:,:,1+(size(New4D,4)*(i-1)):size(New4D,4)*i);
       end
    
       out.data=New4D;
       out.sz=size(New4D);
       clear New4D
       %Update fields in output structure:
        out.dims.ky=4;
        out.dims.averages=3;
        out.averages=out.sz(out.dims.averages);
        out.rawAverages=out.sz(out.dims.averages);
   end 

   if(strcmp(averagingMode,'shortTerm'))
       %Now reshape mrsi_fc so that we have the separate averages and
        %shots dimensions:
        out=mrsi_fc;
        out.data=reshape(out.data,[in.sz(in.dims.t),in.sz(in.dims.coils),in.sz(in.dims.averages),in.sz(in.dims.ky)])
        out.sz=size(out.data);

        %Update fields in output structure:
        out.dims.ky=4;
        out.dims.averages=3;
        out.averages=out.sz(out.dims.averages);
        out.rawAverages=out.sz(out.dims.averages);

   end
    
  end
end

function [out,fs,phs]=op_CSIalignAverages(in,tmax,med,initPars)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.averages==0||in.dims.averages==1
    %DO NOTHING
    disp('WARNING:  No averages found.  Returning input without modification!');
    out=in;
    fs=0;
    phs=0;

else

    parsFit=[0,0];
    
    if nargin<3
        med='n'
        if nargin<2
            %if tmax is not specified, find the time at which the SNR
            %drops below 5
            disp('tmax not supplied.  Calculating tmax....');
            sig=abs(in.fids);
            noise=std(real(in.fids(ceil(0.75*end):end,:,:)),[]);
            noise=mean(mean(mean(noise,2),3),4);
            snr=sig/noise;
            
            for n=1:(numel(snr)/size(snr,1))
                N=find(snr(:,n)>5);
                tmax_est(n)=in.t(N(end));
            end
            tmax=median(tmax_est);
            disp(['tmax = ' num2str(tmax*1000) 'ms.']);
        end
    end
    
    if (strcmp(med,'r') || strcmp(med,'R'))
        if nargin<4
            error('ERROR:  If using the ''r'' option for input variable ''med'', then a 4th input argument must be provided');
        end
    else
        if nargin<4
            ref=struct();
        end
    end
    
    if in.dims.subSpecs==0
        B=1;
    else
        B=in.sz(in.dims.subSpecs);
    end
    
    fs=zeros(in.sz(in.dims.averages),B);
    phs=zeros(in.sz(in.dims.averages),B);
    fids=zeros(in.sz(in.dims.t),1,B);
    for m=1:B
        if med=='y' || med=='Y'
            disp('Aligning all averages to the median of the averages.');
            base=op_median(in);
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        elseif med=='a' || med=='A'
            disp('Aligning all averages to the average of the averages.');
            base=op_averaging(in);
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        elseif med=='n' || med=='N'
            %Use the first average as the base:
            disp('Aligning all averages to the first average.');
            base=[real(in.fids(in.t>=0 & in.t<tmax,1,m));imag(in.fids(in.t>=0 & in.t<tmax,1,m))];
            fids(:,1,m)=in.fids(:,1,m);
            ind_min=1;
        elseif med=='r' || med=='R'
            disp('Aligning all averages to an externally provided reference spectrum.');
            base=ref;
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        end
        for n=1:in.sz(in.dims.averages)
            if n~=ind_min
                parsGuess=parsFit;
                %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
                parsFit=nlinfit(in.fids(in.t>=0 & in.t<tmax,n,m),base,@op_freqPhaseShiftComplexNest,parsGuess);
                fids(:,n,m)=op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
                fs(n,m)=parsFit(1);
                phs(n,m)=parsFit(2);
                %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
            end
        end
    end
    
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.freqcorrected=1;
    
end


    function y=op_freqPhaseShiftComplexNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        shifted=addphase(fid.*exp(1i*t'*f*2*pi),p);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        y=addphase(fid.*exp(1i*t'*f*2*pi),p);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

end






            



