% op_getcoilcombos_specReg.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% coilcombos=op_getcoilcombos_specReg(file_or_struct,tmin,tmax,point);
% 
% DESCRIPTION:
% This funciton finds the relative coil phases and amplitudes.  Coil phases
% are found by fitting in the time domain. (In the original
% op_getcoilcombos, the phases were determined simply by observation of the
% pointth point in the time domain).  The "Base" receiver channel is
% Determined as the one with the highest signal (all other coil channels
% will be registered to that channel.  
% 
% INPUTS:
% file_or_struct     = this function will accept either a string filename or
%                     the name of a structure.  If the input is a string, 
%                     the program will read in the data corresponding to 
%                     that filename.  If the input is a structure, it will
%                     operate on that structure.
% tmin               = The earliest timepoint in the fid to be used for 
%                     spectral registration.  (Optional.  Default=0 sec);
% tmax               = The latest timepoint in the fid to be used for 
%                     spectral registration.  (Optional.  Default=0.2 sec);
% point              = The index of the datapoint in the fid that is used
%                     for determination of Signal intensity. (Optional.
%                     Default = 1);
%
% OUTPUTS:
% coilcombos         = Structure containing the calculated coil weights and phases. 


function coilcombos=op_getcoilcombos_specReg(file_or_struct,tmin,tmax,point);


if isstr(file_or_struct)
    in=io_loadspec_twix(file_or_struct);
else
    in=file_or_struct;
end

if in.flags.addedrcvrs
    error('ERROR:  must provide data prior to coil combination!!  ABORTING!!');
end

if nargin<4
    point=1;
    if nargin<3
        tmax=0.2;
        if nargin<2
            tmin=0;
        end
    end
end

B=in.sz(in.dims.coils);

coilcombos.ph=zeros(B,1);
coilcombos.sig(:,1)=abs(in.fids(point,:,1,1));
bestSNRindex=find(coilcombos.sig(:,1)==max(coilcombos.sig(:,1)))
phGuess=0;

disp('aligning all coils to the first coil');
base=phase(in.fids(in.t>=tmin & in.t<tmax,bestSNRindex,1,1));
begin=1;
for n=begin:B
    %disp(['Fitting coil number ' num2str(n) 'to first coil']);
    phFit=nlinfit(in.fids(in.t>=tmin & in.t<tmax,n,1,1)/coilcombos.sig(n,1),base,@op_phaseShiftRealNest,phGuess);
    coilcombos.ph(n,1)=(-1*phFit*pi/180);
%     subplot(1,2,1);
%     plot(in.ppm,in.specs(:,1,1,1)/coilcombos.sig(1,1),in.ppm,fftshift(ifft(addphase(in.fids(:,n,1,1),phFit)))/coilcombos.sig(n,1));xlim([4 5.5]);
%     subplot(1,2,2);
%     plot(in.t,in.fids(:,1,1,1)/coilcombos.sig(1,1),in.t,addphase(in.fids(:,n,1,1),phFit)/coilcombos.sig(n,1));xlim([0 tmax*2]);
%     pause;
    
end

coilcombos.ph=coilcombos.ph+(phase(in.fids(point,bestSNRindex,1,1)));
%Now normalize the coilcombos.sig so that the max amplitude is 1;
coilcombos.sig=coilcombos.sig/max(coilcombos.sig);


    function y=op_phaseShiftRealNest(pars,input);
        p=pars(1);     %Phase Shift [deg]
        
        fid=input(:);
        y=phase(addphase(fid,p));
          
    end
end


