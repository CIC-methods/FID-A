% op_peakFit.m
% Jamie Near, McGill University 2017.
% 
% USAGE:
% [fit,parsFit,residual]=op_peakFit(in,ppmmin,ppmmax,parsGuess);
% 
% DESCRIPTION:
% Perform a Voigt lineshape fit to a prominent singlet resonance within the
% provided ppm window of a proton MRS dataset.
% 
% INPUTS:
% in	    = input data in matlab stucture format.
% ppmmin	= lower frequency limit for fitting.
% ppmmax	= upper frequency limit for fitting.
% parsGuess = Initial parameter guesses:
%               (Amplutide [arb units]
%                linewidth [Hz]
%                frequency [ppm]
%                phase [degrees])
%
% OUTPUTS:
% fit       = Voigt lineshape fit in simplified FID-A structure format. 
% parsFit   = Fit parameters (same format as parsGuess).
% residual  = The fit residual in simplified FID-A structure format.

function [fit,parsFit,residual]=op_peakFit(in,ppmmin,ppmmax,parsGuess);

%if in.flags.isISIS
%    error('ERROR:  must have combined subspecs in order to do this!  ABORTING');
%end

% if ~in.flags.averaged
%     error('ERROR:  must have averaged in order to do this!  ABORTING');
% end

% if ~in.flags.addedrcvrs
%     error('ERROR:  must have added receivers in order to do this!  ABORTING');
% end

in_range=op_freqrange(in,ppmmin,ppmmax);
specs=in_range.specs;

ppm=in_range.ppm;
spec=in_range.specs;

plot(ppm,spec);

if nargin<4
    parsGuess=zeros(1,4);
    parsGuess(1)=max(abs(spec)); %Amplitude
    parsGuess(2)=5;
    parsGuess(3)= in_range.ppm(abs(in_range.specs)==max(abs(in_range.specs))); %Frequency [ppm];
    %parsGuess(4)=(real(spec(1))-real(spec(end)))/(ppm(1)-ppm(end));  %Baseline slope
    %parsGuess(5)=real(spec(1))-(parsGuess(:,4)*ppm(1));   %Baseline Offset
    parsGuess(4)=0; %phase [degrees];
end




yGuess=op_voigt_linbas_complex_nest(parsGuess,ppm);
parsFit=nlinfit(ppm,real(spec'),@op_voigt_linbas_real_nest,parsGuess);
yFit=op_voigt_linbas_complex_nest(parsFit,ppm);

plot(ppm,spec,ppm,yGuess,':',ppm,yFit);
legend('data','guess','fit');
parsFit;
parsGuess;

area=parsFit(1).*parsFit(2);
for n=1:size(area,1)
    disp(['Area under the ' num2str(n) 'th fitted curve is: ' num2str(area(n))]);
end
area=sum(area);
disp(['Area under the fitted curve is: ' num2str(area)]);

fit = in;
fit.specs=op_voigt_linbas_complex_nest(parsFit,in.ppm);
fit.specs=fit.specs.';
if mod(size(fit.specs,1),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fit.fids=fft(fftshift(fit.specs,1),[],1);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fit.fids=fft(circshift(fftshift(fit.specs,1),1),[],1);
end
residual=op_subtractScans(in,fit);

    function [y]=op_voigt_linbas_real_nest(pars,ppm)
        
        A=pars(1);     %Amplitude  (Use only positive Numbers)
        lw=pars(2);     %full width at half max
        ppm0=pars(3);  %centre frequency of peak
        
        %figure out time domain params
        deltappm=abs(ppm(1)-ppm(2));  %Frequnecy resolution in [ppm]
        deltaf=deltappm*42.577*in.Bo;  %Frequnecy resolution in [Hz]
        tacq=1/deltaf; %FID duration in [sec]
        t=linspace(0,tacq,length(ppm));  %Time scale in [sec]
        t2=1/(pi*lw);  %T2 relaxation constant [sec]
        
        %figure out f0 relative to the new frequency window:
        centreppm=(ppm(1)+ppm(end))/2;  
        f0=(ppm0-centreppm) * 42.577 * in.Bo;
        
        %m=pars(4);    %baseline slope
        %y0=pars(5);    %baseline offset
        theta=pars(4); %Phase shift of peak
        gam=lw/42.577/in.Bo;      %gamma parameter [ppm]
        %Define a complex lorentzian for each set of pars;
        y=zeros(length(A),length(ppm));
        
        %lor = sqrt(2/pi) *(gam(n)+i*(ppm-ppm0(n))) ./ (gam(n)^2 + (ppm-ppm0(n)).^2);
        fid = exp(-t/t2) .* exp(1i * f0 * 2 * pi * t);
        lor = fftshift(ifft(fid));
        %now scale it, add baseline, phase it by theta, and take the real part
        y=real(addphase(lor/max(abs(lor))*A,theta));
        
        
    end
    function [y]=op_voigt_linbas_complex_nest(pars,ppm)
        
        A=pars(1);     %Amplitude  (Use only positive Numbers)
        lw=pars(2);     %full width at half max
        ppm0=pars(3);  %centre frequency of peak
        
        %figure out time domain params
        deltappm=abs(ppm(1)-ppm(2));  %Frequnecy resolution in [ppm]
        deltaf=deltappm*42.577*in.Bo;  %Frequnecy resolution in [Hz]
        tacq=1/deltaf; %FID duration in [sec]
        t=linspace(0,tacq,length(ppm));  %Time scale in [sec]
        t2=1/(pi*lw); %T2 relaxation constant [sec];
        
        %figure out f0 relative to the new frequency window:
        centreppm=(ppm(1)+ppm(end))/2;  %centre of ppm window [ppm]
        f0=(ppm0-centreppm) * 42.577 * in.Bo;  %frequency offset of peak [Hz]
        
        %m=pars(4);    %baseline slope
        %y0=pars(5);    %baseline offset
        theta=pars(4); %Phase shift of peak
        gam=lw/42.577/in.Bo;      %gamma parameter [ppm]
        %Define a complex lorentzian for each set of pars;
        y=zeros(length(A),length(ppm));
        
        fid = exp(-t/t2) .* exp(1i * f0 * 2 * pi * t);
        lor = fftshift(ifft(fid));
        y=addphase(lor/max(abs(lor))*A,theta);
        
            
        
    end

end




