% op_lorentz.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% y=op_lorentz(pars,ppm);
% 
% DESCRIPTION:  
% Generate a parametrized lorentzian peak.
% This function is fed into fitting tools to enable fitting of peaks using
% lorentzian lineshapes.
% 
% INPUTS:
% pars       = The parameters of the Lorentizan function.  This is a five 
%             element vector consisting of the following fields:
%             [ Amplitude, 
%               FWHM, (In Hz) 
%               Centre Freq,  (In ppm)
%               baseline offset, (in Amplitude units)
%               Phase shift];  (in degrees)
% ppm        = frequency axis vector (in ppm);
%
% OUTPUTS:
% y          = Output vector specifying a lorentzian lineshape.  



function [y]=op_lorentz(pars,ppm)

A=pars(:,1);     %Amplitude  (Use only positive Numbers)
w=pars(:,2);     %full width at half max
ppm0=pars(:,3);  %centre frequency of peak

if size(pars,2)==5
    y0=pars(:,4);    %baseline offset
    theta=pars(:,5); %Phase shift of peak
    gam=w/2;      %gamma parameter
    %Define a complex lorentzian for each set of
    y=zeros(length(A),length(ppm));
    for n=1:length(A)
        y(n,:) = sqrt(2/pi) *(gam(n)-i*(ppm-ppm0(n))) ./ (gam(n)^2 + (ppm-ppm0(n)).^2);
        %now scale it, add baseline, phase it by theta, and take the real part
        y(n,:)=real(addphase(y(n,:)/max(abs(y(n,:)))*A(n)+y0(n),theta(n)));
    end
    y=sum(y,1);
elseif size(pars,2)==4
    y0=pars(:,4);    %baseline offset
    gam=w/2;      %gamma parameter
    %assume zero phase:
    theta=zeros(size(pars,1),1);
    %Define a complex lorentzian for each set of
    y=zeros(length(A),length(ppm));
    for n=1:length(A)
        y(n,:) = sqrt(2/pi) *(gam(n)-i*(ppm-ppm0(n))) ./ (gam(n)^2 + (ppm-ppm0(n)).^2);
        %now scale it, add baseline, phase it by theta, and take the real part
        y(n,:)=real(addphase(y(n,:)/max(abs(y(n,:)))*A(n)+y0(n),theta(n)));
    end
    y=sum(y,1);
elseif size(pars,2)==3
    gam=w/2;      %gamma parameter
    %assume zero phase:
    theta=zeros(size(pars,1),1);
    %assume zero baseline offset:
    y0=zeros(size(pars,1),1);
    %Define a complex lorentzian for each set of
    y=zeros(length(A),length(ppm));
    for n=1:length(A)
        y(n,:) = sqrt(2/pi) *(gam(n)-i*(ppm-ppm0(n))) ./ (gam(n)^2 + (ppm-ppm0(n)).^2);
        %now scale it, add baseline, phase it by theta, and take the real part
        y(n,:)=real(addphase(y(n,:)/max(abs(y(n,:)))*A(n)+y0(n),theta(n)));
    end
    y=sum(y,1);
end




    

