%sim_lcmrawbasis.m
%Jamie Near, McGill University 2014.
%Modified (greatly simplified) for new FID-A spin system definitions, 2018.
%
% USAGE:
% [RF,out]=sim_lcmrawbasis(n,sw,Bfield,linewidth,metab,tau1,tau2,addref,makeraw,seq)
% 
% DESCRIPTION:
% Generate an LCModel .RAW file to be used as an individual metabolite basis 
% spectrum in an LCModel basis set.  The relevant characteristics of the
% acquisition can be specified (pulse sequence, number of points, spectral
% width, etc)
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% tau1      = first echo time in [ms] (if seq='st' or 'l', tau1 = TE)
% tau2      = second echo time in [ms].  (Used in Press, but not used in SE or LASER.
%             If seq='st', tau2=TM).
% addref    = add reference at 0ppm (for use in LCModel makebasis) ['y' or 'n']
% makeraw   = make output file for lcmodel ['y' or 'n']
% seq       = pulse sequence ['se' for Spin Echo, 'p' for Press, 'st' for Steam, or 'l' for LASER]
% metab     = one of the following choices
%   'H2O'    = Water
%   'Ala'    = Alanine
%   'Asp'    = Aspartate
%   'PCh'    = PhosphoCholine
%   'Cr'     = Creatine
%   'PCr'    = PhosphoCreatine
%   'GABA'   = Gamma-aminobutyric acid (kaiser)
%   'Gln'    = Glutamine
%   'Glu'    = Glutamate
%   'GSH'    = Glutathione
%   'Gly'    = Glycine
%   'Ins'    = Myo-inositol
%   'Lac'    = Lactate
%   'NAA'    = N-acetyl aspartate
%   'Scyllo' = Scyllo-inositol
%   'Tau'    = Taurine
%   'Asc'    = Ascorbate (Vitamin C)
%   'bHB'    = beta-Hydroxybutyrate
%   'bHG'    = beta-Hydroxyglutarate
%   'Glc'    = Glucose
%   'NAAG'   = N-acetyl aspartyl glutamate
%   'GPC'    = Glycero-phosphocholine
%   'PE'     = Phosphoryl ethanolamine
%   'Ser'    = Serine
%
% OUTPUTS:
% RF        = not used.
% out       = Simulated basis spectrum in FID-A structure format.  

function [RF,out]=sim_lcmrawbasis(n,sw,Bfield,linewidth,metab,tau1,tau2,addref,makeraw,seq)


load('spinSystems.mat');

eval(['sys=sys' metab ';']);
spins=0;
for k=1:length(sys)
    spins=spins+length(sys(k).shifts);
end
disp(['simulating metabolite ' metab ' with ' num2str(spins) ' spins...  Please Wait...']);
switch seq
    case 'se'
        out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
    case 'p'
        out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
    case 'st'
        out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
    case 'l'
        out = sim_laser(n,sw,Bfield,linewidth,sys,tau1);
    otherwise
        disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
end

if addref=='y'||addref=='Y'
    addref=1;
elseif addref=='n'||addref=='N'
    addref=0;
end

if addref
    sysRef.J=0;
    sysRef.shifts=0;
    sysRef.name='Ref_0ppm';
    sysRef.scaleFactor=1;
    switch seq
        case 'se'
            ref = sim_spinecho(n,sw,Bfield,linewidth,sysRef,tau1);
        case 'p'
            ref = sim_press(n,sw,Bfield,linewidth,sysRef,tau1,tau2);
        case 'st'
            ref = sim_steam(n,sw,Bfield,linewidth,sysRef,tau1,tau2);
        case 'l'
            ref = sim_laser(n,sw,Bfield,linewidth,sysRef,tau1);
        otherwise
            disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
    end
    out=op_addScans(out,ref);
end

if makeraw=='y'||makeraw=='Y'
    makeraw=1;
elseif makeraw=='n'||makeraw=='N'
    makeraw=0;
end

if makeraw
    RF=io_writelcmraw(out,[metab '.RAW'],metab);
else
    RF=[];
end


   