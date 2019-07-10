% run_simExampleBasisSet.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% Script to generate simulated basis spectra for all metabolites of interest
% in the human brain.  The script will generate an LCModel format .RAW file 
% for each metabolite basis spectrum, which can then be passed into 
% LCModel's "makebasis" function to generate a complete LCModel basis set.
% 
% INPUTS:
% To run this script, edit the following parameters as desired and then click
% run:
%  lb       = linewidth (Hz)
%  np       = Spectral points
%  sw       = Spectral width (Hz)
%  Bo       = Magnetic Field Strength (Tesla)
%  te1      = First PRESS echo time, or SPECIAL echo time (ms)
%  te2      = Second PRESS echo time (if applicable) (ms).
%  seq      = Pulse sequence ('se'= SPECIAL, 'p'=press, 'st'=steam);
%  ref      = Add reference peak at 0ppm (used in LCModel, y or n);
%
% OUTPUTS:
% H2O       = Simulated water spectrum
% Ala       = Simulated alanine spectrum
% Asp       = Simulated aspartate spectrum
% PCh       = Simulated phosphocholine spectrum
% Cr        = Simulated creatine spectrum
% PCr       = Simulated phosphochreatine spectrum
% GABA      = Simulated GABA spectrum
% Gln       = Simulated glutamine spectrum
% Glu       = Simulated glutamate spectrum
% GSH       = Simulated glutathione spectrum
% Gly       = Simulated glycine spectrum
% Ins       = Simulated myo-inositol spectrum
% Lac       = Simulated Lactate spectrum
% NAA       = Simulated NAA spectrum
% Scyllo    = Simulated scyllo-inositol spectrum
% Tau       = Simulated taurine spectrum
% Asc       = Simulated ascorbate spectrum
% bHB       = Simulated beta-hydroxybutyrate spectrum
% bHG       = Simulated 2-hydroxyglutyrate spectrum
% Glc       = Simulated glucose spectrum
% NAAG      = Simulated N-acetylaspartylglutamate spectrum
% GPC       = Simulated glycerophosphocholine spectrum
% PE        = Simulated phosphoethanolamine spectrum
% Ser       = Simulated serine spectrum
% EtOH      = Simulated ethanol spectrum
%
% 
% ************INPUT PARAMETERS**********************************
 lb=2;         %linewidth (Hz)
 np=8192;      %Spectral points
 sw=4000;      %Spectral width (Hz)
 Bo=7;         %Magnetic Field Strength (Tesla)
 te1=10;       %First PRESS echo time, or SPECIAL echo time (ms)
 te2=125;      %Second PRESS echo time (if applicable) (ms).
 seq='p'       %Pulse sequence ('se'= SPECIAL, 'p'=press, 'st'=steam, 'l'=laser);
 ref='n'       %Add reference peak at 0ppm (used in LCModel, y or n);
% *************END OF INPUT PARAMETERS**************************

    [RF,H2O]=sim_lcmrawbasis(np,sw,Bo,lb,'H2O',te1,te2,ref,'y',seq);
    [RF,Ala]=sim_lcmrawbasis(np,sw,Bo,lb,'Ala',te1,te2,ref,'y',seq);
    [RF,Asp]=sim_lcmrawbasis(np,sw,Bo,lb,'Asp',te1,te2,ref,'y',seq);
    [RF,PCh]=sim_lcmrawbasis(np,sw,Bo,lb,'PCh',te1,te2,ref,'y',seq);
    [RF,Cr]=sim_lcmrawbasis(np,sw,Bo,lb,'Cr',te1,te2,ref,'y',seq);
    [RF,PCr]=sim_lcmrawbasis(np,sw,Bo,lb,'PCr',te1,te2,ref,'y',seq);
    [RF,GABA]=sim_lcmrawbasis(np,sw,Bo,lb,'GABA',te1,te2,ref,'y',seq);
    [RF,Gln]=sim_lcmrawbasis(np,sw,Bo,lb,'Gln',te1,te2,ref,'y',seq);
    [RF,Glu]=sim_lcmrawbasis(np,sw,Bo,lb,'Glu',te1,te2,ref,'y',seq);
    [RF,GSH]=sim_lcmrawbasis(np,sw,Bo,lb,'GSH',te1,te2,ref,'y',seq);
    [RF,Gly]=sim_lcmrawbasis(np,sw,Bo,lb,'Gly',te1,te2,ref,'y',seq);
    [RF,Ins]=sim_lcmrawbasis(np,sw,Bo,lb,'Ins',te1,te2,ref,'y',seq);
    [RF,Lac]=sim_lcmrawbasis(np,sw,Bo,lb,'Lac',te1,te2,ref,'y',seq);
    [RF,NAA]=sim_lcmrawbasis(np,sw,Bo,lb,'NAA',te1,te2,ref,'y',seq);
    [RF,Scyllo]=sim_lcmrawbasis(np,sw,Bo,lb,'Scyllo',te1,te2,ref,'y',seq);
    [RF,Tau]=sim_lcmrawbasis(np,sw,Bo,lb,'Tau',te1,te2,ref,'y',seq);
    [RF,Asc]=sim_lcmrawbasis(np,sw,Bo,lb,'Asc',te1,te2,ref,'y',seq);
    [RF,bHB]=sim_lcmrawbasis(np,sw,Bo,lb,'bHB',te1,te2,ref,'y',seq);
    [RF,bHG]=sim_lcmrawbasis(np,sw,Bo,lb,'bHG',te1,te2,ref,'y',seq);
    [RF,Glc]=sim_lcmrawbasis(np,sw,Bo,lb,'Glc',te1,te2,ref,'y',seq);
    [RF,NAAG]=sim_lcmrawbasis(np,sw,Bo,lb,'NAAG',te1,te2,ref,'y',seq);
    [RF,GPC]=sim_lcmrawbasis(np,sw,Bo,lb,'GPC',te1,te2,ref,'y',seq);
    [RF,PE]=sim_lcmrawbasis(np,sw,Bo,lb,'PE',te1,te2,ref,'y',seq);
    [RF,Ser]=sim_lcmrawbasis(np,sw,Bo,lb,'Ser',te1,te2,ref,'y',seq);
    [RF,EtOH]=sim_lcmrawbasis(np,sw,Bo,lb,'EtOH',te1,te2,ref,'y',seq);
  

%LEGEND:
%   'Ala'    = Alanine
%   'Asp'    = Aspartate
%   'PCh'    = PhosphoCholine
%   'Cr'     = Creatine
%   'PCr'    = PhosphoCreatine
%   'GABA'   = Gamma-aminobutyric acid
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
%   'EtOH'   = Ethanol
