%io_readRF.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf]=io_readRF(filename)
% 
% DESCRIPTION:
% Read a Varian/Agilent .RF file into matlab.  The resulting RF matrix will have 
% 3 columns specifying magnitude, phase, and duration.  This function simply calls
% Martyn Klassen's readrfvnmr.m function. 
% 
% INPUTS:
% filename   = filename of the .RF file to read in. 
%
% OUTPUTS:
% rf        = Input rf pulse waveform saved as a matlab array with 3
%               columns (magnitude, phase, duration).

function [rf]=io_readRF(filename)

rf=readrfvnmr(filename);
   