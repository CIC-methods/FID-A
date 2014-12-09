%rf_readRF.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf,info]=rf_readRF(filename)
% 
% DESCRIPTION:
% Read a Varian/Agilent .RF file into matlab.  The resulting RF matrix will have 
% 3 columns specifying magnitude, phase, and duration.  This function simply calls
% Martyn Klassen's readrfvnmr.m function. 
% 
% INPUTS:
% filename   = filename of the .RF file to read in.  

function [rf]=rf_readRF(filename)

rf=readrfvnmr(filename);
   