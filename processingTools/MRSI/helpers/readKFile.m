% readKFile.m
% Brenden Kadota, Sunnybrook Research Institute 2021.
% Jamie Near, Sunnybrook Research Institute 2024.
%
% USAGE: 
% [kTable,kArray] = readKFile(kFile)
%
% DESCRIPTION:
%This function reads in a k-space coordinate file and returns it as a 
%MATLAB table, and as 2D array of k-space coordinates.  The k-space 
%coordinate file is a text file with comma-separated values (CSV format) 
%that describes the k-space coordinates of an acquired rapid MRSI 
%acquisition with arbitrary trajectory.  The first line of the kFile should 
%be formatted as follows: "TR,Kx,Ky,time".  The subsequent rows of the 
%kFile should specify, for each k-space sample in the trajectory, the TR 
%index number, the Kx coordinate in [mm^-1], the Ky coordinate in [mm^-1], 
%and the time of sampling relative to the start of the ADC.  For repeating
%k-space trajectories (i.e. most cases), the kFile need only specify the 
%coordinates for a single pass (to avoid redundancy).
%
%Note:  This function was previously called "readKSpaceFile.m",
%but the name was changed to "readKFile.m" for conciseness.  
%
%INPUTS:
% kFile      = The k-space coordinate file specifying all 
%                         k-space coordinates.
%
%OUTPUTS
% kTable:    = A Matlab table specifying the contents of the
%                         k-file.
% kArray:    = An 2D array of kx and ky coordinates.  

function [kTable, kArray] = readKFile(kFileName)
    
    %Read the CSV file into a MATLAB table:
    kTable = readtable(kFileName);

    %Make an array of Kx and Ky values:
    kArray = [kTable.Kx, kTable.Ky];
    
end