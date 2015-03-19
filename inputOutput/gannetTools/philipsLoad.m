% philipsLoad.m
% gannetmrs.blogspot.com
%
% USAGE:
% data=philipsLoad(filename);
% 
% DESCRIPTION:
% Reads in philips MRS data (.spar and .sdat files) using code adapted from 
% PhilipsRead.m, provided as part of the Gannet software package by Richard 
% Edden (gabamrs.blogspot.com).
% 
% INPUTS:
% filename   = filename of Philips .sdat file to be loaded.

function data = philipsLoad(filename)

sparname = [filename(1:(end-4)) 'spar']
sparheader = textread(sparname, '%s');
sparidx=find(ismember(sparheader, 'samples')==1);
da_xres = str2num(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'rows')==1);
da_yres = str2num(sparheader{sparidx+2});

data = SDATreadMEGA(filename, da_xres, da_yres);


end

