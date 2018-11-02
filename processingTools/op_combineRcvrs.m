% op_combineRcvrs.m
% Jamie Near, McGill University 2017.
% 
% USAGE:
% [out,outw,out_presum,outw_presum,weights]=op_combineRcvrs(in,inw);
% 
% DESCRIPTION:
% Perform weighted coil recombination on both water suppressed and water 
% unsuppressed MRS data acquired with a reciever coil array.
% 
% INPUTS:
% in            = Water suppressed input spectrum in matlab structure format.
% inw           = Water unsuppressed input spectrum in matlab structure format.
%
% OUTPUTS:
% out           = Water suppressed output following combination of RF 
%                 channels.  
% outw          = Water unsuppressed output following combination of RF 
%                 channels.
% out_presum    = Water suppressed output with RF channels aligned but not 
%                 combined. 
% outw_presum   = Water unsuppressed output with RF channels aligned but not
%                 combined.
% weights       = Structure containing the coil weights and phases that were
%                 applied.

function [out,outw,out_presum,outw_presum,weights]=op_combineRcvrs(in,inw);

%first find the weights using the water unsuppressed data:
weights=op_getcoilcombos(inw,2,'h');

%now apply the weights to both the water unsuppressed and water suppressed
%data, but don't combine the averages:
out_presum=op_alignrcvrs(in,2,'h',weights);
outw_presum=op_alignrcvrs(inw,2,'h',weights);

%now apply the weights and combine the averages:
out=op_addrcvrs(in,2,'h',weights);
outw=op_addrcvrs(inw,2,'h',weights);

