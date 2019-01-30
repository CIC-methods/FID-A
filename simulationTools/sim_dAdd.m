%sim_dAdd.m
%Jamie Near, 2018.
%
% USAGE:
% d_out = sim_dAdd(d1,d2)
% 
% DESCRIPTION:
% Add together two density matrices.  This function is necessary becuase
% the density matrix is a cell array, so cannot be combined using simple
% addition. 
% 
% INPUTS:
% d1        = first input density matrix to be added.
% d2        = second input density matrix to be added.
% factor    = 1 for sum (default), -1 for diff
%
% OUTPUTS:
% d_out     = sum of d1 and d2.

function d_out = sim_dAdd(d1,d2,factor)

if nargin<3
    factor=1;
end

%If d1 is an empty cell, make d_out = d2; otherwise, add them:
if isempty(d1)
    d_out=d2;
else
    if size(d1) ~= size(d2)
        error('ERROR:  can only add density matrices of the same size!!  ABORTING!!');
    end
    
    d_out=cell(size(d1));
    for m=1:length(d1)
        d_out{m}=d1{m}+d2{m}.*factor;
    end
end



