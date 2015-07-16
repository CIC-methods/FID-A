function [par] = getparam(fpath)

% GETPARAM reads in parameter used in FID acquisition
%
%   PAR = GETPARAM(FPATH) reads the parameters in the FID directory
%   specified by FPATH and stores them in the structure PAR.
%
% see also FIDREAD, SETPARAM


%
% History
% 2003/03/14  Written by L Martyn Klassen
%

if nargin < 1
   error('GETPARAM requires one input argument.');
end

par = readprocpar(fpath);

% Calculate the number of complex pairs as a favor since it is an often required
par.nc = par.np/2;

% Some variables are sometimes set to zero when they really should be one
% Even VNMR does this for reasons it never fully explains. It basically treats
% zeros values as implied one values, so we just have to make that explicit so 
% that downstream user don't have to be continually checking for zero values
if par.ni < 1, par.ni = 1; end;
if par.nv < 1, par.nv = 1; end;
if par.nv2 < 1, par.nv2 = 1; end;
if par.nv3 < 1, par.nv3 = 1; end;
if par.ne < 1, par.ne = 1; end;
if par.ne1 < 1, par.ne1 = par.ne; end;


% Check that the slice info is consistent. People seem to like to change
% this without ensuring they stick to VNMR specification.
if length(par.pss) ~= par.ns & par.ns ~= 1
    warning('PROCPAR contains invalid slice information. Assuming uncompressed slices')
    par.seqcon(2) = 's';
    par.ns = 1;
end

% It makes no sense to have par.ns > 1 without compression of slice
% assume somehow the procpar got corrupted
if par.ns > 1
    par.seqcon(2) = 'c';
end

% Calculate some sequence specific information
if length(par.seqfil) >= 5
    if par.seqfil(1:5) == 'gedse'
        par.ne1 = par.ne - par.ne2;
        timespin = par.te/2 + [-par.esp*(par.ne2-1)/2:par.esp:par.esp*(par.ne2-1)/2];
        par.time = par.te/2 + [-1*timespin(par.ne1:-1:1) timespin];
    end
end

return
