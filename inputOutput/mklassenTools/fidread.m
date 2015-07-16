% FIDREAD allows one to read in the Varian 4T FID data acquired with:
%
% [par, img, k, kraw] = fidread(fpath, nodc)
%
% par      - important parameters as read by GETPARAM
% img      - complex image data
% k        - k-space data
% kraw     - raw k-space data (2D matrix of one column per block)
%
% Input (optional)
% fpath    - name of FID directory to read
% nodc     - do not apply DC correction

function [par, img, k, kraw] = fidread(fpath,nodc);

% Get the FID directory to read
if nargin < 1
   fpath = fidpath;
else
   fpath = fidpath(fpath, 1);
end

% Turn on DC correction by default
if nargin < 2
   nodc = 0;
end

% Get MR parameters from procpar file and place in par structure
par = getparam(fpath);
[kraw, hdr] = readfid(fpath,par,~nodc);

[m,n] =size(kraw);
ntraces = hdr.ntraces;
clear hdr;

if par.ni > 1
   nimage = par.arraydim/par.ni;
   k = permute(reshape(kraw, [m/ntraces, ntraces, nimage, n/nimage]), [1 4 3 2]);
else
   k = reshape(kraw, [m/ntraces, ntraces, n]);
end

% Preform the Fourier Transform into complex image space
[m,n,o] = size(k);
img = k;
r1 = floor(m/2);
r2 = floor(n/2);
img(:,:,:) = fft2(k([r1+1:m 1:r1], [r2+1:n 1:r2],:));
r1 = ceil(m/2);
r2 = ceil(n/2);
img(:,:,:) = img([r1+1:m 1:r1], [r2+1:n 1:r2],:);
