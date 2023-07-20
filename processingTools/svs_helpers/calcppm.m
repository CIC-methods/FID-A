function ppm=calcppm(sw,np,Bo,gamma)
% Helper function to calculate ppm
% Peter Truong, Sunnybrook Research Institute, 2023
%
% INPUTS:
% sw    = spectral width
% np    = number of points
% Bo    = magnet strength (T)
% gamma = gyromagnetic ratio (MHz/T)
% OUTPUTS:
% ppm = the frequency range in ppm
% 

f=((-sw/2)+(sw/(2*np)):sw/np:(sw/2)-(sw/(2*np)));
ppm=f/(Bo*gamma);

% if 1H, shift by 4.65
if gamma==42.577
    ppm=ppm+4.65;
end
end