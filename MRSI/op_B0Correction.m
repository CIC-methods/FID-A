%op_B0Correction.m
%Jamie Near, McGill University 2019.
%Edits from Brenden Kadota, 2019.
%
% USAGE:
% [B0Map, phaseMap]=io_loadspec_twix(twix, (optional) phaseMap);
% 
% DESCRIPTION:
% Corrects B0 shift by fitting all the curves to a water peak. Outputs the
% B0Map and phaseMap. Can accept B0 map to apply corrections from the map
% given.
% 
% INPUTS:
% twix   = twix struct
%
% OUTPUTS:
% B0Map        = 2D array of shift intensity at specific voxel
% phaseMap     = 2D array of phase shift intensity at specific voxel
function [out, phaseMap, freqMap] = op_B0Correction(in, freqMap)
    if(in.flags.addedrcvrs == 0)
        error("add coils first");
    end
       
    [x,y] = findMaxCoord(in);
    baseTwix = CSItoMRS(in, x, y);
    baseTwix.flags.averaged = 1;
    baseTwix.flags.isISIS = 0;
    
    %applying the phase map if added as an argument
    if(nargin == 2)
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
               alignTwix = CSItoMRS(in, x, y);
               alignTwix = op_freqshift(alignTwix, freqMap(x, y));
               in.specs(:,x,y) = alignTwix.specs;
            end
        end
        
    else
        phaseMap = zeros([in.sz(in.dims.x), in.sz(in.dims.y)]);
        freqMap = zeros([in.sz(in.dims.x), in.sz(in.dims.y)]);
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
                alignTwix = CSItoMRS(in, x, y);
                alignTwix.flags.averaged = 1;
                alignTwix.flags.isISIS = 0;
                [alignTwix, phase, freq] = op_alignScans(alignTwix, baseTwix, in.t(end));
                in.specs(:,x,y) = alignTwix.specs;
                phaseMap(x,y) = phase;
                freqMap(x,y) = freq;
            end
        end
        figure;
        subplot(2,1,1);
        imagesc(in.k_XCoordinates, in.k_YCoordinates, phaseMap)
        title("phaseMap")
        subplot(2,1,2);
        imagesc(in.k_XCoordinates, in.k_YCoordinates, freqMap)
        title("freqMap")

    end
    out = in;
           
end

function [x,y] = findMaxCoord(in)
    maxNum = -realmax;
    x = -1;
    y = -1;
    for i = 1:in.sz(in.dims.x)
        for j = 1:in.sz(in.dims.y)
            tempMax = max(abs(in.specs(:,i,j)));
            if tempMax > maxNum
                maxNum = tempMax;
                x = i;
                y = j;
            end
        end
    end
end
