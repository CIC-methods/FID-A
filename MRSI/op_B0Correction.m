%op_B0Correction.m
%Jamie Near, McGill University 2019.
%Edits from Brenden Kadota, 2019.
%
% USAGE:
% [out, phaseMap, freqMap]=op_B0Correction(in, (optional) phaseMap);
% 
% DESCRIPTION:
% Corrects B0 shift by fitting all the curves to a water peak. Outputs the
% B0Map, phaseMap, and freqMap. Can accept B0 map to apply corrections from the map
% given.
% 
% INPUTS:
% in   = MRSI struct
% phaseMap = phaseMap matrix used to apply corrections
%
% OUTPUTS:
% B0Map        = 2D array of shift intensity at each specific voxel
% phaseMap     = 2D array of phase shift intensity at each specific voxel
% freqMap      = 2D array of frequency shift intensity at each specific voxel 
function [out, phaseMap, freqMap] = op_B0Correction(in, freqMap)
    
    %only correct if coils have be combined
    if(in.flags.addedrcvrs == 0)
        error("add coils first");
    end
      
    %find the coordinate of the highest intensity
    [x,y] = findMaxCoord(in);
    
    %create twix object from coordinate with highest intensity
    baseTwix = CSItoMRS(in, x, y);
    baseTwix.flags.averaged = 1;
    baseTwix.flags.isISIS = 0;
    
    %initalize phse map
    phaseMap = zeros([in.sz(in.dims.x), in.sz(in.dims.y)]);
    
    
    if(nargin == 2)
        %applying the freqMap map if added as an argument
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
               %create twix object from the coordinates
               alignTwix = CSItoMRS(in, x, y);
               %align twix object with op_freqshift to the equivalent
               %freqmap coordinate.
               alignTwix = op_freqshift(alignTwix, freqMap(x, y));
               %add the specs to MRSI object
               in.specs(:,x,y) = alignTwix.specs;
            end
        end
    else
        %aligning frequencies if no freqMap map provided
        freqMap = zeros([in.sz(in.dims.x), in.sz(in.dims.y)]);
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
                %get twix object at x and y coordinate
                alignTwix = CSItoMRS(in, x, y);
                alignTwix.flags.averaged = 1;
                alignTwix.flags.isISIS = 0;
                %aligns twix object to the maximum intensity twix object
                [alignTwix, phase, freq] = op_alignScans(alignTwix, baseTwix, in.t(end));
                %adds the aligned twix to the specs of MRSI object
                in.specs(:,x,y) = alignTwix.specs;
                %add phase and frequency to maps
                phaseMap(x,y) = phase;
                freqMap(x,y) = freq;
            end
        end
        
        %plot maps
        figure;
        subplot(2,1,1);
        imagesc(in.k_XCoordinates, in.k_YCoordinates, phaseMap)
        title("phaseMap")
        subplot(2,1,2);
        imagesc(in.k_XCoordinates, in.k_YCoordinates, freqMap)
        title("freqMap")

    end
    %return the modified MRSI object
    out = in;
           
end


%finds the (x,y) coordinate with the highest intensity spectrum.
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
