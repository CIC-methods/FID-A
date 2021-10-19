%op_B0Correction.m
%Jamie Near, McGill University 2021.
%Edits from Brenden Kadota, 2021.
%
% USAGE:
% [out, phaseMap, freqMap]=op_B0Correction(in, (optional) phaseMap);
% 
% DESCRIPTION:
% Corrects for slight B0 shifts by aligning all voxels to the voxel with
% the highest water peak. (Possibly not correct output as highest peak may
% no have the correct ppm shifts)
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
    [x,y,a] = findMaxCoord(in);
    
    %create twix object from coordinate with highest intensity
    baseTwix = CSItoMRS(in, x, y, a);
    baseTwix.flags.averaged = 1;
    baseTwix.flags.isISIS = 0;
    
    %initalize phse map
    
    
    phaseMap = zeros([in.sz(in.dims.x), in.sz(in.dims.y)]);
    if(in.dims.averages ~= 0)
        phaseMap = repmat(phaseMap, [1, 1, in.sz(in.dims.averages)]);
    end
    
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
        if(in.dims.averages ~= 0)
            freqMap = repmat(freqMap, [1, 1, in.sz(in.dims.averages)]);
        end
        permute_dims = nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.averages]);
        permute_back(permute_dims) = 1:length(permute_dims);
        specs = permute(in.specs, permute_dims);
        if(in.dims.averages == 0)
            ave_dims = 1;
        else
            ave_dims = in.sz(in.dims.averages);
        end
        for a = 1:ave_dims
            for x = 1:in.sz(in.dims.x)
                for y = 1:in.sz(in.dims.y)
                    %get twix object at x and y coordinate
                    alignTwix = CSItoMRS(in, x, y, a);
                    alignTwix.flags.averaged = 1;
                    alignTwix.flags.isISIS = 0;
                    %aligns twix object to the maximum intensity twix object
                    [alignTwix, phase, freq] = op_alignScans(alignTwix, baseTwix, in.t(length(in.t)));
                    
                    %adds the aligned twix to the specs of MRSI object.
                    %Have to flip along the first dimension. (not really
                    %sure why but it turns out backwards if you don't)
                    specs(:,x,y,a) = flip(alignTwix.specs,1);
                    %add phase and frequency to maps
                    phaseMap(x,y,a) = phase;
                    freqMap(x,y,a) = freq;
                end
            end
        end
        in.specs = permute(specs, permute_back);
        
        %plot maps
        
        for i = 1:size(phaseMap,3)
            figure;
            subplot(2,1,1);
            imagesc(in.xCoordinates, in.yCoordinates, phaseMap(:,:,i)')
            title("phaseMap")
            subplot(2,1,2);
            imagesc(in.xCoordinates, in.yCoordinates, freqMap(:,:,i)')
            title("freqMap")
        end

    end
    %return the modified MRSI object
    out = in;
           
end


%finds the (x,y) coordinate with the highest intensity spectrum.
function [x,y,a] = findMaxCoord(in)
    [~, i] = max(in.specs,[], 'all', 'linear');
    out = cell(size(in.sz));
    [out{:}] = ind2sub(in.sz, i);
    x = out{in.dims.x};
    y = out{in.dims.y};
    if(in.dims.averages)
        a = out{in.dims.averages};
    else
        a = 1;
    end
    

end
