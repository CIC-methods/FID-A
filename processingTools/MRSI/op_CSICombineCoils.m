%op_CSICombineCoils.m
%
% Combines the coils of MRSI data. First, coil's phases are adjusted to
% match the phase of the first coil. Then, sensetivity weights for coils are
% calculated by the formula w(i) = S(i)/sqrt(sum(S^2)) where S is intensity of the
% first point of the fids in the center of k space and i is the coil
% number. Coils are then summed by F = sum(W(i)*f(i)) where W(i) is coil's
% sensetivity weight at coil i, f(i) is the fids at coil i and F is the
% final summed signal.
% 
% USAGE:
% in = op_CSICombineCoils(in)
% 
% INPUT:    
% in            = input MRSI object for coils to be combined
% samplePoint(optional)   = samplePoint from fid to be used for phase
%               correction. Default to first point.
% phaseMap(optional)   = 3D matrix of phases for all coils and points. Dimensions are
%                       [num coils, num x, num y]
% weightMap(optional)   = 3D matrix of weights for all coils and points.
%                         Dimensions are the same as above.
%
% OUTPUT:
% in            = MRSI object with combined coils in fids and specs.

function [in, phases, weights] = op_CSICombineCoils(in,samplePoint, phaseMap, weightMap)
arguments
    in (1,1) struct
    samplePoint (1,1) double = 1
    phaseMap (:,:) double = []
    weightMap (:,:) double = []
end

            
    %some pre condition checks and setting default values
    
    if in.flags.addedrcvrs == 1
        error('coils already combined!')
    end
    
    if in.flags.spatialFT ~= 1 || in.flags.spectralFT ~= 0
        error('Please us op_CSIFourierTransform along the spatial dimension')
    end
    
    
    if isempty(phaseMap)
        %phase map arranged by coils, x coordinate, y coordinate
        phaseArray = zeros(in.sz(nonzeros([in.dims.coils, in.dims.x, in.dims.y, in.dims.averages])));
        specs = permute(in.specs, nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, in.dims.averages]));
        if(in.dims.averages == 0)
            ave_dims = 1;
        else
            ave_dims = in.sz(in.dims.averages);
        end
        for a = 1:ave_dims
            for x = 1:in.sz(in.dims.x)
                for y = 1:in.sz(in.dims.y)
                    for c = 1:in.sz(in.dims.coils)
                        %calculate the phase of each each coordinate for all the
                        %coils
                        phaseArray(c,x,y,a) = phase(specs(samplePoint,x,y,c,a));
                    end
                end
            end
        end
        
        phases = phaseArray; 
    else
        phases = phaseMap;
    end
    %permute dims
    permute_dims = nonzeros([in.dims.coils, in.dims.x, in.dims.y, in.dims.averages, in.dims.t]);
    
    %dims used to permute back once permuted
    permute_back(permute_dims) = 1:length(permute_dims);
    
    %permute dims to allow for phase multiplication
    specs = permute(in.specs, permute_dims);
    %mutliply phases
    specs = specs.*exp(-1i*phaseArray);
    %permute back to original dims
    in.specs = permute(specs, permute_back);
    
    if isempty(weightMap)
        %initalizing variables for the weighting functions
        weights = zeros(in.sz(nonzeros([in.dims.coils, in.dims.x, in.dims.y, in.dims.averages])));
        specs = permute(in.specs, nonzeros([in.dims.t, in.dims.coils, in.dims.x, in.dims.y, in.dims.averages]));
        %getting all the first points at each voxel
        if(in.dims.averages == 0)
            ave_dims = 1;
        else
            ave_dims = in.sz(in.dims.averages);
        end
        for a = 1:ave_dims
            for x = 1:in.sz(in.dims.x)
                for y = 1:in.sz(in.dims.y)
                    for c = 1:in.sz(in.dims.coils)
                        %calculate the phase of each each coordinate for all the
                        %coils
                        weights(c,x,y,a) = specs(samplePoint,c,x,y,a);
                    end
                end
            end
        end
        %Get the root sum squared value of coils
        weight_sum = squeeze(rssq(weights, 1));
        %permute for easy multiplication
        permute_dims = [2,3,4,1];
        back(permute_dims) = 1:length(permute_dims);
        weights = permute(weights, permute_dims);
        %element wise multiplication
        weights_norm = weights./weight_sum;
        %permute back
        weights_norm = permute(weights_norm, back);
        %getting the weights factor which is first fid
        
        weights = weights_norm;
    else
        weights = weightMap;
    end
    %adding weights to each coil
    permute_dims = nonzeros([in.dims.coils,in.dims.x,in.dims.y, in.dims.averages, in.dims.t]);
    permute_back(permute_dims) = 1:length(permute_dims);
    
    permuted_specs = permute(in.specs, permute_dims);
    combinedSpecs = permuted_specs.*weights;
    combinedSpecs = permute(combinedSpecs, permute_back);
    combinedSpecs = squeeze(sum(combinedSpecs, in.dims.coils));
    
    
    fn = fieldnames(in.dims);
    %updating MRSI parameters
    for names = 1:length(fn)
        if(in.dims.(fn{names}) > in.dims.coils)
            in.dims.(fn{names}) = in.dims.(fn{names}) - 1;
        end
    end
    in.dims.coils = 0;
    in.flags.addedrcvrs = 1;
    in.specs = combinedSpecs;
    in.sz = size(in.specs);

end
