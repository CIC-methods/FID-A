function in = op_combineCSICoils(in, samplePoint)
% input:    in = input twix object for coils to be combined
%           samplePoint = samplePoint from fid to be used for phase
%           correction. Default to first point.
%           samplePoint = point in in.fids to be used for phase correction
%           of the coils 
            
    %some pre condition checks and default variables
    if ~exist('samplePoint', 'var')
        samplePoint = 1;
    end
    if in.flags.addedrcvrs == 1
        error('coils already combined!')
    end
    
    %phase map arranged in coils, x coordinate, y coordinate
    phaseArray = zeros([in.sz(in.dims.coils), in.sz(in.dims.x), in.sz(in.dims.y)]);
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            for i = 1:in.sz(in.dims.coils)
                phaseArray(i,x,y) = phase(in.fids(samplePoint,i,x,y));
            end
        end
    end
    %[~,index] = max(phaseArray,[],'all', 'linear');
    %[coilIndex, ~, ~] = ind2sub(size(phaseArray), index);
    
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            phaseArray(:,x,y) = phaseArray(:,x,y) - phaseArray(1,x,y);           
        end
    end

    %multiplying the phase array to the fids and specs          
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            for i = 1:in.sz(in.dims.coils)
                if isfield(in, 'specs')
                    in.specs(:,i,x,y) = in.specs(:,i,x,y) .* exp(-1i * phaseArray(i,x,y));
                end
                in.fids(:,i,x,y) = in.fids(:,i,x,y) .* exp(-1i * -phaseArray(i,x,y));
            end
        end
    end
    
    %initalizing variables for the weighting functions 
    weights = zeros([1, in.sz(in.dims.coils)]);
    weightsFactor = zeros([1, in.sz(in.dims.coils)]);
    
    %getting all the first points of the fids from the center of k space
    weights(i) = in.fids(1,i,floor(end/2),floor(end/2));
    
    %getting the weights factor which is first fid 
    for i = 1:in.sz(in.dims.coils)
        weightsFactor(i) = weights(i)/sqrt(sum(weights.^2));
    end

    %adding weights to each coil 
    combinedSpecs = ones([in.sz(in.dims.t), in.sz(in.dims.x), in.sz(in.dims.y)]);
    combinedFids = ones([in.sz(in.dims.t), in.sz(in.dims.x), in.sz(in.dims.y)]);
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            if isfield(in, 'specs')
                combinedSpecs(:,x,y) = in.specs(:,:,x,y)*weightsFactor';
            end
            combinedFids(:,x,y) = in.fids(:,:,x,y)*weightsFactor';
        end
    end
    
    %updating twix parameters
    in.fids = combinedFids;
    in.dims.x = 2;
    in.dims.y = 3;
    in.dims.coils = 0;
    in.sz = size(in.fids);
    in.flags.addedrcvrs = 1;
    if isfield(in, 'specs')
        in.specs = combinedSpecs;
    end
end