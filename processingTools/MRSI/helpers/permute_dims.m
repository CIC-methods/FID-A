% permute_dims.m
%
% Permute dimensions based on the input of dims. Non inputted diimensions are appeneded to the end
% INPUT
% in: MRSI strucutre
% dims: dimension array
% OUTPUT
% out: output MRSI structure with permuted dims
% prev: Previous permutation to use to permute back to original
function [out, prev_permute, prev_size] = permute_dims(in, dims)
    out = in;
    permute_dims = [dims, setdiff(1:ndims(in.fids), dims)];
    out.fids = permute(in.fids, permute_dims);
    
    prev_size = size(out.fids);
    dims_size = num2cell(prev_size(1:length(dims)));
    out.fids = reshape(out.fids, dims_size{:}, []);

    prev_permute = zeros(1, length(permute_dims));
    for i = 1:length(permute_dims)
       prev_permute(permute_dims(i)) = i;
    end

    if(isfield(out, 'specs'))
        out.specs = permute(in.specs, permute_dims);
        out.specs = reshape(out.specs, dims_size{:}, []);
    end
end
