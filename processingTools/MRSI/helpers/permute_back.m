% Used to permute back from pervious permutation from permute_dims
% INPUT
% in: MRSI structure that was permuted
% prev: prev array from premute_dims
% 
% OUTPUT
% out: permuted MRSI structure
function out = permute_back(in, prev_permute, prev_size)
    out = in;
    out.fids = reshape(out.fids, prev_size);
    out.fids = permute(in.fids, prev_permute);
    if(isfield(out, 'specs'))
        out.specs = reshape(out.specs, prev_size);
        out.specs = permute(in.specs, prev_permute);
    end

end
