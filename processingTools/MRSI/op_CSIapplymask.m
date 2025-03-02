% op_CSIapplymask.m
% Emile Kadalie, Sunnybrook 2025.
%
% USAGE:
% out=op_CSIB0Correction_v2(in);
%
% DESCRIPTION:
% Applies precomputed brain mask to MRSI data to spatially remove
% surrounding lipids.
%
% INPUTS:
% in   = MRSI struct
% 
% OUTPUTS:
% out          = MRSI struct with masked data
function MRSIStruct = op_CSIapplymask(MRSIStruct)
    if isfield(MRSIStruct,'mask')
        MRSIStruct.data=MRSIStruct.data.*permute(repmat(MRSIStruct.mask.brainmasks,1,1,MRSIStruct.sz(1)),[3 1 2]);
    else
        error("No mask found. Please segment");
    end
end