function OutArray = myrepmat(OutArray,DesiredSize,DimensionCorrespondence)
%
% myrepmat Replicate matrix to desired size
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The normal MATLAB repmat function cannot add a dimension at the beginning of the array. E.g., if you have an array, size(array) = [64 64] and
% you want it to replicate it to size [32 64 64] there is no easy way to do that. If you do repmat(array, [32 1 1]) you get a array of size
% [64*32 64] = [2048 64]. This function gets you the desired [32 64 64] result. So it is a generalized repmat function, the result will always
% be your desired size (except if DesiredSize(i) < InArray(DimensionCorrespondance(1,i), because that wouldnt be replicating the matrix but
% cutting it down)
%
%
% OutArray = myrepmat(InArray,DesiredSize,DimensionCorrespondence)
%
% Input: 
% -  InArray                   ...    The input array that should be replicated. For memory reasons: InArray = OutArray;
% -  DesiredSize               ...    The desired size of the replicated array. size(InArray) must occur somewhere in this var.
% -  DimensionCorrespondence   ...    Matrix with numel(DesiredSize) elements, which tells the function which dimension of the
%                                     InArray should correspond to which index of the OutArray. zeros mean
%                                     that the dimension, given by the place in DimensionCorrespondence, should be created totally new.
%                                     If you dont replicate the array in a dimension that already exists (e.g. from
%                                     [64 64 8] to [64 64 16] replicates in an already existing dim, but
%                                     [64 64 8] to [32 64 64 8 2048] or [32 64 2 64 8] does not)
%                                     you can omit DimensionCorrespondence.
%                                     Example: 
%                                     size(InArray) = [64 64 16]; DesiredSize = [32 64 64 32 2048]; 
%                                     DimensionCorrespondence = [0 1 2 3 0];
%                                     That means: DesiredSize(1) = 32 should be repmatted from the scratch (so just take
%                                     the InArray, replicate it 32 times, and bring the new dimension to the right position).
%                                     The first index of size(InArray) should be replicated to the second index of the DesiredSize
%                                     (since they are the same, nothing has to be done). etc. etc.
%
% Output:
% -  OutArray                  ...    The replicated output array
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 




%% 0. Preparations, Definitions



% 0.1 Preparations

if(numel(size(OutArray)) == numel(DesiredSize) && sum(size(OutArray) ~= DesiredSize) == 0 )
	return;
end

if(~exist('DimensionCorrespondence', 'var'))
    OutArray = squeeze(OutArray);
end

Size_InArray = size(OutArray);
if(Size_InArray(end) == 1)
	Size_InArray(end) = [];
end
if(numel(Size_InArray) == 2 && Size_InArray(1) == 1)
	Size_InArray(1) = [];
end

if(~exist('DesiredSize', 'var'))
    return
end



if(~exist('DimensionCorrespondence', 'var'))
	
	DimensionCorrespondence = zeros([1 numel(DesiredSize)]);
	Size_InArray_dum = Size_InArray;
	% Create DimensionCorrespondence by checking index by index if Size_InArray is somehow distributed in DesiredSize
	for curoutind=1:numel(DesiredSize)
		
		% Check if DesiredSize(curoutind) is in Size_InArray_dum
		IndexMatch = Size_InArray_dum == DesiredSize(curoutind);
		if(sum(IndexMatch) > 0)
			IndexOfInArrayWhichMatched = find(IndexMatch);
			DimensionCorrespondence(curoutind) = IndexOfInArrayWhichMatched(1);
			Size_InArray_dum(IndexOfInArrayWhichMatched(1)) = 0;
		end
	end
	
	% If not all of the indices of Size_InArray were found in DesiredSize, the automatic calculation of DimensionCorespondence doesnt work. Display error.  
	if(sum(Size_InArray_dum) > 0)
       display([ char(10) 'Error: All elements of size(InArray) must occur in DesiredSize.' char(10) 'Automatic calculation of DimensionCorrespondence failed. Manual input needed.' ])
       return
	end
	
	clear Size_InArray_dum
	
end



%% 1. Exchange dims if DimensionCorrespondence is not ordered

DimensionCor_NoZeros = DimensionCorrespondence(DimensionCorrespondence>0);
if(numel(DimensionCor_NoZeros) > 1)
	OutArray = permute(OutArray, DimensionCor_NoZeros);
	Size_InArray = size(OutArray);
end



%% 2. Replicate


ReshapeTo = DesiredSize; ReshapeTo(DimensionCorrespondence == 0) = 1;
RepmatTo = DesiredSize; RepmatTo(DimensionCorrespondence ~= 0) = 1;

OutArray = reshape(OutArray,ReshapeTo);

if( sum(RepmatTo ~= 1) > 0 )
	OutArray = repmat(OutArray,RepmatTo);
end



