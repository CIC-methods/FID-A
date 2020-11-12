function mask_ShrunkOrGrown = MaskShrinkOrGrow(mask,NoOfVoxels,MaskGrow_flag,SliceBySlice_flag)
%
% mask_ShrunkOrGrown Shrink or Grow Mask
%
% This function was written by Bernhard Strasser.
%
%
% The function makes a mask bigger or smaller by a certain amount of voxels.
%
%
% mask_ShrunkOrGrown = MaskShrinkOrGrow(mask,NoOfVoxels,MaskGrow_flag,SliceBySlice_flag)
%
% Input: 
% -         mask                       ...     Input mask which should be shrunk or enlarged.
% -         NoOfVoxels                 ...     By how many voxels should the mask be shrunk or enlarged?
% -         MaskGrow_flag              ...     If true, the mask will be enlarged/grown. If false it will be shrunk.
% -         SliceBySlice_flag          ...     If true, each slice of the mask will be enlarged or shrunk independently. 
%                                              If false, the 3D-mask will be enlarged or grown as a whole.
% Output:
% -         mask_ShrunkOrGrown         ...     Output data.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
        



%% 0. Assign standard values

if(~exist('mask','var'))
    mask_ShrunkOrGrown = 0;
	fprintf('Input needed.')
	return
end

if(~exist('MaskGrow_flag','var'))
    MaskGrow_flag = true;
end
   
if(~exist('NoOfVoxels','var') || NoOfVoxels < 1)
    NoOfVoxels = 1;
end
if(~exist('SliceBySlice_flag','var'))
    SliceBySlice_flag = true;
end



%% Grow

% Very slow implementation but works
if(~MaskGrow_flag)
	mask = ~mask;
end
mask_ShrunkOrGrown = mask;	


if(SliceBySlice_flag)
	mask2 = EllipticalFilter(ones(size(mask(:,:,1))), [1 2], [1 1 1 NoOfVoxels],1);
else
	mask2 = EllipticalFilter(ones(size(mask)), [1 2 3], [1 1 1 NoOfVoxels],1);
end

for z = 1:size(mask,3)
	if(SliceBySlice_flag)
		mask4 = mask(:,:,z);
	end
	for x = 1:size(mask,1)
		for y = 1:size(mask,2)
			if(mask(x,y,z) == 1)
				continue;
			end			

			% create a second mask with radius NoOfVoxels
			ShiftBy = [x-(ceil(size(mask2,1)/2))-1, y-(ceil(size(mask2,1)/2))-1, z-(ceil(size(mask2,1)/2))-1];
			ShiftBy = ShiftBy(1:numel(size(mask2)));
			mask3 = NonCircShift(mask2,ShiftBy);
			if(~isempty(find(mask4(logical(mask3)),1)))
				mask_ShrunkOrGrown(x,y,z) = 1;
			end

		end
	end
end

	
if(~MaskGrow_flag)
	mask_ShrunkOrGrown = ~mask_ShrunkOrGrown;
end


   