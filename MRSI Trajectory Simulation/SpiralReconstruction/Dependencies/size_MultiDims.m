function size_vec = size_MultiDims(A,GetSizeAlong)
%
% template_1_0 Do nothing specific
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     This is output A
% -         B                           ...     This is output B
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('A','var'))
	error('\nInput needed.')
end
if(~exist('GetSizeAlong','var'))
	GetSizeAlong = 1:numel(size(A));
end


% 0.3 Definitions
    






%% 1. Compute A

size_vec = ones([1 numel(GetSizeAlong)]);
for CurDim = 1:numel(GetSizeAlong)
	size_vec(CurDim) = size(A,GetSizeAlong(CurDim));
end






%% 3. Postparations

% fclose(fid)






