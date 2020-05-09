function out_array = squeeze_single_dim(in_array,squeeze_dims)
%
% SearchHistory Searches the command history. 
%
% This function was written by Bernhard Strasser, July 2014.
%
%
% The function searches the command history (file prefdir/history.m) for SearchString.
%
%
% SearchHistory(SearchString,CaseInsensitive_flag)
%
% Input: 
% -     SearchString                   ...    String you want to search for.
% -     CaseInsensitive_flag           ...    If you want to search case insensitive, set this flag to true or 1.
%
% Output:
% None
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 




%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations

sz = size(in_array);                                 % size of the array


if(any(numel(sz) < squeeze_dims))
    TooHighDims = squeeze_dims(numel(sz) < squeeze_dims);
    fprintf('\nWarning in %s: You asked me to squeeze dimension(s) %s, but your input array has only %d dimensions.', mfilename,mat2str(TooHighDims),numel(sz))
    fprintf(' Do nothing on those dimensions.\n')
    squeeze_dims = setxor(squeeze_dims,TooHighDims);
end

if(any(sz(squeeze_dims) > 1))
    NonSingDims = squeeze_dims(sz(squeeze_dims) > 1);
    fprintf('\nWarning in %s: The dimension(s) %s is/are non-singleton, but has/ve size %s.',mfilename, mat2str(NonSingDims),mat2str(sz(NonSingDims)))
    fprintf(' Do nothing on this/these dimension(s).\n')
    squeeze_dims = setxor(squeeze_dims,NonSingDims);
% 	regy = cell([1 numel(sz)]);
% 	regy(:) = {':,'};
% 	regy{squeeze_dims} = '1,';
% 	regy{end} = regy{end}(1);
% 	
% 	% if it is the last index, its done.
% 	if(squeeze_dims == numel(sz))
% 		out_array = eval(['in_array(' [regy{:}] ');']);
% 		return;
% 	else
% 		in_array = eval(['in_array(' [regy{:}] ');']);		
% 	end
		
end
	
if(isempty(squeeze_dims))
    out_array = in_array;
    return
end
     

% 0.2 Declarations


% 0.3 Definitions
    




%% 1. Squeeze

new_dimensions = setxor(1:numel(sz), squeeze_dims);  % numel(size(..)) gives the number of dimensions of in_array; so 1:numel... gives [1 2 3 ... dimensionality(in_array)]; setxor: removes squeeze_dims out of this vec
new_size = sz(new_dimensions);                       % gives a vector with the same size as sz but where the sz(squeeze_dims) does not occur

out_array = reshape(in_array, new_size);


