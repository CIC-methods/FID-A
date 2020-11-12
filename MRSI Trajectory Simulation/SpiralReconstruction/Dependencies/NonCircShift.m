function OutArray = NonCircShift(InArray,ShiftArray)
%
% NonCircShift Shift Array like circshift but set wrapped elements to zero.
%
% This function was written by Bernhard Strasser, November 2014.
%
%
%
% OutArray = NonCircShift(InArray,ShiftArray)
%
% Input: 
% -         InArray                     ...     Obvious
% -         ShiftArray                 ...      Shift Matrix by this vector
%
% Output:
% -         OutArray                    ...     Obvious
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!




%% 0. Assign standard values

if(~exist('InArray','var'))
    OutArray = 0;
	fprintf('Input needed.')
	return
end

   
if(~exist('ShiftArray','var') || numel(ShiftArray) < 1 || numel(ShiftArray~=0) < 1)
    OutArray = InArray;
	fprintf('Nothing to do.')
	return
end

if(numel(ShiftArray) > numel(size(InArray)))
	fprintf('\nWarning: numel(ShiftArray) > numel(size(InArray)). Cut end of ShiftArray.')
	ShiftArray = ShiftArray(1:numel(size(InArray)));
elseif(numel(ShiftArray) < numel(size(InArray)))
	fprintf('\nWarning: numel(ShiftArray) < numel(size(InArray)). Append zeros.')
	ShiftArray = cat(  2, reshape(ShiftArray,[1 numel(ShiftArray)]) , zeros([1 numel(size(InArray))-numel(ShiftArray)])  );
end



%% Shift Array

OutArray = circshift(InArray,ShiftArray);



%% Set entries to zero which are not wanted


EvalString_temp = repmat(':,',[1 numel(size(ShiftArray))]);
EvalString_temp(end) = [];

for dim = 1:numel(size(ShiftArray))
	
	EvalString = EvalString_temp;
	
	if(ShiftArray(dim) == 0)
		continue;
	elseif(ShiftArray(dim) < 0)
		SetZeroLeftBorder = size(InArray,dim) + ShiftArray(dim) + 1;		% + ShiftArray(dim) subtracts, because its negative
		SetZeroRightBorder = size(InArray,dim);
	else
		SetZeroLeftBorder = 1;
		SetZeroRightBorder = ShiftArray(dim);		
	end
	
	EvalString = strcat(EvalString(1:(dim-1)*2),num2str(SetZeroLeftBorder), ':' ,num2str(SetZeroRightBorder), EvalString((dim-1)*2+2:end));
	eval([ 'OutArray(' EvalString ') = 0;' ]);
	
	
end
















