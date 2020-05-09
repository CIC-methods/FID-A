function FieldExists_flag = isfield_recursive(InStruct,RecursiveFieldsToCheck)
%
% Recursively check if fields exist in structure
%
% This function was written by Bernhard Strasser.
%
%
% The function computes an exponential filter in Hertz
%
%
% FieldExists_flag = isfield_recursive(InStruct,RecursiveFieldsToCheck)
%
% Input: 
% -         InStruct                            ...    Structure which should be tested for fields
% -         RecursiveFieldsToCheck              ...    String of field-names which should be tested. The subfields are separated by '.'
%
% Output:
% -         FieldExists_flag                    ...    Bool saying if all subfields exist or not
%
%
% Example:
% Test.Operator.Input = 1; isfield_recursive(Test,'Operator.Input') % Result: 1
% isfield_recursive(Test,'Operator.In') % Result: 0
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy:

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% Define variables if not defined






%% 1. Check for 

FieldExists_flag = true;

RecursiveFieldsToCheck_cell = regexp(RecursiveFieldsToCheck,'\.','split');

CurStruct = InStruct;
for CurFieldNo = 1:numel(RecursiveFieldsToCheck_cell)
	CurField = RecursiveFieldsToCheck_cell{CurFieldNo};
	
	if(isfield(CurStruct,CurField))
		CurStruct = CurStruct.(CurField);
	else
		FieldExists_flag = false;
		break
	end
	
end


