function [id] = getDicomPrivateTag(info,strGroupId,tagName)
% Search for a private DICOM tag.
%
% Inputs:
%
% info - structure returned by dicominfo().
% strGroupId - the DICOM group as a string of hexidecimal characters
%            - e.g. '0029'
%            - N.B. This value should be an ODD number according to the
%              DICOM standard.
% tagName - String containing the name of the private tag to be found.
%
% Returns a string with the leading digits of the key.

% Copyright Chris Rodgers, University of Oxford, 2011.
% $Id$

myFn = fieldnames(info);

myFnMatches = regexp(myFn,['^Private_' strGroupId '_([0-9]{2})xx_Creator$'],'tokens','once');

for myDx = 1:numel(myFn)
    if isempty(myFnMatches{myDx})
        continue
    end
    
%     fprintf('%s = ''%s''\n',myFn{myDx},info.(myFn{myDx}));
    if strcmp(tagName,info.(myFn{myDx}))
        id = myFnMatches{myDx}{1};
        return
    end
end

id = '';