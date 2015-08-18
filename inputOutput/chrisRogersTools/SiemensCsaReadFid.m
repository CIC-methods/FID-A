function [signal, info] = SiemensCsaReadFid(info, bDeleteFromInfo, conjMode, debugMode)
% Read FIDs from a Siemens MR spectroscopy data set
%
% [fids] = SiemensCsaReadFid(info, bDeleteFromInfo, conjMode)
%
% info - DICOM headers read in with Matlab dicominfo() function
% bDeleteFromInfo - if true, remove private field from info to save RAM
% conjMode - if 'conj' then return the conjugate of the stored spectrum
% debugMode - if 'debug' then output extra logging information
%
% For improved performance, decode the CSA header with SiemensCsaParse
% prior to calling this function.

% Copyright: Chris Rodgers (University of Oxford), 2008-11.
% All rights reserved.

% $Id: SiemensCsaReadFid.m 6154 2013-02-11 16:08:19Z crodgers $

if nargin < 4 || (islogical(debugMode) && ~debugMode) || (ischar(debugMode) && ~strcmp(debugMode,'debug'))
    debugMode = false;
else
    debugMode = true;
end

if nargin < 3 || (islogical(conjMode) && ~conjMode) || (ischar(conjMode) && ~strcmp(conjMode,'conj'))
    conjMode = false;
else
    conjMode = true;
end

% Read CSA block if it is not supplied...
if ~isfield(info,'csa')
    info = SiemensCsaParse(info);
end

% Interpolated data dimensions
ipolDataPointColumns=info.csa.DataPointColumns;
ipolPhaseRows=info.csa.Rows;
ipolPhaseColumns=info.csa.Columns;
ipolNumberOfFrames=info.csa.NumberOfFrames;

if debugMode
fprintf('Interpolated dimensions from the DICOM file (DataPointColumns = %d, PhaseColumns = %d, PhaseRows = %d, NumberOfFrames = %d).\n',...
    ipolDataPointColumns,ipolPhaseColumns,ipolPhaseRows,ipolNumberOfFrames)
end

% Raw (acquired) data dimensions
rawDataColumns = info.csa.SpectroscopyAcquisitionDataColumns;
rawPhaseRows = info.csa.SpectroscopyAcquisitionPhaseRows;
rawPhaseColumns = info.csa.SpectroscopyAcquisitionPhaseColumns;
rawOutofplanePhaseSteps = info.csa.SpectroscopyAcquisitionOutofplanePhaseSteps;

if debugMode
fprintf('Raw dimensions from the DICOM file (DataColumns = %d, PhaseColumns = %d, PhaseRows = %d, OutofplanePhaseSteps = %d).\n',...
    rawDataColumns, rawPhaseColumns, rawPhaseRows, rawOutofplanePhaseSteps)
end

% Check whether these have been interpolated?
if ipolDataPointColumns ~= rawDataColumns || ...
    ipolPhaseRows ~= rawPhaseRows || ...
    ipolPhaseColumns ~= rawPhaseColumns || ...
    ipolNumberOfFrames ~= rawOutofplanePhaseSteps
    warning('SiemensCsaReadFid:Interpolated','These spectra have been spatially interpolated.')
end

%% Read the signal in as a complex FID
% Find appropriate private tag
specTag = ['Private_7fe1_' getDicomPrivateTag(info,'7fe1','SIEMENS CSA NON-IMAGE') '10'];

tmp = typecast(info.(specTag),'single');

if ~conjMode
    % Return in an array DATA POINT, VOXEL COLUMN, VOXEL ROW, OUT OF PLANE (Z)
    signal = reshape([tmp(1:2:end)+1i*tmp(2:2:end)],ipolDataPointColumns,ipolPhaseColumns,ipolPhaseRows,ipolNumberOfFrames);
else % conjMode
    % Return in an array DATA POINT, VOXEL COLUMN, VOXEL ROW, OUT OF PLANE (Z)
    signal = reshape([tmp(1:2:end)-1i*tmp(2:2:end)],ipolDataPointColumns,ipolPhaseColumns,ipolPhaseRows,ipolNumberOfFrames);
end

% Economise on memory by wiping the FID raw data from the DICOM struct if
% requested and if the wiped value can actually be assigned somewhere.
if bDeleteFromInfo && nargout >= 2
    info.(specTag) = [];
end
