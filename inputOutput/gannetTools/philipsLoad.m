% philipsLoad.m
% Georg Oeltzschner, Johns Hopkins University 2018.
%
% USAGE:
% data=philipsLoad(filename);
% 
% DESCRIPTION:
% Reads in philips MRS data (.spar and .sdat files) using code adapted from 
% PhilipsRead.m, provided as part of the Gannet software package by Richard 
% Edden (gabamrs.com).
% 
% INPUTS:
% filename   = filename of Philips .sdat file to be loaded.
%
% OUTPUTS:
% data       = FID array.
% header     = Structure containing header information from SPAR file.

function [data, header] = philipsLoad(filename)

% Get the .spar filename
sparname = [filename(1:(end-4)) 'spar'];

% Populate the header information from the SPAR file
% Look for regular expression separated by a colon
fid_spar = fopen(sparname);
while ~feof(fid_spar)
    tline = fgets(fid_spar); % get first line
    [tokens,matches] = regexp(tline,'([\w\[\].]*)\s*:\s*([\w\s.\"\\:\.]*)','tokens','match');
        % When a matching string is found, parse the results into a struct
        if length(tokens) == 1
            fieldname = regexprep(tokens{1}{1}, '\[|\]',''); % delete invalid characters
            
            % Convert numbers to doubles, leave strings & empty lines alone
            if ~isnan(str2double(tokens{1}{2}))
                value = str2double(tokens{1}{2});
            else
                value = tokens{1}{2};
            end

            header.(fieldname) = value;

        end
end
fclose(fid_spar);

% Now open the SDAT file and read in the data
fid_sdat    = fopen(filename,'r','ieee-le');
sdat_length = length(fread(fid_sdat));
fclose(fid_sdat);
fid_sdat    = fopen(filename,'r','ieee-le');
data        = freadVAXG(fid_sdat,sdat_length,'float32');
fclose(fid_sdat);
% Reshape and save
data        = reshape(data, [2 header.samples header.rows]);    
data        = squeeze(data(1,:,:)+1i*data(2,:,:));

end

