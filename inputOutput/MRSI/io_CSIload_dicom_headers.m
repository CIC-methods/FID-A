function dicomHeaders = MRSI_read_dicom_headers(folderName)
    arguments
        folderName (1, :) char {mustBeFolder}
    end
    dicomFolder = dir(folderName);
    % remove . and .. files from dicomFolder
    dicomFolder = dicomFolder(~ismember({dicomFolder.name}, {'.', '..'}));

    dicomPath = arrayfun(@(listing) fullfile(listing.folder, listing.name), dicomFolder, 'UniformOutput', false);

        
    dicomHeaders = spm_dicom_headers(dicomPath);
end