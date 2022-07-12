function nifti = io_CSIload_nifti(niftiFileName)
    arguments
        niftiFileName (1, :) char {mustBeFile}
    end
    
    matlabError = [];
    
    try 
        fileID = fopen(niftiFileName, 'r');
        size = fread(fileID, 1, '*int32')';
        
    catch matlabError
    end

    if(~isempty(matlabError))
        throw(matlabError)
    end

    fclose(fileID);

    if(size == 348)
        disp('Nifi-1 format detected. Loading...')
        nifti = io_CSIload_nifi1(niftiFileName);
    elseif(size == 540)
        disp('Nifti-2 format detected. Loading...')
        nifti = io_CSIload_nifti2(niftiFileName);
    end

end