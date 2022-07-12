function checkSpectralFT(MRSIStruct) 
    if(~getFlags(MRSIStruct, 'spectralFT'))
        error('Please Fourier transform along the spectral dimension')
    end
end
