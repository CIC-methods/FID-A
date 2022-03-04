function MRSIStruct = setFlags(MRSIStruct, flagName, value)
    flags = MRSIStruct.flags;
    switch lower(flagName)
        case 'writtentostruct'
            flags.writtentostruct = value; 
        case 'gotparams'
            flags.gotparams = value; 
        case 'leftshifted'
            flags.leftshifted = value; 
        case 'filtered'
            flags.filtered = value; 
        case 'zeropadded'
            flags.zeropadded = value; 
        case 'freqcorrected'
            flags.frefqcorrected = value; 
        case 'phasecorrected'
            flags.phasecorrected = value; 
        case 'averaged'
            flags.averaged = value; 
        case 'addedrcvrs'
            flags.addedrcvrs = value; 
        case 'subtracted'
            flags.subtracted = value; 
        case 'writtentotext'
            flags.writtentotext = value; 
        case 'downsampled'
            flags.downsampled = value; 
        case 'spatialft'
            flags.spatialFT = value; 
        case 'spectralft'
            flags.spectralFT = value; 
        case 'coilcombined'
            flags.coilCombined = value; 
        case 'isfoursteps'
            flags.isfoursteps = value; 
        case 'iscartesian'
            flags.iscartesian = value; 
        case 'apodized'
            flags.apodized = value;
        otherwise
            error('could not find the flag %s', flagName)
    end
    MRSIStruct.flags = flags;
end
