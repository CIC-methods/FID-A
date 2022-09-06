function flagValue = getFlags(MRSIStruct, flagName)
    flags = MRSIStruct.flags;
    switch lower(flagName)
        case 'writtentostruct'
            flagValue = flags.writtentostruct;
        case 'gotparams'
            flagValue = flags.gotparams;
        case 'leftshifted'
            flagValue = flags.leftshifted;
        case 'filtered'
            flagValue = flags.filtered;
        case 'zeropadded'
            flagValue = flags.zeropadded;
        case 'frefqcorrected'
            flagValue = flags.freqcorrected;
        case 'phasecorrected'
            flagValue = flags.phasecorrected;
        case 'averaged'
            flagValue = flags.averaged;
        case 'addedrcvrs'
            flagValue = flags.addedrcvrs;
        case 'subtracted'
            flagValue = flags.subtracted;
        case 'writtentotext'
            flagValue = flags.writtentotext;
        case 'downsampled'
            flagValue = flags.downsampled;
        case 'spatialft'
            flagValue = flags.spatialFT;
        case 'spectralft'
            flagValue = flags.spectralFT;
        case 'coilcombined'
            flagValue = flags.coilCombined;
        case 'isfoursteps'
            flagValue = flags.isFourSteps;
        case 'iscartesian'
            flagValue = flags.isCartesian;
        otherwise
            error('could not find the flag %s', flagName)
    end
end