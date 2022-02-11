

function rotationDirection = getRotationDirection(signal, startIndex, endIndex)
    signalPhase = phase(signal(startIndex:endIndex));
    if(signalPhase(end) - signalPhase(1) > 0)
        rotationDirection = 1;
    else
        rotationDirection = -1;
    end
end