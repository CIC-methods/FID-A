function gamma = getGamma(argumentVariables)
    arguments 
        argumentVariables.overTwoPi (1, 1) logical = 1;
    end
    if(argumentVariables.overTwoPi)
        gamma = 42.577478518e6;
    else
        gamma = 267.52218744e6;
    end
end