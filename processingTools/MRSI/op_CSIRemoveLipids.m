% op_CSIRemoveLipids.m
%
% Removes lipids from CSI data using L2 regularization. This minimizes the
% equation norm(x - x_0, 2) + beta * norm(W'x, 2)
%
% INPUT:
% MRSIStruct        = MRSI structure used in FID-A
% lipidComponents   = number of lipid spectra in the lipid basis
% lineWidthRange    = range of linewidth used for building lipid basis
% ppmRange          = ppm range used in building lipid basis
% beta              = regularization term
% plotBasis         = plot basis spectra
%                                                                                                                                                                                                                                                                          
%
% OUTPUT:
% MRSIStruct        = MRSI structure with lipids removed

function [MRSIStruct] = op_CSIRemoveLipids(MRSIStruct, basisArguments, plottingArguments) 
    arguments
        MRSIStruct (1, 1) struct
        basisArguments.lipidComponenets (1, 1) double = 1000
        basisArguments.lineWidthRange (1, 2) double = [1 80]
        basisArguments.lipidPPMRange (1, 2) double = [0.3 1.9000]
        basisArguments.beta (1, 1) double = 1
        plottingArguments.plotBasis (1, 1) logical = false
    end
    % extract arguments from name value pairs
    lipidComponents = basisArguments.lipidComponenets;
    lineWidthRange = basisArguments.lineWidthRange;
    lipidPPMRange = basisArguments.lipidPPMRange;
    beta = basisArguments.beta;

    spectraSize = getSizeFromDimensions(MRSIStruct, {'t'});

    % calculate lipid basis used for L2 regularization
    lipidBasis = createLipipBasis(MRSIStruct, lipidComponents, lineWidthRange, lipidPPMRange);
    if(plottingArguments.plotBasis)
        figure
        plot(MRSIStruct.ppm, flip(real(lipidBasis),1));
    end
    % calculate solution from the basis
    L2Solution = inv(eye(spectraSize) + beta * (lipidBasis * lipidBasis'));
    
    [MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t'});
    data = getData(MRSIStruct);
    
    data = L2Solution * data;
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);
end

% calculate lipid basis
function lipidBasis = createLipipBasis(MRSIStruct, lipidComponents, lineWidthRange, lipidPPMRange)
    spectralWidth = getSpectralWidth(MRSIStruct);

    spectralPoints = getSizeFromDimensions(MRSIStruct, {'t'});
    fidBasis = zeros(spectralPoints, lipidComponents);
    
    lipidStructure = load('Lip.mat', 'sysLip');
    lipidStructure = lipidStructure.sysLip;
    
    for iSpectra = 1:lipidComponents
        
        fidBasis(:, iSpectra) = getRandomLipidFids(spectralPoints, spectralWidth, ...
                                                  lineWidthRange, lipidPPMRange, ...
                                                  lipidStructure);
    end
    lipidBasis = fftshift(fft(fidBasis, [], 1), 1);
end

% calculate the single lipid spectra for the basis
function lipidFids = getRandomLipidFids(spectralPoints, spectralWidth, lineWidth, ...
                                        lipidPPMRange, lipidSystem)
    [randomLineWidth, randomPPM] = getRandomLineWidthandPPM(lineWidth, lipidPPMRange);
    lipidSystem.shifts = randomPPM;

    simulatedSignal = sim_onepulse(spectralPoints, spectralWidth, 3, randomLineWidth, lipidSystem);
    simulatedSignal = op_complexConj(simulatedSignal);
    simulatedSignal = addRandomPhase(simulatedSignal);
    %simulatedSignal = scaleSpectra(simulatedSignal, randomPPM, lipidPPMRange);
    lipidFids = simulatedSignal.fids;
end


% pick a random number from lower bounds and upper bounds
function randomNumber = randomNumberInRange(lowerBounds, upperBounds)
    difference = upperBounds - lowerBounds;
    randomNumber = lowerBounds + difference * rand(1);
end
                
function simulatedSignal = addRandomPhase(simulatedSignal)
    sepctralPhase = randomNumberInRange(-180, 180);
    simulatedSignal = op_addphase(simulatedSignal, sepctralPhase, 0, 4.65, 1);
end

% scale spectra based on a normal distribution. Signal near the center of lipid
% range will be scaled high and signal near the edges scaled lower.
function simulatedSignal = scaleSpectra(simulatedSignal, lipidPPM, lipidRange)
    normalProbability = normpdf(lipidPPM, mean(lipidRange), diff(lipidRange)/4);
    simulatedSignal = op_ampScale(simulatedSignal, normalProbability);
    simulatedSignal = op_ampScale(simulatedSignal, 10);

end

function [randomLineWidth, randomPPM] = getRandomLineWidthandPPM(lineWidth, lipidPPMRange)
    randomLineWidth = randomNumberInRange(lineWidth(1), lineWidth(2));
    randomPPM = randomNumberInRange(lipidPPMRange(1), lipidPPMRange(2));
end