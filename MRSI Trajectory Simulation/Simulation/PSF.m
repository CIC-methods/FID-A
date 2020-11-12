% PSF.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% PSF(inTrajectory)
% 
% DESCRIPTION:
% This function takes in a trajectory strucutre and displays a graph of the
% point spread function. Furhter, it returns values for full width half
% maximum and signal leakage. The singal leakage is the sum of absolute
% value of the signal outside of the main signal peak.
%
% INPUTS:
% TrajectoryStruct = FID-A Trajectory structure from CreateTrajectory.m
%
% OUTPUTS:
% lineWidth   = the full width half maximum
% leakage     = leakage of point spread function signal outside of the main
%               peak



function [lineWidth,leakage] = PSF(TrajectoryStruct)


    [x,y] = meshgrid(-0.05:.001:0.05, -0.05:.001:0.05);
    pointSpreadCoordinates(:,:) = [x(:), y(:)];
    
    sftOperator = sft2_Operator(TrajectoryStruct.k_space_trajectory, pointSpreadCoordinates, 0);
    kSpace = ones(1,size(TrajectoryStruct.k_space_trajectory,1));
    pointSpreadFunction = sftOperator*kSpace(:);
    pointSpreadFunction = reshape(pointSpreadFunction, [size(x,1), size(x,2)]);
    peak = max(pointSpreadFunction, [], 'all');
    pointSpreadFunction = pointSpreadFunction./peak;

    figure;
    subplot(2,2,1), plot(inTrajeTrajectoryStruct.k_space_trajectoryctory(:,1), TrajectoryStruct.k_space_trajectory(:,2)), title('kSpace trajectory');
    
    subplot(2,2,2), mesh(-0.05:.001:0.05, -0.05:.001:0.05, ...
    real(pointSpreadFunction)), axis square, title('Point Spread Function')
    
    pointSpread2D = real(pointSpreadFunction(:,(floor(end/2)+1)));
    subplot(2,2,3), plot(-0.05:.001:0.05, pointSpread2D), title('2D point spread through the middle');
    
    lineWidth = fwhm(-0.05:.001:0.05, pointSpread2D);
    
    
    middle = find(pointSpread2D == max(pointSpread2D));
    startOfLeakage = find(pointSpread2D(middle:end) < 0, 1);
    leakage = sum(sqrt(pointSpread2D(startOfLeakage:end).^2));
    
    
end