% rf2d_getRFPeakPower.m
% Jamie Near, 2025.
%
% USAGE:
%   [rfPeakPower, rfPeakB1] = rf2d_getRFPeakPower(rf2d,params,TxAmpPwrRating,sysB1max)
%
% DESCRIPTION:
% DESCRIPTION:
% Calculate the instantaneous peak RF power of a given RF waveform given
% its duration, and some properties of the MR system transmitter: the
% transmitter amplifier power rating (TxAmpPwrRating, in kilowatts), and the 
% system peak B1 amplitude (sysB1max, in microtesla).  This function first 
% calculates the required peak B1 amplitude (rfPeakB1) for the desired RF 
% pulse and flip angle.  The function then determines what fraction of the 
% system peak B1 strength (f_B1) this represents.  Finally, the function 
% calculates the instantaneous peak power (rfPeakPower) by multiplying the system 
% transmitter amplifier power rating by the square of the peak B1 fraction, 
% f_B1.  For example, the body coil of a Siemens Prisma scanner is driven 
% with a 50kW amplifier, and delivers a peak B1 intensity of 26 microtesla.
% If a pulse requires a peak B1 intensity of 13 microtesla, the instantaneous
% peak power of the RF pulse will be 50 kW x (1/2)^2 = 12.5 kW.
%
%
% FORMULA:
%   rfPeakPower = (rfPeakB1 / sysB1max)^2 * TxAmpPwrRating
%
% INPUTS:
%   rf2d          - Structure with fields:
%                     • flip_angle  : flip angle [rad]
%                     • tp          : pulse duration [ms] 
%                   (In this implementation, hardware constants are
%                    hard-coded; see NOTES to externalize them.)
%
% OUTPUTS:
%   rfPeakPower   - Instantaneous peak transmit power [kW].
%   rfPeakB1      - Echo of input b1_peak [T].
%
% WARNINGS:
%   Prints a warning if rfPeakB1/sysB1max > 1 (required B1 exceeds system).

function [rfPeakPower,rfPeakB1]=rf2d_getRFPeakPower(rf2d)
tp = rf2d.tp;
rfPeakB1 = rf2d.b1_peak;
flipAngle = rf2d.flip_angle;
%Gyromagnetic ratio
gamma=42577000;  %[Hz/T].... assuming Proton

rfPeakB1=rfPeakB1*(rad2deg(flipAngle)/90);
%Find out the required peak B1 for the input RF pulse:
%Calculate the fraction of available system B1 amplitude:
f_B1=rfPeakB1/sysB1max;

%Warning if pulse exceeds B1 requirements:
if f_B1>1
    disp('WARNING: Pulse exceeds system B1 amplitude limit!!  Calculating anyway.');
end

%Calculate the instantaneous peak power of the RF pulse:
rfPeakPower = (f_B1^2)*TxAmpPwrRating;