% rf_getRFPeakPower.m
% Jamie Near, McGill University 2021.
%
% USAGE:
% [rfPeakPower,rfPeakB1]=rf_getRFPeakPower(RF,tp,TxAmpPwrRating,sysB1max,flipAngle);
% 
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
% INPUTS:
% RF                = RF pulse definition structure
% tp                = pulse duration in [ms]
% TxAmpPwrRating    = System transmitter amplifier power rating [kW]
% sysB1max          = System peak B1 amplitude [microtesla]
% flipAngle         = Desired RF pulse flip angle [Degrees] (Optional.
%                     Default = 90 deg for pulses of type 'exc'; 180 deg for
%                     pulses of type 'inv' and 'ref'.)
%
% OUTPUTS:
% rfPeakPower       = Instantaneous peak RF power of the RF pulse [kW]
% rfPeakB1          = Instantaneous peak B1 amplitude of the RF pulse [microtesla]


function [rfPeakPower,rfPeakB1]=rf_getRFPeakPower(RF,tp,TxAmpPwrRating,sysB1max,flipAngle);

%Gyromagnetic ratio
gamma=42577000;  %[Hz/T].... assuming Proton

%Find out the required peak B1 for the input RF pulse:
rfPeakB1=RF.tw1/gamma/(tp/1000)/1e-6;

%The above required B1 value is correct for pulses of type 'exc' with a
%flip angle of 90 degrees, and for pulses of type 'ref' or 'inv' with a
%flip angle of 180 degrees.  However, if a different flip angle was
%specified, we need to adjust the B1_required:
if nargin>4
    if strcmp(RF.type,'exc')
        rfPeakB1=rfPeakB1*flipAngle/90;
    elseif strcmp(RF.type,'inv') || strcmp(RF.type,'ref')
        rfPeakB1=rfPeakB1*flipAngle/180;
    end
end

%Calculate the fraction of available system B1 amplitude:
f_B1=rfPeakB1/sysB1max;

%Warning if pulse exceeds B1 requirements:
if f_B1>1
    disp('WARNING: Pulse exceeds system B1 amplitude limit!!  Calculating anyway.');
end

%Calculate the instantaneous peak power of the RF pulse:
rfPeakPower = (f_B1^2)*TxAmpPwrRating;



