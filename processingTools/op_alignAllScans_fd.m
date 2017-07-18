% op_alignAllScans_fd.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% [out,ph,frq]=op_alignAllScans_fd(in,fmin,fmax,tmax,ref,mode);
% 
% DESCRIPTION:
% Use spectral registration to align many separate scans.  In this version 
% only a limited frequncy range is used for fitting;
% 
% INPUTS:
% in        = cell array of inputs (spectra all to be registered)
% fmin      = Minimum frequency for spectral alignment (ppm).
% fmax      = Maximum frequency for spectral alignment (ppm).
% tmax      = Maximum time (s) in time domain to use for registration.
% ref       = Align to what? (optional)
%                  'f' - Align to the first input?  (default)
%                  'a' - Align to average of inputs. 
% mode      = (optional)'f' - Frequency align only
%                  'p' - Phase align only
%                  'fp or pf' - Frequency and phase align (default)
%
% OUTPUTS:
% out       = Cell array of multiple datasets after alignment.
% ph        = Vector of phase shifts (in degrees) used for alignment.
% frq       = Vector of frequency shifts (in Hz) used for alignment.


function [out,ph,frq]=op_alignAllScans_fd(in,fmin,fmax,tmax,ref,mode);

%Make sure input is a cell array with length greater than 2;
if ~iscell(in)  || length(in)<2
    error('ERROR:  The input must be a cell array containing two or more MRS datasets in FID-A structure format.  ABORTING!!');
end

%Provide default values if last two arguments are not given.
if nargin<6
    mode='fp';
    if nargin<5
        ref='f';
    end
end

%figure out what will be the reference spectrum and align spectra:
switch ref
    case 'f'
        in1=in{1};
        out{1}=in{1};
        ph(1)=0;
        frq(1)=0;
        for n=2:length(in)
            [out{n},ph(n),frq(n)]=op_alignScans_fd(in{n},in1,fmin,fmax,0.5);
        end
    case 'a'
        in1=in{1};
        for n=2:length(in)
            in1=op_addScans(in1,in{n});
        end
        in1=op_ampScale(in1,1/length(in));
        
        for n=1:length(in)
            [out{n},ph(n),frq(n)]=op_alignScans_fd(in{n},in1,fmin,fmax,0.5);
        end
        
    otherwise
        error('ERROR: Reference spectrum not recognized;');
end
         
