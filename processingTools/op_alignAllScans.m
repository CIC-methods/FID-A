% op_alignAllScans.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% [out,ph,frq]=op_alignAllScans(in,tmax,ref,mode);
% 
% DESCRIPTION:
% Use spectral registration to align many separate scans;
% 
% INPUTS:
% in        = cell array of inputs (spectra all to be registered)
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

function [out,ph,frq]=op_alignAllScans(in,tmax,ref,mode);

%Make sure input is a cell array with length greater than 2;
if ~iscell(in)  || length(in)<2
    error('ERROR:  The input must be a cell array containing two or more MRS datasets in FID-A structure format.  ABORTING!!');
end

%Provide default values if last two arguments are not given.
if nargin<4
    mode='fp';
    if nargin<3
        ref='f';
    end
end

%figure out what will be the reference spectrum and align spectra:
switch ref
    case 'f'
        in1=in{1};
        ampRef=max(abs(in1.fids));
        out{1}=in{1};
        ph(1)=0;
        frq(1)=0;
        
        for n=2:length(in)
            amp=max(abs(in{n}.fids));
            dummyIn=op_ampScale(in{n},ampRef/amp);
            [dummyOut,ph(n),frq(n)]=op_alignScans(dummyIn,in1,0.5);
            out{n}=op_freqshift(op_addphase(in{n},ph(n)),frq(n));
        end
    case 'a'
        in1=in{1};
        for n=2:length(in)
            in1=op_addScans(in1,in{n});
        end
        in1=op_ampScale(in1,1/length(in));
        ampRef=max(abs(in1.fids));
        
        for n=1:length(in)
            amp=max(abs(in{n}.fids));
            dummyIn=op_ampScale(in{n},ampRef/amp);
            [dummyOut,ph(n),frq(n)]=op_alignScans(dummyIn,in1,0.5);
            out{n}=op_freqshift(op_addphase(in{n},ph(n)),frq(n));
        end
        
    otherwise
        error('ERROR: Reference spectrum not recognized;');
end
         
