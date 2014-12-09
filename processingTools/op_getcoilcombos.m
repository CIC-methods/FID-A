%op_getcoilcombos.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% coilcombos=op_getcoilcombos(file_or_struct,point);
% 
% DESCRIPTION:
% This funciton finds the relative coil phases and amplitudes.  Coil phases
% are found by determining the phase and amplitude of the pointth point in 
% the time domain.
% 
% INPUTS:
% file_or_struct     = this function will accept either a string filename or
%                     the name of a structure.  If the input is a string, 
%                     the program will read in the data corresponding to 
%                     that filename.  If the input is a structure, it will
%                     operate on that structure.
% point              = The index of the datapoint in the fid that is used
%                     for determination of signal intensity and phase.

function coilcombos=op_getcoilcombos(file_or_struct,point,mode);


if isstr(file_or_struct)
    in=op_loadspec_twix(file_or_struct);
else
    in=file_or_struct;
end

if nargin<3
    mode='w';
end

if in.flags.addedrcvrs
    error('ERROR:  must provide data prior to coil combination!!  ABORTING!!');
end

coilcombos.ph=zeros(in.sz(in.dims.coils),1);
coilcombos.sig=zeros(in.sz(in.dims.coils),1);

for n=1:in.sz(in.dims.coils);
    coilcombos.ph(n)=phase(in.fids(point,n,1,1));
    switch mode
        case 'w'
            coilcombos.sig(n)=abs(in.fids(point,n,1,1));
        case 'h'
            S=max(abs(in.fids(:,n,1,1)));
            N=std(in.fids(end-100:end,n,1,1));
            coilcombos.sig(n)=(S/(N.^2));
    end
end

