%op_combinesubspecs.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=op_combinesubspecs(in,mode);
% 
% DESCRIPTION:
% Combine the subspectra in an acquisition either by addition or
% subtraction.
% 
% INPUTS:
% in        = input data in matlab structure format.
% mode      = -"diff" adds the subspectra together (this is counter intuitive,
%             but the reason is that many "difference editing" sequences use phase
%             cycling of the readout ADC to achieve "subtraction by addition".
%             -"summ" performs a subtraction of the subspectra.

function out=op_combinesubspecs(in,mode);

if in.flags.subtracted
    error('ERROR:  Subspectra have already been combined!  Aborting!');
end
if in.flags.isISIS
    error('ERROR:  MEGA-SPECIAL data must first be converted using makesubspecs.m!  Aborting!');
end

% if ~in.flags.freqcorrected
%     disp('WARNING:  Frequency correction has not yet been performed!');
% end
% if ~in.flags.phasecorrected
%     disp('WARNING:  Phase correction has not yet been performed!');
% end


if mode=='diff'
    %add the spectrum along the subSpecs dimension;
    fids=sum(in.fids,in.dims.subSpecs);
    fids=fids/2;  %divide by 2 so that this is an averaging operation;
elseif mode=='summ'
    %subtract the spectrum along the subSpecs dimension;
    fids=diff(in.fids,1,in.dims.subSpecs);
    fids=fids/2;  %divide by 2 so that this is an averaging operation;
end

fids=squeeze(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims.t=in.dims.t;
dims.coils=in.dims.coils;
dims.averages=in.dims.averages;
dims.subSpecs=0;

%re-calculate the sz variable
sz=size(fids);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;
out.subspecs=1;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.subtracted=1;
