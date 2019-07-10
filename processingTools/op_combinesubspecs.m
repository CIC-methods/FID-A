% op_combinesubspecs.m
% Jamie Near, McGill University 2014.
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
%
% OUTPUTS:
% out       = Output following combination of subspectra.  

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
    fids=fids/in.sz(in.dims.subSpecs); %divide by number of subspecs so that this is an averaging operation;
elseif mode=='summ'
    %subtract the spectrum along the subSpecs dimension;
    fids=diff(in.fids,1,in.dims.subSpecs);
    fids=fids/in.sz(in.dims.subSpecs); %divide by nymber of subspecs so that this is an averaging operation;
end

fids=squeeze(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
if in.dims.t>in.dims.subSpecs
    dims.t=in.dims.t-1;
else
    dims.t=in.dims.t;
end
if in.dims.coils>in.dims.subSpecs
    dims.coils=in.dims.coils-1;
else
    dims.coils=in.dims.coils;
end
if in.dims.averages>in.dims.subSpecs
    dims.averages=in.dims.averages-1;
else
    dims.averages=in.dims.averages;
end
dims.subSpecs=0;
if in.dims.extras>in.dims.subSpecs
    dims.extras=in.dims.extras-1;
else
    dims.extras=in.dims.extras;
end

%re-calculate the sz variable
sz=size(fids);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;
out.subspecs=1;
out.averages=in.averages/2;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.subtracted=1;
