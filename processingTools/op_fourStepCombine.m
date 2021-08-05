% op_fourStepCombine.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_fourStepCombine(in,addInd);
% 
% DESCRIPTION:
% Combine dimensions in datasets with 4 subspectra to produce
% a new dataset with only 2 subspecs.  Common applications of this function 
% are to combine the ISIS sub-spectra of the MEGA-SPECIAL sequence (but not 
% (the Edit-on/edit-off subspectra), or to combine the A-B and C-D pairs of 
% a hadamard encoded MRS acqusition (HERMES or HD-SPECIAL).   
% 
% INPUTS:
% in         = input data in matlab structure format.
% mode       = (optional) 
%                   0 (default) - row indices [1 2] and [3 4] will be added.  
%                   1 - row indices [1 2] and [3 4] will be subtracted.
%                   2 - row indices [1 3] and [2 4] will be added.
%                   3 - row indices [1 3] and [2 4] will be subtracted. 
%
% OUTPUTS:
% out        = Output following combination of subspectra.  

function out=op_fourStepCombine(in,mode);



if ~in.flags.isFourSteps
    error('ERROR:  requires a dataset with 4 subspecs as input!  Aborting!');
end
if in.sz(end)~=4
    error('ERROR: final matrix dim must have length 4!!  Aborting!');
end


if nargin<2
    mode=0;
end

%now make subspecs and subfids (This doesn't do anything to MEGA-PRESS
%data, but it combines the SPECIAL iterations in MEGA-SPECIAL).
sz=in.sz;
fids=in.fids;
reshapedFids=reshape(fids,prod(sz(1:end-1)),sz(end));
ind1=[1 2];
ind2=[3 4];
ind3=[1 3];
ind4=[2 4];
sz(end)=sz(end)-2;
if mode==0
    reshapedFids(:,1)=sum(reshapedFids(:,ind1),2);
    reshapedFids(:,2)=sum(reshapedFids(:,ind2),2);
elseif mode==1
    reshapedFids(:,1)=diff(reshapedFids(:,ind1),1,2);
    reshapedFids(:,2)=diff(reshapedFids(:,ind2),1,2);
elseif mode==2
    reshapedFids(:,1)=sum(reshapedFids(:,ind3),2);
    reshapedFids(:,2)=sum(reshapedFids(:,ind4),2);
elseif mode==3
    reshapedFids(:,1)=diff(reshapedFids(:,ind3),1,2);
    reshapedFids(:,2)=diff(reshapedFids(:,ind4),1,2);
else
    error('ERROR: mode not recognized. Value must be 0, 1, 2 or 3');
end
fids=reshape(reshapedFids(:,ind1),sz);
fids=fids/2;  %Divide by 2 so that this is an averaging operation;


%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.subspecs=out.sz(out.dims.subSpecs);

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.isFourSteps=0;

