% op_ISIScombine.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_ISIScombine(in,addInd);
% 
% DESCRIPTION:
% Combine dimensions corresponding to ISIS on/off acquisitions to produce
% fully localized MRS volumes.  Mostly used for MEGA-SPECIAL data to
% combine the ISIS subspectra but not the Edit-on/edit-off subspectra.
% 
% INPUTS:
% in         = input data in matlab structure format.
% addInd     = (optional) If add==1, then row indices [1 2] and [3 4]
%             will be added.  Otherwise, row indices [1 2] and [3 4] will 
%             be subtracted.
%
% OUTPUTS:
% out        = Output following combination of ISIS subspectra.  

function out=op_ISIScombine(in,addInd);



if ~in.flags.isISIS
    error('ERROR:  requires a MEGA-SPECIAL dataset as input!  Aborting!');
end
if in.sz(end)~=4
    error('ERROR: final matrix dim must have length 4!!  Aborting!');
end


if nargin<2
    addInd=0;
end

%now make subspecs and subfids (This doesn't do anything to MEGA-PRESS
%data, but it combines the SPECIAL iterations in MEGA-SPECIAL).
sz=in.sz;
fids=in.fids;
reshapedFids=reshape(fids,prod(sz(1:end-1)),sz(end));
ind1=[1 2];
ind2=[3 4];
sz(end)=sz(end)-2;
if addInd==0
    reshapedFids(:,1)=sum(reshapedFids(:,ind1),2);
    reshapedFids(:,2)=sum(reshapedFids(:,ind2),2);
else
    reshapedFids(:,1)=diff(reshapedFids(:,ind1),1,2);
    reshapedFids(:,2)=diff(reshapedFids(:,ind2),1,2);
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
out.flags.isISIS=0;

