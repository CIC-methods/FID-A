% op_dsv2timecourse.m
% Matthew David Robson, Oxford University, 2008
% Jamie Near, Sunnybrook Research Institute, 2026
%
% USAGE:
% output_data = op_dsv2timecourse(dsv);
%
% DESCRIPTION:
% Uncompress a siemens dsv waveform (read using io_readdsv.m) into a
% timecourse (basically unpack the dsv data).  
%
% INPUTS: 
% dsv           = dsv structure, obtained by running io_readdsv.m on the raw
%                 Siemens .dsv file. 
% 
% OUTPUTS: 
% output_data   = Uncompressed, digital waveform describing the pulse
%                 seuqence channel of interest.  

function output_data=op_dsv2timecourse(dsv);

in=dsv.output_data;
N=length(in);


n=1;
out(1)=in(n);
n=n+1;


while n < length(dsv.output_data)-1
    if in(n)==in(n-1)
        out(end+1:end+1+in(n+1))=ones(1,in(n+1)+1)*in(n);
        n=n+2;
        out(end+1)=in(n);
        n=n+1;
        %disp(num2str(n));
    else
        out(end+1)=in(n);
        n=n+1;
        %disp(num2str(n));
    end
end

output_data=cumsum(out);
        


