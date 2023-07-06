function out=FIDAfft(in,dim,in_domain)
%Helper function: performs fft & fftshift on your data that's in the time 
%("t") or frequency domain ("f")
%Peter Truong, 2023
%   INPUTS:
%           in       :   data, either FID or SPEC
%           dim      :   dimension on matrix where you will transform
%           in_domain:   "t" or "f"; domain of "in" - time or frequency,
%                     respectively
%
%   OUTPUT:
%           out   :   the transformed data
    switch in_domain
        case 't'
            out = fftshift(fft(in,[],dim),dim);
        case 'f'
            out = ifft(fftshift(in, dim),[],dim);
        otherwise
            disp('Invalid input for "domain". Please either input "t" or "f" for domain of your input data, respectively');
            return
    end
end