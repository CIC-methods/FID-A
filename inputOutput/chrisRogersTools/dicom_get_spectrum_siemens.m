% A special function to read Siemens' spectroscopy FIDs.
% Copyright Chris Rodgers, University of Oxford.
function complex_fid = dicom_get_spectrum_siemens(fd) 

% advance the field to the appropriate place
% in the file
field_length = dicom_move(fd, '7FE1', '1010');

field_size = field_length / 4;

% we can use fread to read data in as floats
[fid, fid_size] = fread(fd, field_size, 'float32','ieee-le');
 
if( fid_size ~= field_size )
     fprintf('\nWarning: field size was %d and %d elements were read.', field_size, fid_size);
end

real_part = zeros(length(fid)/2, 1);
imag_part = real_part;

% sort into two columns or make complex
k = 1;
for n = 1:1:length(fid)
    if mod(n,2)
        real_part(k) = fid(n);
    else
        imag_part(k) = fid(n);
        k = k + 1;
    end
end

complex_fid = real_part + j*imag_part;
