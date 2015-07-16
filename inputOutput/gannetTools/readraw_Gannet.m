function [d] = readraw_Gannet(filename, type, ds, endian,offset)

if( nargin == 3 )
	endian = 'b';
end

fp=fopen(filename, 'rb', endian);
%fseek(fp, 0, -1);
%fseek(fp, 262144, -1);
fseek(fp, offset, -1);
d = fread(fp, prod(ds), type);
fclose(fp);
d = reshape(d, ds);
