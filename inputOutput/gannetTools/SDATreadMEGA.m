function A = SDATreadMEGA(filename,da_xres, da_yres)
%function A = SDATread(filename,da_xres)
%i = sqrt(-1);
% Open file to read reference scan data.
fid = fopen(filename,'rb', 'ieee-le');
if fid == -1
   err_msg = sprintf('Unable to locate File %s', filename)
   return;
end
%%Set up a structure to take the data:
A=zeros(da_xres,1);
 totalframes=1;
 totalpoints = totalframes*da_xres*da_yres*2;

    raw_data = freadVAXG(fid, totalpoints, 'float'); %NP26-11-13 freadVAXG not found in Gannet2.0 folder. Still necessary to have original Gannet added to path
 
    TempData=reshape(raw_data,[2 da_xres da_yres]);
    A=squeeze(TempData(1,:,:)+1i*TempData(2,:,:));
 fclose(fid);
end