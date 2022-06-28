%io_writeRFbruk.m
%Jamie Near, Ved Hatolkar, Sunnybrook Research Institute 2022.
%
% USAGE:
% RF=io_writeRFbruk(rf,outfile);
% 
% DESCRIPTION:
% Write a matlab RF pulse structure (containing N x 3 waveform array field 
% with rf.waveform(:,1)= phase, rf.waveform(:,2)=amplitude, and rf.waveform(:,3)=timestep), 
% to a Bruker format rf pulse file.
% 
% INPUTS:
% rf         = matlab RF pulse.
% outfile    = name of the output .pta file to be written.
%
% OUTPUTS:
% RF            = Same as input.  Not used.  The primary output of this
%                   function is a text file in Siemens .pta format. 

function RF=io_writeRFbruk(rf,outfile);

B1INT=sum(rf.waveform(:,2));

l=length(rf.waveform(:,2));
RF.waveform(:,1)=100*rf.waveform(:,2)/max(rf.waveform(:,2));
RF.waveform(:,2)=mod(rf.waveform(:,1),360);

size(RF.waveform)


minx=min(RF.waveform(1,:));
maxx=max(RF.waveform(1,:));
miny=min(RF.waveform(2,:));
maxy=max(RF.waveform(2,:));


%We need to fill in the INTEGFAC value, but we are not sure how to
%calculate this yet, so we will enter it manually for now.  

IFcalc=input('Would you like to perform an automatic INTEGFAC Calculation?(y or n)  ', 's');

if IFcalc=='y' || IFcalc=='Y'
    %FIX THIS LATER
    INTEGFAC=sqrt((sum(real(RF_complex)).^2)+(sum(imag(RF_complex)).^2));
elseif IFcalc=='n' || IFcalc=='N'
    % If the pulse is not a plane wave (Adiabatic pulses, off resonance pulses,
    % or dual Band Pulses), you should calculate the INTEGFAC
    % yourself and input it manually.  For methods of calculating INTEGFAC,
    % see bruker manual??
    INTEGFAC=input('Now you must input the INTEGFAC Manually.  Lets have it?  ');
else
    error('ERROR!  You must input a value for INTEGFAC!  ABORTING!');
end


%Reference gradient is only used in slice selective pulses.  This refers to
%the gradient required to excite a 10mm slice if the pulse duration is
%5.12ms  Given that the Time-bandwidth product is included in the rf
%structure, the refgrad value can be calculated here, but the option is
%provided to enter it manually, in case a different refgrad value is
%desired.  

refgrad=1;
refg=input('Would you like to input a reference gradient (REFGRAD)? ','s');
if refg=='y' || refg=='Y'
    disp('INFORMATION:  REFGRAD is the gradient amplitude in [mT/m] required ');
    disp('to excite a 10mm slice with a pulse duration of 5.12 ms. ');
    man=input('Calculate refgrad automatically (y or n) ? ','s');
    if man=='y' || man=='Y'
        %refgrad [mT/m] = (1000 [mT]/[1T]) * tbw [unitless] ) / ( (0.00512 [s]) * (42577000 [Hz]/[T]) * (0.01 [m]) );
        refgrad = 1000 * rf.tbw/(0.00512 * 42577000 * 0.01);
    elseif man=='n' || man=='N'
        refgrad=input('Manually enter the reference gradient (mT/m):  ');
    else
        error('ERROR:  Response not recognized!!');
    end
end

todaysdate=date;
currenttime=datestr(now,'HH:MM:SS.fff');

%write to bruker format rf pulse file
fid=fopen(outfile,'w+');
fprintf(fid,['##TITLE= ' outfile]);
fprintf(fid,'\n##JCAMP-DX= 5.00 Bruker JCAMP library');
fprintf(fid,'\n##DATA TYPE= Shape Data');
fprintf(fid,'\n##ORIGIN= Bruker 7T');
fprintf(fid,'\n##OWNER= jn_vh');
% fprintf(fid,'\n');
fprintf(fid,['\n##DATE= ' todaysdate]);
fprintf(fid,['\n##TIME= ' currenttime]);
fprintf(fid,'\n##MINX= %2.6e',minx);
fprintf(fid,'\n##MAXX= %2.6e',maxx);
fprintf(fid,'\n##MINY= %2.6e',miny);
fprintf(fid,'\n##MAXY= %2.6e',maxy);
if (rf.type=='inv')
    fprintf(fid,'\n##$SHAPE_EXMODE= Inversion');
    fprintf(fid,'\n##$SHAPE_TOTROT= 18.00000e+1');
elseif (rf.type=='exc')
    fprintf(fid,'\n##$SHAPE_EXMODE= Excitation');
    fprintf(fid,'\n##$SHAPE_TOTROT= 90.00000e0');
elseif (rf.type=='ref')
    fprintf(fid,'\n##$SHAPE_EXMODE= Refocusing');
    fprintf(fid,'\n##$SHAPE_TOTROT= 18.00000e+1');
end 
fprintf(fid,'\n##$SHAPE_BWFAC= %2.6e', rf.tbw);
fprintf(fid,'\n##$SHAPE_INTEGFAC= %2.6e', INTEGFAC);
fprintf(fid,'\n##$SHAPE_REPHFAC=0');
fprintf(fid,'\n##$SHAPE_TYPE=adiabatic');
fprintf(fid,'\n##$SHAPE_MODE=0');
fprintf(fid,'\n##NPOINTS= %d ', size(rf.waveform,1));
fprintf(fid,'\n##XYPOINTS= (XY..XY)');
n=0;
fprintf(fid,'\n');
fprintf(fid,'%1.6e, %1.6e\n',RF.waveform(:,1:2)');
fclose(fid);