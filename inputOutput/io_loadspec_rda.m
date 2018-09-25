%io_loadspec_rda.m
%Jay Hennessy, McGill University 2016.
%
% USAGE:
% [out]=io_loadspec_rda(rda_filename);
% 
% DESCRIPTION:
% Reads in siemens rda data (.rda file).
% 
% op_loadspec_rda outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Siemens rda data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function [out] = io_loadspec_rda(rda_filename)

fid = fopen(rda_filename);

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';

tline = fgets(fid);

while (isempty(strfind(tline , head_end_text)))
    
    tline = fgets(fid);
    
    if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)
        
        
        % Store this data in the appropriate format
        
        occurence_of_colon = findstr(':',tline);
        variable = tline(1:occurence_of_colon-1) ;
        value    = tline(occurence_of_colon+1 : length(tline)) ;
        
        switch variable
        case { 'PatientID' , 'PatientName' , 'StudyDescription' , 'PatientBirthDate' , 'StudyDate' , 'StudyTime' , 'PatientAge' , 'SeriesDate' , ...
                    'SeriesTime' , 'SeriesDescription' , 'ProtocolName' , 'PatientPosition' , 'ModelName' , 'StationName' , 'InstitutionName' , ...
                    'DeviceSerialNumber', 'InstanceDate' , 'InstanceTime' , 'InstanceComments' , 'SequenceName' , 'SequenceDescription' , 'Nucleus' ,...
                    'TransmitCoil' }
            eval(['rda.' , variable , ' = value; ']);
        case { 'PatientSex' }
            % Sex converter! (int to M,F,U)
            switch value
            case 0
                rda.sex = 'Unknown';
            case 1
                rda.sex = 'Male';
            case 2
                
                rda.sex = 'Female';
            end
            
        case {  'SeriesNumber' , 'InstanceNumber' , 'AcquisitionNumber' , 'NumOfPhaseEncodingSteps' , 'NumberOfRows' , 'NumberOfColumns' , 'VectorSize' }
            %Integers
            eval(['rda.' , variable , ' = str2num(value); ']);
        case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
                     'SliceThickness' ,  'FoVHeight' , 'FoVWidth' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol'}
            %Floats 
            eval(['rda.' , variable , ' = str2num(value); ']);
        case {'SoftwareVersion[0]' }
            rda.software_version = value;
        case {'CSIMatrixSize[0]' }
            rda.CSIMatrix_Size(1) = str2num(value);    
        case {'CSIMatrixSize[1]' }
            rda.CSIMatrix_Size(2) = str2num(value);    
        case {'CSIMatrixSize[2]' }
            rda.CSIMatrix_Size(3) = str2num(value);    
        case {'PositionVector[0]' }
            rda.PositionVector(1) = str2num(value);    
        case {'PositionVector[1]' }
            rda.PositionVector(2) = str2num(value);     
        case {'PositionVector[2]' }
            rda.PositionVector(3) = str2num(value);    
        case {'RowVector[0]' }
            rda.RowVector(1) = str2num(value);    
        case {'RowVector[1]' }
            rda.RowVector(2) = str2num(value);       
        case {'RowVector[2]' }
            rda.RowVector(3) = str2num(value);    
        case {'ColumnVector[0]' }
            rda.ColumnVector(1) = str2num(value);     
        case {'ColumnVector[1]' }
            rda.ColumnVector(2) = str2num(value);       
        case {'ColumnVector[2]' }
            rda.ColumnVector(3) = str2num(value);    
            
        otherwise
            % We don't know what this variable is.  Report this just to keep things clear
            %disp(['Unrecognised variable ' , variable ]);
        end
        
    else
        % Don't bother storing this bit of the output
    end
    
end

%
% So now we should have got to the point after the header text
% 
% Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
%

bytes_per_point = 16;
complex_data = fread(fid , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) *rda.CSIMatrix_Size(1) *rda.VectorSize * 2 , 'double');  
%fread(fid , 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
fclose(fid);

% Now convert this data into something meaningful

 %Reshape so that we can get the real and imaginary separated
 hmm = reshape(complex_data,  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );
 
 %Combine the real and imaginary into the complex matrix
 fids = complex(hmm(1,:,:,:,:),hmm(2,:,:,:,:));
 fids = fids';
 %Remove the redundant first element in the array
% Time_domain_data = reshape(hmm_complex, rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3));
% Time_domain_data = fft(fft(fftshift(fftshift(padarray(fftshift(fftshift(ifft(ifft(Time_domain_data,[],2),[],3),2),3),[0 16 16],'both'),2),3),[],2),[],3);

% get the spectrum from the fid
specs = fftshift(ifft(fids,[],1),1);

% make calculations for the output mrs structure
sz = rda.VectorSize;
dwelltime = rda.DwellTime/1000000;
spectralwidth=1/dwelltime;
txfrq = rda.MRFrequency*1000000;
dims.t = 1;
dims.subSpecs = 0;
Bo = rda.MagneticFieldStrength;
rawAverages = rda.NumberOfAverages;
averages = 1;
subspecs =1;
rawSubspecs = 'na';
date = rda.StudyDate;
seq = rda.SequenceDescription;
TE = rda.TE;
TR = rda.TR;
pointsToLeftShift = 'N/A';

%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.6082;

t=[0:dwelltime:(sz(1)-1)*dwelltime];

%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=TE;
out.tr=TR;
out.pointsToLeftshift=pointsToLeftShift;


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end
      
end
