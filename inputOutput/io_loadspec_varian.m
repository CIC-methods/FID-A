%io_loadspec_varian.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_loadspec_varian(filename);
% 
% DESCRIPTION:
% Reads in varian .fid data using the readfid.m and readprocpar.m functions 
% from Martyn Klassen (mklassen@robarts.ca).
% 
% io_loadspec_varian outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Varian .fid data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_varian(filename);

%read in the data using read_meas_dat
%dOut=op_read_meas_dat2(filename);
%[par,img,k]=fidread(filename);
[fids,hdr,block_hdr]=readfid(filename);
par=readprocpar(filename);

%get the echo time
te=par.te*1000;

%get the repetition time
tr=par.tr*1000;

%get the sequence id;
sequence=par.pslabel{1};

fids=squeeze(fids);
sz=size(fids);

%NOW WE HAVE TO SORT OUT WHAT DIMENSIONS ARE WHICH.  FOR NOW, THIS WILL BE
%DONE USING SOME USER INTERACTION:

fprintf(['\nTHE SIZE OF THE DATA ARRAY IS:  ']);
k=size(fids);
while length(k)>0
    if length(k)>1
        fprintf([num2str(k(1)) ' x ' ]);
        k=k(2:end);
    else
        fprintf([num2str(k(1)) ' \n\n' ]);
        k=[];
    end
end

fprintf('Now please identify each of the data dimensions.  Note, if any two \n');
fprintf('data dimensions have the same array index, you will also need to \n');
fprintf('specify whether they are interleaved or not.  For example, if averages \n');
fprintf('and coils are both indexed in the 2nd diemension of the array, please \n');
fprintf('answer ''2'' below for both averages and coils, and then answer ''y'' or ''n'' to \n');
fprintf('the question about whether they are interleaved or not. \n\n');

dims.t=input('Which is the time dimension? (Usually it is ''1''):  ');
dims.coils=input('Which is the coils dimension? (''0'' for none):  ');
dims.averages=input('Which is the averages dimension? (''0'' for none):  ');
dims.subSpecs=input('Which is the subspecs dimension? (''0'' for none):  ');
dims.extras=input('Any extra dimensions not listed above?  (''0'' for none):  ');

%Now check if any dimensions are stored along the same array index.  If so,
%then separate them.
dimensions=[dims.t dims.coils dims.averages dims.subSpecs dims.extras];
dimensionNames={'time','coils','averages','subSpecs','extras'};
nonZeroDimIndices=(dimensions>0);
nonZeroDimensionNames=dimensionNames(nonZeroDimIndices);
nonZeroDimensions=dimensions(nonZeroDimIndices);
[nonZeroDimensionsSorted,I]=sort(nonZeroDimensions);
diffs=diff(nonZeroDimensionsSorted);
if ~all(diffs)
    fprintf('One data array index contains more than one dimension of data!!!\n\n');
    zeroIndices=find(diffs==0);
    first=I(zeroIndices(1));
    second=I(zeroIndices(1)+1);
    arrayIndexToBeSeparated=nonZeroDimensionsSorted(zeroIndices(1));
    firstDoubledDim=nonZeroDimensionNames{first};
    secondDoubledDim=nonZeroDimensionNames{second};
    disp(['Specifically, array dimension number ' num2str(arrayIndexToBeSeparated) ' contains both the ''' firstDoubledDim ''' and ''' secondDoubledDim ''' dimensions.']);
    interleaved=input(['Please indicate whether ' firstDoubledDim ' and ' secondDoubledDim ' are interleaved (''y'') or not (''n''):  '],'s');

    
    %NOW FIGURE OUT THE SIZE OF EACH OF THESE TWO DIMENSIONS
    A=[];
    B=[];
    if strcmp(firstDoubledDim,'coils') || strcmp(secondDoubledDim,'coils')
        A=par.nrcvrs;
        dA='coils';
    end
    if strcmp(firstDoubledDim,'averages') || strcmp(secondDoubledDim,'averages')
        if isempty(A)
            A=max(par.nt,length(par.nt));
            dA='averages';
        else
            if length(par.nt)>1
                B=length(par.nt);
            elseif length(par.nt)==1
                B=par.nt;
            end
            dB='averages';
        end
    end
    if isempty(B)
        if strcmp(firstDoubledDim,'subSpecs') || strcmp(secondDoubledDim,'subSpecs')
            if isempty(A)
                A=input('Please input the number of subSpectra that are stored in this dataset:  ');
                dA='subSpecs';
            else
                B=input('Please input the number of subSpectra that are stored in this dataset:  ');
                dB='subSpecs';
            end
        end
        if isempty(B)
            if strcmp(firstDoubledDim,'extras') || strcmp(secondDoubledDim,'extras')
                B=input('Please input the number of extra dimensions that are stored in this dataset:  ');
                dB='extras';
            end
        end
    end
    
    %NOW CHECK IF A x B is equal to the size of the double dimension:
    if ~(A*B==size(fids,arrayIndexToBeSeparated));
        error('ERROR:  Array dimension size does not match header, or information given by user.');
    end

    %NOW SEPARATE THE DIMENSIONS THAT ARE STORED TOGETHER:
    if arrayIndexToBeSeparated==2
        if max(dimensions)==2
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],n)=fids(:,[n:min(A,B):sz(2)]);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],n)=fids(:,[(n-1)*max(A,B)+1:n*max(A,B)]);
                end
            end
            if A<B
                eval(['dims.' dA '=3;']);
            else
                eval(['dims.' dB '=3;']);
            end
        elseif max(dimension)==3
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],:,n)=fids(:,[n:min(A,B):sz(2)],:);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],:,n)=fids(:,[(n-1)*max(A,B)+1:n*max(A,B)],:);
                end
            end
            if A<B
                eval(['dims.' dA '=4;']);
            else
                eval(['dims.' dB '=4;']);
            end
        elseif max(dimensions)==4
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],:,:,n)=fids(:,[n:min(A,B):sz(2)],:,:);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,[1:1:max(A,B)],:,:,n)=fids(:,[(n-1)*max(A,B)+1:n*max(A,B)],:,:);
                end
            end
            if A<B
                eval(['dims.' dA '=5;']);
            else
                eval(['dims.' dB '=5;']);
            end
        end
    elseif arrayIndexToBeSeparated==3
        if max(dimension)==3
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,:,[1:1:max(A,B)],n)=fids(:,:,[n:min(A,B):sz(3)]);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,:,[1:1:max(A,B)],n)=fids(:,:,[(n-1)*max(A,B)+1:n*max(A,B)]);
                end
            end
            if A<B
                eval(['dims.' dA '=4;']);
            else
                eval(['dims.' dB '=4;']);
            end
        elseif max(dimensions)==4
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,:,[1:1:max(A,B)],:,n)=fids(:,:,[n:min(A,B):sz(3)],:);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,:,[1:1:max(A,B)],:,n)=fids(:,:,[(n-1)*max(A,B)+1:n*max(A,B)],:);
                end
            end
            if A<B
                eval(['dims.' dA '=5;']);
            else
                eval(['dims.' dB '=5;']);
            end
        end
    elseif arrayIndexToBeSeparated==4
        if max(dimensions)==4
            if strcmp(interleaved,'y');
                for n=1:min(A,B)
                    tempfids(:,:,:,[1:1:max(A,B)],n)=fids(:,:,:,[n:min(A,B):sz(4)]);
                end
            else
                for n=1:min(A,B)
                    tempfids(:,:,:,[1:1:max(A,B)],n)=fids(:,:,:,[(n-1)*max(A,B)+1:n*max(A,B)]);
                end
            end
            if A<B
                eval(['dims.' dA '=5;']);
            else
                eval(['dims.' dB '=5;']);
            end
        end
    end
    fids=squeeze(tempfids);
    %fids=fids([end:-1:1],:,:);
end



sz=size(fids);

if isa(fids,'int32')
    fids=double(fids);
end
specs=fftshift(ifft(fids,[],dims.t),dims.t);

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
spectralwidth=par.sw;
dwelltime=1/spectralwidth;

    
%Get TxFrq
txfrq=par.sfrq*1e6;


%Get Date
date=par.date;



%GET THE NUMBER OF RAW AVERAGES AND SUBSPECTRA:

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if length(par.nt)>1
    averages=length(par.nt);
    rawAverages=length(par.nt);
else
    averages=1;
    rawAverages=par.nt;
end


%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end


%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(txfrq/1e6);
ppm=ppm+4.65;

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
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.date=date;
out.dims=dims;
out.Bo=out.txfrq/42.577/1e6;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
if dims.averages>0
    out.flags.averaged=0;
else
    out.flags.averaged=1;
end
if dims.coils>0
    out.flags.addedrcvrs=0;
else
    out.flags.addedrcvrs=1;
end
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.isISIS=0;



%DONE
