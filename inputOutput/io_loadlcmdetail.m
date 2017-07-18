%io_loadlcmdetail.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [metabs,corrMatrix]=io_loadlcmdetail(filename);
% 
% DESCRIPTION:
% This function loads in the "detailed output" of LCModel and returns the
% matrix of metabolite correlation coefficients.
% 
% INPUTS:
% filename   = Filename of the lcmodel detailed output file.  
%
% OUTPUTS:
% metabs     = A listing of the metabolites included in the correlation
%               coefficients table.
% corrMatrix = A matrix of correlation coefficients between metabolites,
%               with indices specified by the 'metabs' variable.  

function [metabs,corrMatrix]=io_loadlcmdetail(filename);

%First open the detailed output file and search for the beginning of the
%Correlation Coefficients table
fid=fopen(filename);
line=fgets(fid);
search_index=findstr(line,'Correlation coefficients');

while isempty(search_index)
    line=fgets(fid);
    search_index=findstr(line,'Correlation coefficients');
end

%Now that we've found the top of the correlation coefficients table, we
%need to record the names of all of the metabolites that are listed in the
%table.  Start by going down two lines to where the metabolite names are
%listed:
line=fgets(fid);
line=fgets(fid);

%Read the first line of metabolite names and then go down one line.  The 
%table of correlation coefficients can only fit 17 metabolites per row.
%Therefore if there are more than 17 metabolites to read, then we will have
%to go down another row and read again.  A simple trick to know whether to
%go down another row is to check and see if the first token in the next row
%is equal to the second token from the first row.  If so, then we're done
%reading metabolite names.  If not, then we need to keep reading.  
metabstr=line;
[temp,remain]=strtok(line);
secondMetabName=strtok(remain);
line=fgets(fid);

while ~strcmp(strtok(line),secondMetabName)
    metabstr=[metabstr line];
    line=fgets(fid);
end

%We are now finished reading all of the metabolite names.  Now make a cell 
%array containing the names of all of the meatbolites.  Note that at this
%point in the code, the final element in metabs will be a blank space,
%since the last metabolite name is missing from the top row of the correlation
%coefficient table in the LCModel detailed output.
remain=metabstr;
n=1;
while ~isempty(remain)
    [metabs{n},remain]=strtok(remain);
    n=n+1;
end

%Now calculate the total number of metabolites, and initialize the 
%correlation coefficient matrix:
numMetabs=numel(metabs);
corrMatrix=zeros(numMetabs);
fillLine=1;

%Now start populating the correlation coefficient matrix.  The first token
%of each line is a metabolite name, so we have to remove it.  The rest of
%the line is a row in the correlation matrix.  The first 17 lines are
%simple to read becuase they all have 17 elements or less, so they don't
%wrap around the page (remember there are only 17 columns in the page).  NOTE:
%The correlation matrix is symmetric about the diagonal, but we start by only
%constructing the lower half of the diagonal (The upper half will be all 
%zeros).  The diagonal elements themselves will be set to 0.5.  Then later, 
%the full matrix, will be completed by adding the matrix to it's transpose,
%thus completing the symmetry, and making all of the diagonals ones.  We
%will also rename "metabs" as we go, becuase the full names of each
%metabolite are given with greater precision on the left column of the
%table in the lcmodel detailed output (compared to the top row of the table).
corrMatrix(1,1)=0.5;
if numMetabs>17
    upto=18;
else
    upto=numMetabs;
end
while fillLine<upto
    [name,corrs]=strtok(line);
    metabs{fillLine+1}=name;
    corrs=[corrs(1:end-1) ' 0.5'];
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end

    

line=fgets(fid);

%For the 18th to 34th lines, there is one line of 'wrap around'.
if numMetabs>34
    upto=35;
else
    upto=numMetabs;
end
while fillLine<upto
    [name,corrs]=strtok(line);
    metabs{fillLine+1}=name;
    line=fgets(fid);
    corrs=[corrs(1:end-1) ' ' line(1:end-1) ' 0.5'];
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end


%For the 35th to 52nd lines, there are two lines of 'wrap around'.  Assume
%that there will be no more than 52 elemens in the correlation matrix.
while fillLine<numMetabs
    [name,corrs]=strtok(line);
    metabs{fillLine+1}=name;
    line1=fgets(fid);
    line2=fgets(fid);
    corrs=[corrs(1:end-1) ' ' line1(1:end-1) ' ' line2(1:end-1) ' 0.5'];
    size(str2num(corrs));
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end
metabs{end}=name;

%Now add the correlation matrix (the lower have of which is currently only
%populated) to its transpose to complete the symmetry.
corrMatrix=corrMatrix+corrMatrix';

fclose(fid);

% imagesc(corrMatrix);
% impixelinfo;


    







