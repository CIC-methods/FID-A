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

function [metabs,corrMatrix]=io_loadlcmdetail(filename);


fid=fopen(filename);
line=fgets(fid);
search_index=findstr(line,'Correlation coefficients');



while isempty(search_index)
    line=fgets(fid);
    search_index=findstr(line,'Correlation coefficients');
end


line=fgets(fid);
line=fgets(fid);
metabstr=line;
line=fgets(fid);
metabstr=[metabstr line];
line=fgets(fid);
metabstr=[metabstr line];

remain=metabstr;
n=1;
while ~isempty(remain)
    [metabs{n},remain]=strtok(remain);
    n=n+1;
end

corrMatrix=zeros(numel(metabs)-1);
line=fgets(fid);
fillLine=1;

corrMatrix(1,1)=0.5;
while fillLine<18
    [name,corrs]=strtok(line);
    corrs=[corrs(1:end-1) ' 0.5'];
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end

line=fgets(fid);

while fillLine<35
    [name,corrs]=strtok(line);
    line=fgets(fid);
    corrs=[corrs(1:end-1) ' ' line(1:end-1) ' 0.5'];
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end

while fillLine<numel(metabs)
    [name,corrs]=strtok(line);
    line1=fgets(fid);
    line2=fgets(fid);
    corrs=[corrs(1:end-1) ' ' line1(1:end-1) ' ' line2(1:end-1) ' 0.5'];
    size(str2num(corrs));
    corrMatrix(fillLine+1,1:length(str2num(corrs)))=str2num(corrs);
    line=fgets(fid);
    fillLine=fillLine+1;
end
metabs{end}=name;

corrMatrix=corrMatrix+corrMatrix';

fclose(fid);

% imagesc(corrMatrix);
% impixelinfo;


    







