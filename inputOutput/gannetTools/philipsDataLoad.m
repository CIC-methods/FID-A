function [FullData,WaterData] = philipsDataLoad( fname, fname_water )
% RE/CJE Parse SPAR file for header info
% 110825

%ADH 20130701 - modify such that if the water reference scan is loaded in
%the same file, it works. for this case, will only have one input fname
%because all the raw data is in the same place. Next step will be to put
%a flag to see if there is water data in a separate scan, in the same scan
%as a water ref (as called in SDAT) or if its just not been collected at all

% will need to get all the indicies first because will need to know how big
% the waterdata is before loading in the GABAdata.



% work out data header name
sparname = [fname(1:(end-4)) 'list'];
sparheader = textread(sparname, '%s');

%ADH - decide if there is water data as ref data included in the data and
%if so, set a flag to pull it out properly...
sparidx=find(ismember(sparheader, 'number_of_mixes')==1);
nu_mixes = str2num(sparheader{sparidx+2});
%   MRS_struct.nu_mixes = nu_mixes;

if nu_mixes == 2
    % ADH 20130708 - if have 2 nu_mixes - that mean there is water data
    % within the data and need to separate it our properly
    Reference_compound = 'H2O'; % ADH because actually is water, its jsut not in a separate file
    
    sparidx=find(ismember(sparheader, 'F-resolution')==1, 1, 'first'); % gaba_data - ADH 20130701
    npoints = str2num(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'F-resolution')==1, 1, 'last'); %water_data- ADH 20130701
    npoints_water = str2num(sparheader{sparidx+2});
    
    sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1, 1, 'first');
    nrows = str2num(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1, 1, 'last');
    nrows_water = str2num(sparheader{sparidx+2}); % this might not be correct because we dont collect all the water data - ADH 20130701
    nrows_water = 1; % ADH 20130701
    % set = 1 because only collect 1 water scan but the
    % header thinks there will be as many rows as in the
    % gaba data
    
    
    %ADH not completely sure why Navg changes across the structure of data
    %reads and others don't but for now keep this characteristic... may
    %screw up later though - something to keep an eye out for. - ADH 20130701
    sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1, 1, 'first');
    Navg = str2num(sparheader{sparidx+2})*nrows;
    sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1, 1, 'last');
    Navg_water = str2num(sparheader{sparidx+2}); %Trial SDAT might be average not sum.
    %Need to determine the offset of the data - i.e. how many channels are
    %there...
    sparidx=find(ismember(sparheader, 'NOI')==1); %ADH skip past noise
    coil_channels=size(sparidx,1)-2;
    sparidx=find(ismember(sparheader, 'STD')==1);
    
    %MRS_struct.p.ptr_offset_water=str2num(sparheader{sparidx(3)+20});
    %MRS_struct.p.ptr_offset_data=str2num(sparheader{sparidx(3+MRS_struct.Navg_water(MRS_struct.ii))+20});
    ptr_offset_water=str2num(sparheader{sparidx(3)+20});
    ptr_offset_data=str2num(sparheader{sparidx(3+Navg_water)+20});
    %ADH the actual reading of the data happens with readraw_Gannet - just
    % takes in the real data to analyse so need to decide on the offsets
    % proior to this. will need to decide it for both the water data and the
    % GABA data.
    % Use 2 readraw_Gannet - one for water and one for data - need to sort
    % out the skips of data...- ADH 20130701
    
    %Need to skip rows associated with the '
    data = readraw_Gannet(fname, 'float', ...
        [2 npoints 1 Navg/nrows nrows], 'l',ptr_offset_data);
    %  Make data complex.
    data = squeeze(data(1,:,:,:,:)+ 1i*data(2,:,:,:,:));
    
    
    % then the same for the water data- ADH 20130701
    data_water = readraw_Gannet(fname, 'float', ...
        [2 npoints_water 1 Navg_water nrows_water], 'l',ptr_offset_water);
    %  Make data complex.- ADH 20130701
    data_water = squeeze(data_water(1,:,:,:)+ 1i*data_water(2,:,:,:));
    
    % % this part is from the original PhilipsRead_data needs to be
    % % removed because not a 2nd waterfile- ADH 20130701
    
    %%'Duplicated' code added in to handle .data water files
    % work out water data header name
    %
    %    if nargin >2
    %        sparnameW = [fname_water(1:(end-4)) 'list'];
    %        sparheaderW = textread(sparnameW, '%s');
    %        sparidx=find(ismember(sparheaderW, 'number_of_signal_averages')==1);
    %        MRS_struct.p.Nwateravg(MRS_struct.ii) = str2num(sparheaderW{sparidx+2}); %Trial SDAT might be average not sum.
    %
    %
    %        %Need to skip rows associated with the '
    %
    %        MRS_struct.fids.data_water = readraw_Gannet(fname_water, 'float', [2 MRS_struct.npoints MRS_struct.p.Nwateravg(MRS_struct.ii) 1], 'l', MRS_struct.p.ptr_offset);
    %
    %        %  Make data complex.
    %        MRS_struct.fids.data_water = squeeze(MRS_struct.fids.data_water(1,:,:,:,:)+ 1i*MRS_struct.fids.data_water(2,:,:,:,:));
    %    end
    %
    
    
    
    %undo time series phase cycling to match GE
    corrph = ones(size(data));
    for jj=1:size(data,3)
        corrph(:,:,jj) = corrph(:,:,jj) * (-1).^(jj+1);
    end
    
    data = data .* corrph;
    %MRS_struct.fids.data = MRS_struct.fids.data .*repmat(conj(MRS_struct.fids.data(1,:))./abs(MRS_struct.fids.data(1,:)),[MRS_struct.p.npoints 1]);
    %Philips data appear to be phased already (ideal case)
    data = conj(data); %RE 110728 - empirical factor to scale 'like GE'
    data = reshape(data,[size(data,1) size(data,2)*size(data,3)]);
    % Depending on the ordering of OFF and ON, the minus sign here may have
    % to be removed.... not sure how to automate/fix... raee 4/9/12
    data = data .*repmat(conj(data(1,:))./abs(data(1,:)),[npoints 1]);
    
    % if nargin >2 %- ADH 20130701 don't need this if because the Ref data is
    % in the same file
    data_water = conj(data_water); %RE 110728 - empirical factor to scale 'like GE'
    data_water = data_water .*repmat(conj(data_water(1,:))./abs(data_water(1,:)),[npoints 1]);
    % Apply spectral registration to water, use 3-end water scans
    %[dummy, MRS_struct] = Spectral_Registration(MRS_struct,1);
    tempwater = data_water;
    if(size(data_water,2)>3)
        data_water = squeeze(mean(data_water(:,3:end),2)).'; %adh 20130208 sdat is mean of .data - not sum
    else
        data_water = squeeze(mean(data_water(:,end),2)).'; %adh 20130208 sdat is mean of .data - not sum
    end
    
    Nwateravg = 1; %adh 20130208 sdat is mean of .data - not sum - Nwateravg now = to sdat
    % not sure if this is needed but will keep
    % in for now. - ADH 20130701
    % end
    
    
    
else  % ADH 20130708 this is the PhilipsRead_data with only 1 mix (so may be separate water files or no water reference
    % work out data header name
    sparname = [fname(1:(end-4)) 'list'];
    sparheader = textread(sparname, '%s');
    sparidx=find(ismember(sparheader, 'F-resolution')==1);
    npoints = str2num(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1);
    nrows = str2num(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1);
    Navg = str2num(sparheader{sparidx+2});
    Navg = Navg*nrows;
    %Need to determine the offset of the data - i.e. how many channels are
    %there...
    sparidx=find(ismember(sparheader, 'NOI')==1);
    coil_channels=size(sparidx,1)-2;
    sparidx=find(ismember(sparheader, 'STD')==1);
    ptr_offset=str2num(sparheader{sparidx(3)+20});
    
    %Need to skip rows associated with the '
    data = readraw_Gannet(fname, 'float', [2 npoints 1 Navg/nrows nrows], 'l',ptr_offset);
    
    %  Make data complex.
    data = squeeze(data(1,:,:,:,:)+ 1i*data(2,:,:,:,:));
    
    %%'Duplicated' code added in to handle .data water files
    % work out water data header name
    
    if nargin >2
        sparnameW = [fname_water(1:(end-4)) 'list'];
        sparheaderW = textread(sparnameW, '%s');
        sparidx=find(ismember(sparheaderW, 'number_of_signal_averages')==1);
        Nwateravg = str2num(sparheaderW{sparidx+2}); %Trial SDAT might be average not sum.
        
        
        %Need to skip rows associated with the '
        
        data_water = readraw_Gannet(fname_water, 'float', [2 npoints Nwateravg 1], 'l', ptr_offset);
        
        %  Make data complex.
        data_water = squeeze(data_water(1,:,:,:,:)+ 1i*data_water(2,:,:,:,:));
    end
    
    
    
    
    %undo time series phase cycling to match GE
    corrph = ones(size(data));
    for jj=1:size(data,3)
        corrph(:,:,jj) = corrph(:,:,jj) * (-1).^(jj+1);
    end
    
    data = data .* corrph;
    %MRS_struct.fids.data = MRS_struct.fids.data .*repmat(conj(MRS_struct.fids.data(1,:))./abs(MRS_struct.fids.data(1,:)),[MRS_struct.p.npoints 1]);
    %Philips data appear to be phased already (ideal case)
    data = conj(data); %RE 110728 - empirical factor to scale 'like GE'
    data = reshape(data,[size(data,1) size(data,2)*size(data,3)]);
    % Depending on the ordering of OFF and ON, the minus sign here may have
    % to be removed.... not sure how to automate/fix... raee 4/9/12
    data = data .*repmat(conj(data(1,:))./abs(data(1,:)),[npoints 1]);
    
    if nargin >2
        data_water = conj(data_water); %RE 110728 - empirical factor to scale 'like GE'
        data_water = data_water .*repmat(conj(data_water(1,:))./abs(data_water(1,:)),[npoints 1]);
        data_water = squeeze(mean(data_water,2)).'; %adh 20130208 sdat is mean of .data - not sum
        
        Nwateravg = 1; %adh 20130208 sdat is mean of .data - not sum - Nwateravg now = to sdat
    end
    
end

FullData=data;
WaterData=data_water;

end




