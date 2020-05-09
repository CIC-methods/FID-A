function [Output, AdditionalOut] = op_ReconstructNonCartMRData(Output,B0,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         ?                     ...     
% -         ?                     ...     
% -         ?             ...     
%
% Output:
% -         ?                      ...     
% -         ?                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations



if(~exist('Settings','var'))
   Settings.Phaseroll_flag = true;
   Settings.DensComp_flag = true;   
    
end
if(~isfield(Settings,'Phaseroll_flag'))
   Settings.Phaseroll_flag = true;    
end
if(~isfield(Settings,'DensComp_flag'))
   Settings.DensComp_flag = true;    
end
if(~isfield(Settings,'DensComp'))
   Settings.DensComp = struct();    
end
if(~isfield(Settings,'ConjInBegin_flag'))
   Settings.ConjInBegin_flag = false;    
end
if(~isfield(Settings,'ConjAtEnd_flag'))
   Settings.ConjAtEnd_flag = true;    
end
if(~isfield(Settings,'Correct4SpatialB0_flag'))
   Settings.Correct4SpatialB0_flag = false;    
end
if(~isfield(Settings,'CircularSFTFoV_flag'))
   Settings.CircularSFTFoV_flag = false;    
end
if(~isfield(Settings,'DensCompAutoScale_flag'))
   Settings.DensCompAutoScale_flag = false;    
end


%% FOV SHIFTs

% NOT YET IMPLEMENTED

% THATS HOW LUKAS HINGERL DOES IT (CODE FROM LUKAS HINGERL):
% %Inplane FOV shift
% kspace_max=sqrt(kx(1,end).^2+ky(1,end).^2);
% yIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Sag; %this should be right!
% xIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Cor;
% zIPFOVshift=ReadInInfo.General.Ascconv.PosVOI_Tra;
% FOV=ReadInInfo.General.Ascconv.FoV_Phase;
% FOVz=ReadInInfo.General.Ascconv.FoV_Partition;
% 
% 
% for circles=1:nc
%     for samplepoint=1:ns
%         csi_k(:,circles,samplepoint,1,:)=csi_k(:,circles,samplepoint,1,:)*exp(1i*kx(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*xIPFOVshift)*exp(1i*ky(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*yIPFOVshift);
%     end
% end



%% Conj in Beginning

if(Settings.ConjInBegin_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
if(Settings.Phaseroll_flag)

    nTI = Output.RecoPar.nTempInt;
    vs = Output.RecoPar.vecSize;
    ns = Output.RecoPar.TrajPts;
    nc = Output.RecoPar.nAngInts;
    nrew = Output.RecoPar.RewPts;
    ncha = size(Output.Data,6);

    timeoffset = 0:(ns-1);
    timeoffset = repmat(transpose(timeoffset),[1 vs]);
    Freq = ((0:vs-1)/vs-0.5)/(nrew + ns)*nTI;
    Freq = repmat(Freq,[ns 1]);    

    % This comes from:
    % timeoffset = (0:(ns-1))*Output.RecoPar.ADC_Dt/10^6;
    % sBW = nTI/((nrew + ns)*Output.RecoPar.ADC_Dt/10^6);
    % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
    % Output.RecoPar.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
    % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
    % the same anyway
    
    phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    % Freq(:,2*end/3)
    phasecorr = myrepmat(phasecorr,size(Output.Data));
%     phasecorr = conj(phasecorr);
%     bla = Output.Data(:,:,:,:,1,:,:);
    Output.Data = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.Data,[],5),5),5),[],5);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.NoiseData,[],5),5),5),[],5);
    end
%     Output.Data = conj(phasecorr).*Output.Data;     % For correcting only one constant frequency (e.g. at 3ppm which is about the metabo region)

%     Output.Data(:,:,:,:,1,:,:) = bla;

    TiltTrajMat = reshape(phasecorr(:,:,1,1,:,1),[Output.RecoPar.TrajPts*Output.RecoPar.nAngInts Output.RecoPar.vecSize]);   
    
else
    
	TiltTrajMat = ones([Output.RecoPar.TrajPts*Output.RecoPar.nAngInts Output.RecoPar.vecSize]);

end

if(nargout > 1)
    AdditionalOut.TiltTrajMat = TiltTrajMat;
end



%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
SizeData_k = size(Output.Data); SizeData_k = cat(2,SizeData_k,ones([1 5-numel(SizeData_k)]));
Output.Data = reshape(Output.Data,[prod(SizeData_k(1:2)) prod(SizeData_k(3:end))]);   
if(isfield(Output,'NoiseData'))
    Output.NoiseData = reshape(Output.NoiseData,[prod(SizeData_k(1:2)) prod(SizeData_k(3:end))]);   
end


sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Output.RecoPar.DataSize(1)),transpose(Output.InTraj.GM(:,:)),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Output.RecoPar.DataSize(1:2)),[1 2],[1 1 1 Output.RecoPar.DataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end

%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(Settings.Correct4SpatialB0_flag)
    t   = (0:Output.RecoPar.TrajPts-1)*Output.RecoPar.ADC_dt/10^9;
    t = repmat(t,[1 1 Output.RecoPar.nAngInts]); t = t(:);
    CurB0 = imresize(B0.B0Map,Output.RecoPar.DataSize(1:2));    
    if(isfield(B0,'Mask'))
        Mask = imresize(MaskShrinkOrGrow(B0.Mask,2,0,1),Output.RecoPar.DataSize(1:2),'nearest');
        CurB0 = CurB0 .* Mask;
    end
    
    B0CorrMat_Spatial = exp(transpose(-2*pi*1i*CurB0(:)) .* t);

    % SpSpice.Spi2Cart.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;
    % SpSpice.Spi2Cart_NoTiltCorr.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;

    sft2_Oper = sft2_Oper .* B0CorrMat_Spatial;
    
end




%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    [Output,AdditionalOut.DCFPreG] = op_CalcAndApplyDensComp(Output,sft2_Oper,Settings.DensComp);
end






%% Apply sft2-Operator (Fourier-Transform from spiral k-space --> Cartesian image 

Output.Data = sft2_Oper' * Output.Data * size(Output.OutTraj.GM(:,:),2);
Output.Data = reshape(Output.Data,[Output.RecoPar.DataSize(1:2) SizeData_k(3:end)]);
if(isfield(Output,'NoiseData'))
    Output.NoiseData = sft2_Oper' * Output.NoiseData * size(Output.OutTraj.GM(:,:),2);
    Output.NoiseData = reshape(Output.NoiseData,[Output.RecoPar.DataSize(1:2) SizeData_k(3:end)]);
end

% Remove fov_overgrid
Output.RecoPar.DataSize(1:2) = Output.RecoPar.DataSize(1:2)/Output.RecoPar.fov_overgrid;
bla = ([size(Output.Data,1) size(Output.Data,2)] - Output.RecoPar.DataSize(1:2))/2+1;
Output.Data = Output.Data(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);
if(isfield(Output,'NoiseData'))
    Output.NoiseData = Output.NoiseData(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);
end

if(nargout > 1)
    AdditionalOut.sft2_Oper = sft2_Oper;
end


%% Perform Reconstruction in Slice and z-dimension

% For now just reshape them. We dont have slices or 3D-measurements for now...
Size = size(Output.Data);
Output.Data = reshape(Output.Data, [Size(1:2) prod(Size(3:4)) Size(5:end)]); 
if(isfield(Output,'NoiseData'))
    Output.NoiseData = reshape(Output.NoiseData, [Size(1:2) prod(Size(3:4)) Size(5:end)]); 
end


%% Conj at End

if(Settings.ConjAtEnd_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end


%% Flip left right

% Output.Data = flip(Output.Data,2);


%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);

