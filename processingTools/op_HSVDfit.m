%op_HSVDfit.m
%Jay Hennessy, McGill University 2017.
%Jamie Near, Sunnybrook Research Institute 2025.
%
% USAGE:
% [model, resid, K, wppm, amp, alpha, ph]=op_HSVDfit(in,ppmlim,Kinit,M,plot_bool);
% 
% DESCRIPTION:
% This function fits an MRS spectrum to a series of lorentzian peaks using 
% the HSVD method described by H. BARKHUIJSEN et al. 1987.
% 
% INPUTS:
% in        = MRS data structure used by FID-a toolkit. Data should be
%             pre-processed, for example by out = run_pressproc(filename)
% ppmlim    = This is the frequency range of the spectrum to be fitted, in
%             ppm. (default = [0.2 4.2])
% Kinit     = The number of frequency components in the data model This parameter
%             might have to be played with. (default is 20).
% M         = M is the integer number of columns in the henkel matrix. Note: L
%             is the number of rows and L+M=N where N is the number of data
%             points. For best results 0.5<=L/M<=2. (default M= .75*length.
% plot_bool = if 1, water fit is plotted (default =1)
%
% OUTPUTS:
% model     = The model spectrum resulting from HSVD fit
% resid     = The residual (differnece betwen input and model)
% K         = The number of frequency components used to fit the data.
% ppms      = The frequencies of the components found [ppm]
% amp       = The amplitudes of the components found [arb units]
% alpha     = The damping factors of the components found [radians/s]
% ph        = The phases of the components found [radians]

function [model, resid, K, ppms, amp, alpha, ph] = op_HSVDfit(in,ppmlim,Kinit,M,plot_bool)

% set default values ( intended for seimens data with 4096 data points)
if nargin<5
    plot_bool =0;
    if nargin<4
        M = floor(in.sz(1)*.75);
        %M = 1500;  % good for data sets of 2048 points
        %M = 3000;  % good for data sets of 4096 points
        if nargin<3
            Kinit=30;   % this value might have to be played with
            if nargin<2
                ppmlim= [0.2 4.2];
            end
        end
    end
end

% set parameters
N =in.sz(1);
dt = in.dwelltime;
t=in.t(1:N);
ppm = in.ppm;

fid = in.fids(1:N)';
spec = in.specs;

% Make a hankel data matrix
H = hankel(fid(1:M),fid(M:end));

% calculate the svd
[U,S,V] = svd(H);

% initialize while loop
amp = 0;
count = 0;

% find the highest number of components
while sum(find(amp==0))>=1 || sum(isnan(amp) >= 1)
    
    %initialize K
    K = Kinit-count;
    
    % truncate the data
    Uk = U(:,1:K);
    Sk = S(1:K,1:K);
    Vk = V(:,1:K);
    
    % get he eigenvalues of the transform matrix
    Utk = Uk(2:end,:);
    Ubk = Uk(1:end-1,:);
    Eh = Utk\Ubk;
    
    Eeigs = eig(Eh');
    
    % convert eigenvalues to poles to get freq and damping factor
    [w_model,alpha_model] = cart2pol(real(Eeigs), imag(Eeigs));
    w = w_model/dt;
    freqs = w/(2*pi);
    alpha = (alpha_model-1)/dt;
    
    % make a model guess using only damping factor and freq
    fid_temp = exp((-alpha + (1i*w))*t);  % this is mostly zero except for 1 very small row

    % do a least square fit of your model guess to the data
    phamp = fid_temp'\fid';
   % NofNan = size(find(phamp==0))


    % convert the eigenvalues to phase and amplitude of the model
    [ph,amp] = cart2pol(real(phamp), imag(phamp));

    count = count+1;
    if K<2
        display('####### Could not find a suitable number of components ########');
        break;
    end
end

% use alpitude, phase, damping factor and frequency to model data
fid_components = (amp.*exp(1i*ph))'*exp((-alpha + (1i*w))*t);
fid_model = sum(fid_components,1);
spec_model=fftshift(ifft(fid_model',[],1),1);


% calculate residual error of fit
er = sum(abs(spec-spec_model).^2)/length(spec);


% Model the water signal 
ppms = -(freqs)/(in.txfrq/1000000)+4.65;
ppminrange = find(ppms>ppmlim(1) & ppms<ppmlim(2)); % find the frequencies components associated with water
% water = find(min(freqs));
% water = find(freqs>72 & freqs<76); % find the frequencies components associated with water
fid_model = (amp(ppminrange).*exp(1i*ph(ppminrange)))'*exp((-alpha(ppminrange) + (1i*w(ppminrange)))*t);

% remove water signal from the data
fid_resid = fid-fid_model;
spec_resid=fftshift(ifft(fid_resid',[],1),1);
spec_model=fftshift(ifft(fid_model',[],1),1);



% plot the results
if plot_bool ==1
    figure;
    subplot(211)
    plot(ppm,spec);
    hold
    %plot(ppm,spec_ws, 'color','red');
    plot(ppm, -(spec_resid-spec), 'color','green');
    legend('Data','Water');
    set(gca,'XDir','reverse');
    xlabel('Freq');
    title('Original Data Spectrum + Water Estimate');
    
    subplot(212)
    plot(ppm,spec_resid);
    xlim([.2 6])
    title('Water Suppressed Spectrum');
    set(gca,'XDir','reverse');
    xlabel('Freq');
end


if ~sum(amp(ppminrange))
    display('##########  The fit did not work. Try reducing the number of components K.  #############');
end

% write residual
resid = in;
resid.fids = fid_resid';
resid.specs = spec_resid;

%write model output
model=in;
model.fids=fid_model;
model.specs=spec_model;


end

