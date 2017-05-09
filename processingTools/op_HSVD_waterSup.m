%op_HSVD_waterSup.m
%Jay Hennessy, McGill University 2017.
%
% USAGE:
% out=op_HSVD_waterSup( out, K, M, plot_bool);
% 
% DESCRIPTION:
% This function removes the water signal from MRS data using HSVD method 
% described by H. BARKHUIJSEN et al. 1987.
% 
% 
% INPUTS:
% out   = MRS data structure used by FID-a toolkit. Data should be
%         pre-processed, for example by out = run_pressproc(filename)
%
% K     = The number of frequency components in the data model. (default is
%         18).
%
% M     = M is the integer number of columns in the henkel matrix. Note: L
%         is the number of rows and L+M=N where N is the number of data
%         points. For best results 0.5<=L/M<=2. (default M= 2048, with 
%         N=4096 assumed)
%
% plot_bool = if 1, water fit is plotted (default =1)
%
%
%

function [ out, Eeigs ] = op_HSVD_waterSup( in, K, M, plot_bool)

% set default values ( intended for seimens data with 4096 data points)
if nargin<4
    plot_bool =1;
    if nargin<3
            M = 3000;   
        if nargin<2
            K=20;
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
fid_temp = exp((-alpha + (1i*w))*t);

% do a least square fit of your model guess to the data
phamp = fid_temp'\fid';

% convert the eigenvalues to phase and amplitude of the model
[ph,amp] = cart2pol(real(phamp), imag(phamp));

% use alpitude, phase, damping factor and frequency to model data
fid_components = (amp.*exp(1i*ph))'*exp((-alpha + (1i*w))*t);
fid_model = sum(fid_components,1);
spec_model=fftshift(ifft(fid_model',[],1),1);


% calculate residual error of fit
er = sum(abs(spec-spec_model).^2)/length(spec);


% Model the water signal 
water = find(freqs<5); % find the frequencies components associated with water
fid_water = (amp(water).*exp(1i*ph(water)))'*exp((-alpha(water) + (1i*w(water)))*t);

% remove water signal from the data
fid_ws = fid-fid_water;
spec_ws=fftshift(ifft(fid_ws',[],1),1);


% plot the results
if plot_bool ==1
    figure;
    subplot(211)
    plot(ppm,spec);
    hold
    %plot(ppm,spec_ws, 'color','red');
    plot(ppm, -(spec_ws-spec), 'color','green');
    legend('Data','Water');
    set(gca,'XDir','reverse');
    xlabel('Freq');
    title('Original Data Spectrum + Water Estimate');
    
    subplot(212)
    plot(ppm,spec_ws);
    xlim([.2 6])
    title('Water Suppressed Spectrum');
    set(gca,'XDir','reverse');
    xlabel('Freq');
end

% set output
out = in;
out.fids = fid_ws';
out.specs = spec_ws;
out.watersupp.damp = alpha(water);
out.watersupp.freq = w(water);
out.watersupp.phase = ph(water);
out.watersupp.amp = amp(water);
out.watersupp.freq = freqs(water);
out.watersupp.freq_all = freqs;
out.watersupp.damp_all = alpha;
out.watersupp.k = K;
out.watersupp.residual_error = er;



end

