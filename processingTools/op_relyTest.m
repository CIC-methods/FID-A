% op_relyTest.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out = op_relyTest(in);
% 
% DESCRIPTION:
% Perform reliability testing on a dataset with multiple averages, as
% described in Slotboom J et al, Meas Sci Technol 20 (2009).  This involves
% calculating the skewness (3rd standard moment) and kurtosis (4th standard
% moment) of any dataset along the "averages" dimension.  The output structure
% contains skewness and kurtosis vectors (as a function of t), as well as 
% the mean of the absolute values skewness and kurtosis for all timepoints 
% where the average FID SNR is greater than 2.  
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Structure containing skewness and kurtosis indices for
%          reliability testing.

function out=op_relyTest(in);

if in.flags.averaged || in.dims.averages==0 || in.averages<2
    error('ERROR:  Averaging has already been performed!  Aborting!');
end

%Find N_snr (the fid index where the SNR is drops below 2):

%First calculate the noise in the last 25% of the fid (this value is taken 
%in the individual averages and then we take the mean value of the noise 
%across all averages:
noise=std(real(in.fids(ceil(0.75*in.sz(in.dims.t):end),:)));
noise=mean(noise);

%Now divide signal by noise to get an estimate of SNR.  Here, to estimate
%the signal, we're using the average of all the averages.

signal=mean(real(in.fids),in.dims.averages);
SNR=signal./noise;
plot(in.t,abs(SNR),[0 max(in.t)],[2 2]);

%N_snr is then defined as the index of the last point in the fid whose SNR
%exceeded 2.  
N_snr=find(SNR>2);
N_snr=N_snr(end);

%Calculate the skewness and kurtosis throughout the FID
skewVect=skewness(real(in.fids),1,in.dims.averages);
kurtVect=kurtosis(real(in.fids),1,in.dims.averages)-3;

%Calculate the k_skewness and k_kurtosis values for all times up to N_snr.
k_skew=mean(abs(skewVect(1:N_snr)));
k_kurt=mean(abs(kurtVect(1:N_snr)));

%Calculate the k_skewness and k_kurtosis values for all times after N_snr.
k_skew_noise=mean(abs(skewVect(N_snr+1:end)));
k_kurt_noise=mean(abs(kurtVect(N_snr+1:end)));

%Find the variance in skewness (var_k_skew) and kurtosis (var_k_kurt):
var_k_skew=mean((abs(skewVect(N_snr+1:end))-k_skew_noise).^2);
var_k_kurt=mean((abs(kurtVect(N_snr+1:end))-k_kurt_noise).^2);

%Now find the standard deviation in skewness and kurtosis:
std_k_skew=sqrt(var_k_skew);
std_k_kurt=sqrt(var_k_kurt);

%Now calculate the intervals for k_skewness and k_kurtosis:
skewInterval=[k_skew-std_k_skew k_skew+std_k_skew];
kurtInterval=[k_kurt-std_k_kurt k_kurt+std_k_kurt];

%Now determine if the data is reliable.  According to Slotboom et al., the
%spectrum is reliable if the skewInterval contains the value 0.3272 and if
%the kurtosisInterval contains the value 0.6101.  Therefore:
if (skewInterval(1) <= 0.3272 && 0.3272 <= skewInterval(2)) && (kurtInterval(1) <= 0.6101 && 0.6101 <= kurtInterval(2))
    isReliable=true;
else
    isReliable=false;
end

out.skewVect=skewVect;
out.k_skew=k_skew;
out.k_skew_noise=k_skew_noise;
out.var_k_skew=var_k_skew;
out.std_k_skew=std_k_skew;
out.skewInterval=skewInterval;

out.kurtVect=kurtVect;
out.k_kurt=k_kurt;
out.k_kurt_noise=k_kurt_noise;
out.var_k_kurt=var_k_kurt;
out.std_k_kurt=std_k_kurt;
out.kurtInterval=kurtInterval;

out.isReliable=isReliable;

