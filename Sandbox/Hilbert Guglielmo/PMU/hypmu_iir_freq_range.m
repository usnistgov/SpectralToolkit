clear
close all
clc

%% Code description
% this code investigates the performance of Hilbert Transform-based PMU in
% static test, non distorted test condition - frequency range
% in this case, we implement HT through an IIR filter approach (ord 61)

%% Acquisition parameters
Fs = 10e3; % sampling frequency, Hz
Ts = 1/Fs; % sampling step, s
Tt = 2; % acquisition time, s
t_axis = 0:Ts:(Tt-Ts); % time axis, s
Ns = length(t_axis); % number of samples

%% Signal parameters
A = 1; % amplitude, pu
f = 55; % frequency, Hz
p = deg2rad(0); % initial phase, rad
s = A*cos(2*pi*f.*t_axis + p); % cosine signal, pu

%% Noise generation
SNRdB = 80; % signal to noise ratio, dB
stdn = (1/sqrt(2))/(10^(SNRdB/20)); % corresponding noise std dev, pu
s = s + stdn*randn(1,numel(s)); % adding noise to signal, pu

%% PMU parameters
Tw = 0.04; % window length, s
Nw = round(Tw/Ts); % window sample length
Rr = Fs; % reporting rate, fps
Nr = round(1/(Rr*Ts)); % reporting period sample length
r_sample = 1:Nr:(Ns-Nw); % first indices of considered windows
r_axis = 0:(1/Rr):(t_axis(end) - Tw); % time reference, s
% NB: we refer the estimated quantities to the window midpoint

%% HT filter implementation
Nf  = 61; % filter order
Tw = 50; % transition width
h = fdesign.hilbert('N,TW', Nf, Tw, Fs); % definition of filter properties
Hd = design(h, 'ellip'); % design of filter through elliptic technique
s_hat = filter(Hd, s); % signal hilbert transform, pu
kmax = 2; % expected bin index of DFT maximum

%% Reference values
ft = f*ones(1, length(r_sample)); % frequency reference values, Hz
At = A*ones(1, length(r_sample)); % amplitude reference values, pu
pht = wrapToPi(2*pi*f.*r_axis + p); % initial phase reference values, rad
fhe = zeros(1, length(r_sample)); % initialize the frequency estimates, Hz
Ahe = zeros(1, length(r_sample)); % initialize the amplitude estimates, pu
phe = zeros(1, length(r_sample)); % initialize the phase estimates, rad

%% Estimation cycle
for i = 10001:length(r_sample)
    x_hat = s_hat(r_sample(i):(r_sample(i)+Nw-1)); % considered window, pu
%     x_hat = s((r_sample(i)):(r_sample(i)+Nw-1)); % considered window, pu
    
    cd biblio
    
    Y_hat = fftshift(fft(x_hat))/Nw; % DFT of analytic signal, pu
    Y_win = Y_hat(1+Nw/2+(-1:1:8)); 
    Y_hat = Y_hat((1+Nw/2):end); % only positive coefficients
    for ii = 1:(length(Y_win) - 2)
        Y_han(ii) = - 0.25*Y_win(ii) + 0.5*Y_win(ii+1) - 0.25*Y_win(ii+2);
    end
    [f3h(i), A3h(i), p3h(i)] = IpDFT(Y_han, Nw, kmax, Fs, '3p', 'hann'); 
    % 3 points - ipDFT, hanning
    cd ..
    
end

%% Estimation error
% evaluation of total vector error, frequency, amplitude and phase errors
[tve3h, fe3h, ae3h, pe3h] = ERR(At, ft, pht, 4*A3h, f3h, p3h); % 3pts hann

note = [max(tve3h(10001:end)) mean(tve3h(10001:end)) ...
    std(tve3h(10001:end)) max(fe3h(10001:end)) mean(fe3h(10001:end)) ...
    std(fe3h(10001:end))]; % statistical properties of TVEs and FEs
   
%% Saving results
ris3h = [tve3h, fe3h, ae3h, pe3h];

% eval(['save ris_fr_f' num2str(f*10) '.mat ris2r ris3r ris2h ris3h'])