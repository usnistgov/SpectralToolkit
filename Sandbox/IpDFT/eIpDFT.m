function [f, A, ph, delta, eps] = eIpDFT(X, N, kmax, fs, n, f_old, ...
    A_old, ph_old, k, interpolator)
%% Enhanced interpolated DFT algorithm
% this function computes the interpolated DFT estimation around the peak
% index kmax, and minimizes the leakage coming from the negative frequency
% components, based on the following inputs:
% - INPUT: 
% X: DFT of the considered signal window with Hanning windowing
% N: Hanning window length
% kmax: index of the DFT peak (indexing start from 0), if unknown put -1
% fs: sampling rate, in Hz
% n: number of iterations for the minimization of leakage from negative
% frequency component
% f_old, A_old, ph_old: estimates of amplitude, frequency and phase from 
% the previous iteration (at first iteration set them = 0)
% k: index of DFT bins where to compute the negative frequency component
% typically, k = 0:(N - 1)
% interpolator: string defining the interpolation type: '2p' or '3p' for
% two points or three points interpolation, respectively
% - OUTPUT:
% f, A, phi: estimated frequency, amplitude and phase, in rad
% delta: fractional bin offset of the true frequency with respect to kmax
% eps: variable to detect which is the second highest bin around kmax, i.e.
% -1 and +1 corresponds to left and right side of kmax
% further details in:
% G. Frigo, A. Derviškadi?, P. A. Pegoraro, C. Muscas and M. Paolone,...
% "Harmonic Phasor Measurements in Real-World PMU-Based Acquisitions,"...
% 2019 IEEE International Instrumentation and Measurement Technology ...
% Conference (I2MTC), Auckland, New Zealand, 2019, pp. 1-6.

%% Positive frequency estimate
if f_old == 0 % at first iteration
    [f, A, ph, delta, eps] = IpDFT(X, N, kmax, fs, interpolator); % perform
    % traditional IpDFT
else % at any other iteration
    % use previous estimate of frequency, amplitude and phase
    f = f_old; A = A_old; ph = ph_old;
end

%% Leakage minimization
for i = 1:n % at each iteration
    Xn = spec_imag(-f, A, -ph, N, fs, k); % recover the negative frequency
    % component spectrum in the bins indicated in k
    Xp = X - Xn; % subtract the negative frequency contribution
    [f, A, ph, delta, eps] = IpDFT(Xp, N, kmax, fs, interpolator); % IpDFT
    % on the data with "less" leakage
end
end

function [DFT, k] = spec_imag(f, A, ph, N, fs, k)
%% Spectral image reconstruction algorithm
% this function recovers the negative frequency component interference in
% bins indicated in k, based on the parameters of the positive frequency
% component (frequency, amplitude, phase)
% - INPUT:
% f, A, ph: frequency, amplitude and phase of the signal to recover
% N: window length
% fs: sampling rate, in Hz
% k: axis of DFT bins, typically k = 0:(N - 1)
% - OUTPUT:
% DFT: DFT of the negative frequency component in the selected bins k
% k: corresponding bin indices 
% please notice this function can be applied on multi-tone signals, where
% input parameters f, A and ph are vectors with the same length

DFT = zeros(1,numel(k)); % initialization of the DFT vector

for i = 1:length(f) % for each frequency of interest
    DFT = DFT + 2*A(i)*HnW_k(k - f(i)/fs*N, N, ph(i)); % reconstruction of
    % the DFT associated to i-th component with Hanning window weighing
end
end

function Wk = HnW_k(k, N, theta)
%% Hanning window
% this function reproduces the effect of a Hanning window based on the bin
% of interest, the window length, and the fractional displacement between
% true frequency and the frequency of the DFT bins
Wk = -0.25*DK_k(k-1, N, theta) + 0.5*DK_k(k, N, theta) - ...
    0.25*DK_k(k+1, N, theta); % Hanning window effect as function of three
    % Dirichlet kernels (rectangular window in the frequency domain)

end

function Dk = DK_k(k, N, theta) 
%% Rectangular window
% this function reproduces the effect of a rectangular window based on the
% bin of interest, the window length, and the fractional displacement btw
% true frequency and the frequency of the DFT bins
Dk = 1/N*sin(pi*k)./(sin(pi*k/N));
for i = 1:length(k)
    if k(i) == 0
        Dk(i) = 1;
    end
end
Dk = Dk.*exp(-1i*pi*k*(N-1)/N).*exp(1i*theta);
end