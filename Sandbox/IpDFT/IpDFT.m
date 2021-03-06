function [f, A, ph, delta, eps] = IpDFT(X, N, kmax, fs, interpolator)
%% Interpolated DFT algorithm
% this function computes the interpolated DFT estimation around the peak
% index kmax, based on the following inputs:
% - INPUT: 
% X: DFT of the considered signal window with Hanning windowing
% N: Hanning window length
% kmax: index of the DFT peak (indexing start from 0), if unknown put -1
% fs: sampling rate, in Hz
% interpolator: string defining the interpolation type: '2p' or '3p' for
% two points or three points interpolation, respectively
% - OUTPUT:
% f, A, phi: estimated frequency, amplitude and phase, in rad
% delta: fractional bin offset of the true frequency with respect to kmax
% eps: variable to detect which is the second highest bin around kmax, i.e.
% -1 and +1 corresponds to left and right side of kmax
% further details are available in:
% G. Frigo, A. Derviškadi? and M. Paolone, "Reduced Leakage Synchro...
% phasor Estimation: Hilbert Transform Plus Interpolated DFT," in ...
% IEEE Transactions on Instrumentation and Measurement, vol. 68, no. 10,...
% pp. 3468-3483, Oct. 2019.

%% Peak defintion
df = fs/N; % DFT resolution, in Hz
if kmax == -1 % in case of unknown kmax
    [~,kmax] = max(abs(X)); % kmax defined as the bin index associated to 
        % DFT module peak
    kmax = kmax - 1; % subtraction to compensate Matlab indexing that 
        % starts from 1 and not 0
end

%% Surrounding bins
if kmax == 0 % DC case
    eps = -1;
    kmax = kmax + 1;
elseif kmax == numel(X) - 1 % Nyquist case
    eps = 1;
    kmax = kmax - 1;
elseif abs(X(kmax + 1 + 1)) >= abs(X(kmax + 1 - 1)) % right-side
    eps = 1; 
elseif abs(X(kmax + 1 + 1)) <= abs(X(kmax + 1 - 1)) % left-side
    eps = -1;
end

indexkmax = kmax + 1;
Xmax = X(indexkmax);

Xepsp = X(indexkmax+eps);
Xepsm = X(indexkmax-eps);

%% Interpolation
% estimation of delta
switch interpolator
    case '2p'
        if abs(Xmax) - abs(Xepsp) == 0
            delta = 0;
        else
            delta = eps*(2*abs(Xepsp) - abs(Xmax))/...
                (abs(Xmax) + abs(Xepsp));
        end
    case '3p'
        if - abs(Xepsm) + 2*abs(Xmax) - abs(Xepsp) == 0
            delta = 0;
        else
            delta = 2*eps*(abs(Xepsp) - abs(Xepsm))/...
                (abs(Xepsm) + 2*abs(Xmax) + abs(Xepsp));
        end
end

% estimation of amplitude
if delta == 0
    A = abs(Xmax)*abs(delta^2-1);
else
    A = abs(Xmax)*abs((delta*pi)/sin(delta*pi))*abs(delta^2-1);
end

% estimation of frequency
f = (kmax + delta)*df;

% estimation of phase
ph = wrapToPi(angle(Xmax) - pi*delta);
% wrapToPi might not be implemented in Matlab basic version, so please
% replace it with an equivalent function
end