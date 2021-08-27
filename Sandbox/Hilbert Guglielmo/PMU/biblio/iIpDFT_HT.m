function [f, A, ph, delta, eps, fi, Ai, phi, deltai, epsi] = iIpDFT_HT(X, N, kmax, fs, window, Q, P0, Pi, interpolator)
 
[f, A, ph, delta, eps] = eIpDFT(X, N, kmax, fs, window, P0, interpolator, 0, 0, 0); % fundamental tone parameters estimation

fi = 0;
Ai = 0;
phi= 0;
deltai = 0;
epsi = 0;

for i=1:Q
    X0 = spec_imag(f, A, ph, N, fs, window); % fundamental tone reconstruction
    Xi = X - X0; % interfering tone estimation
    
    [fi, Ai, phi, deltai, epsi] = eIpDFT(Xi, N, -1, fs, window, Pi, interpolator, fi, Ai, phi); % interfering tone parameters estimation
    Xi = spec_imag(fi, Ai, phi, N, fs, window); % interfering tone reconstruction
    X0 = X - Xi; % fundamental tone estimation
    [f, A, ph, delta, eps] = eIpDFT(X0, N, kmax, fs, window, P0, interpolator, f, A, ph); % fundamental tone parameters estimation
    
end
end

