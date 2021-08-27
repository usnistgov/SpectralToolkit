function [f, A, ph, fi, Ai, phi, fiter, piter] = iIpDFT_GU(X, N, kmax, ...
    fs, window, Q, P0, Pi, interpolator, fittedmodel)

[f, A, ph, delta, eps] = eIpDFT_GU(X, N, kmax, fs, window, P0, ...
    interpolator, 0, 0, 0); % fundamental tone parameters estimation

fi = 0;
Ai = 0;
phi= 0;

fiter = zeros(1, Q); % initialize frequency values
piter = zeros(1, Q); % initialize phase values

for i=1:Q
    X0 = spec_imag_GU(f, A, ph, N, fs, window, -1:1:11); % reconstruction 
    % of spectral fundamental tone contribution, pu
    Xi = X - X0; % once removed only interfering tone is present, pu
    [fi, Ai, phi, deltai, epsi] = eIpDFT_GU(Xi, N, -1, fs, window, Pi, ...
        interpolator, fi, Ai, phi); % interfering tone parameters estimate
    Xi = spec_imag_GU(fi, Ai, phi, N, fs, window, -1:1:11); % recovered 
    % spectral interfering tone contribution
    X0 = X - Xi; % fundamental tone estimation
    [f, A, ph, delta, eps] = eIpDFT_GU(X0, N, kmax, fs, window, P0, ...
        interpolator, f, A, ph); % fundamental tone parameters estimation
    fiter(i) = f; % recursive estimation of fundamental frequency, Hz
    piter(i) = ph; % recursive estimation of initial phase, rad
end
% ph = wrapToPi(ph - feval(fittedmodel, f)); % phase compensation, rad
end

