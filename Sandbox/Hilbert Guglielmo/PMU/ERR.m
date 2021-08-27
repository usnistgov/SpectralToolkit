function [tve, fe, ae, pe] = ERR(AT, fhT, phT, AE, fhE, phE)
%% ERR
% this function enables us to compute total vector error, frequency error,
% amplitude error and phase error based on the true and estimated values of
% fundamental component parameters

load p_model10kHz.mat p_model10kHz % modelled phase response in the
% frequency range of interest

phE = phE - feval(p_model10kHz, fhE)'; % phase shift to observation
phasorT = AT.*exp(1i*phT); % true phasor, pu
phasorE = AE.*exp(1i*phE); % estimated phasor

tve = 100 * (abs(phasorE-phasorT))./(abs(phasorT)); % total vector error, %

fe = fhE - fhT; % frequency error, Hz

ae = 100 * (AE - AT)./AT; % amplitude error, pu

pe = unwrap(phE) - unwrap(phT); % phase error, rad
pe = wrapToPi(pe); % scaled in the range [-pi, +pi], rad

end

