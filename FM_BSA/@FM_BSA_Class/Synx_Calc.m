function [Synx,Freq,ROCOF] = Synx_Calc(obj,estParams)
% using the estimated parameters, calculate the synchrophasor, Frequency
% and ROCOF at the center of the sample window

n = ceil(obj.nSamples/2);
mod_Freq = 2*pi*obj.Fm(obj.phase)*obj.dT;  % prior knowledge
wf = 2*pi*obj.Fm(obj.phase);

Af = (estParams.Modulo_BSA(2)/sqrt(2))*obj.MagCorr(obj.phase);
Theta = estParams.Fcarr_BSA*n...
        + (estParams.dF_BSA/mod_Freq) * sin(mod_Freq*n + estParams.Phim_BSA)...
        + estParams.Phi_BSA(2);
Theta = Theta + obj.DlyCorr(obj.phase)*1e-9*wf;
Synx = Af*exp(1i*Theta);

Freq = estParams.Fcarr_BSA/(2*pi*obj.dT)...
       - (estParams.dF_BSA/mod_Freq)*(obj.Fm(obj.phase)*cos(wf*n+estParams.Phim_BSA));
ROCOF = (estParams.dF_BSA/mod_Freq)*obj.Fm(obj.phase)^2*sin(wf*n+estParams.Phim_BSA)*2*pi;

end