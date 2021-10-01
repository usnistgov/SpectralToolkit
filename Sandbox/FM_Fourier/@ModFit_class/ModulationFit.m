function [Synx,Freqs,ROCOFs, iterations] = ModulationFit (obj)

[~,Fin,~,~,~,~,Fa,Ka,Fx,Kx] = obj.TS.getParamVals(obj.SignalParams);

% Combined modulation assumes the same modulation frequency
if Ka(1) == 0
    Fm = Fx; Km = Kx; %modulation frequency
else
    Fm = Fa;Km = Ka; %modulation frequency, amplitude
end
dt = 1/obj.Fs;

if Km <= 0.375
    [Synx,Freqs,ROCOFs, iterations] = obj.ModFit2Sb(Fin,Fm,Km,real(obj.Window.Data),dt,obj.MagCorr,obj.DlyCorr);
else
    [Synx,Freqs,ROCOFs, iterations] = obj.ModFitNSb(Fin,Fm,Km,real(obj.Window.Data),dt,obj.MagCorr,obj.DlyCorr);
end
end
