function [] = FmFitter( ...
    SignalParams,...
    DelayCorr,...
    MagCorr,...
    F0,...
    AnalysisCycles,...
    SampleRate,...
    Samples ...
    )
% Citation: G. N. Stenbakken, "Calculating combined amplitude and phase
% modulated power signal parameters," 2011 IEEE Power and Energy Society
% General Meeting, Detroit, MI, USA, 2011, pp. 1-7, doi: 10.1109/PES.2011.6039442.
%

[~,Fin,~,~,~,~,Fa,Ka,Fx,Kx,~,~,~,~,~] = getParamVals(SignalParams);

% signal parameters (used for initial estimates
fm = Fa(1); Km = Ka(1);     % Assume  FM
if Km == 0
    fm = Fx(1); Km = Kx(1); % or irt is AM
end

Freqs = Fin;    % fundamental frequency
[NSamples,NPhases] = size(Samples);
Freq = Freqs(1);
ROCOFs = zeros(1,NPhases);


wF = 2*pi*Freq;     % For large Km, this might not be a good intial guess.
% try using Hilbert for a better initial guess
% h=hilbert(Samples);
% hPhi = angle(h);
% hF = gradient(unwrap(hPhi));
% idx = floor(length(hF)/2);
% wF = mean(hF(idx-2:idx+2))*SampleRate;
% disp(wF/(2*pi))

wm = 2*pi*fm;

dt = 1/SampleRate;
tn = (-(NSamples/2-(1/2)):NSamples/2-(1/2))*dt;
MaxIter = 40;
FitCrit = 1e-7;   %dFm min


% pre-allocate results
Ain = zeros(1,NPhases);
Theta = zeros(1,NPhases);
iterations = zeros(1,NPhases);

% method uses fundamental and 4 sidebands

% this will only calculate one phase
    Af = Km; phif = 0;
    Freqs = 0; ROCOFs = 0;
    
    for k = 1:MaxIter
       Afm = abs(Af); phif = angle(Af);
       
       % model matrix
       H = [cos(wF*tn + Afm*sin(wm*tn + phif))', ...
           sin(wF*tn + Afm*sin(wm*tn + phif))', ...
           cos((wF-wm)*tn + Afm*sin(wm*tn + phif))', ...
           sin((wF-wm)*tn + Afm*sin(wm*tn + phif))', ...
           cos((wF+wm)*tn + Afm*sin(wm*tn + phif))', ...
           sin((wF+wm)*tn + Afm*sin(wm*tn + phif))', ...
           ones(1,NSamples)'];
       
       S = (H'*H)\(H'*Samples);
  
       %-------------- DEBUG -----------------------------------
       % plot the model, original and the residual
       fit = sum(S.*H',1);
       figure(1)
       subplot(2,1,1)
       plot(tn,Samples,'b',tn,fit','r');

       subplot(2,1,2)
       plot(tn,(Samples-fit')) 
       
       %--------------
       % also plot an FFT of the original and of the fit
       ts = timeseries(Samples);
       ts = setuniformtime(ts,'Interval',1/SampleRate);
       FT = FourierSeries_class('Timeseries',ts);       
       figure(2)
       subplot(2,1,1)
       %FT.plot('SingleSided','yScale','log')
       FT.plot('SingleSided')
       xlim([0,100])
       
       ts = timeseries(fit');
       ts = setuniformtime(ts,'Interval',1/SampleRate);
       FT = FourierSeries_class('Timeseries',ts);
       subplot(2,1,2)
       %FT.plot('SingleSided','yScale','log')
       FT.plot('SingleSided')
       xlim([0,100])

       pause
       %---------------DEBUG---------------------------------------
      
       AF = sqrt(S(1)^2 + S(2)^2); phiF = atan2(-S(2),S(1));
       AL = sqrt(S(3)^2 + S(4)^2); phiL = atan2(-S(4),S(3));
       AU = sqrt(S(5)^2 + S(6)^2); phiU = atan2(-S(6),S(5));
       
       %LV reverse engineering - ModPhaseFit.vi
       mi = (AU*cos(phiU - phiF) - AL*cos(phiL-phiF)) + 1i*(AU*sin(phiU - phiF) - AL*sin(phiL-phiF));
       su = (AU*cos(phiU - phiF) + AL*cos(phiL-phiF)) + 1i*(AU*sin(phiU - phiF) + AL*sin(phiL-phiF));
       fcos = abs(mi)*cos(angle(mi))/AF;
       fsin = abs(su)*sin(angle(su))/AF;
       acos = abs(su)*cos(angle(su))/AF;
       asin = abs(mi)*sin(angle(mi))/AF;
       dFm = sqrt(fcos^2 + fsin^2);
       ma = sqrt(asin^2 + acos^2);
       phia = atan2(asin,acos);
       phif = atan2(fsin,fcos);
       Af = (Af + dFm*exp(1i*phif));
       if abs(dFm) < FitCrit
           break
       end
       
       
    end




%--------------------------------------------------------------------------
% local function to get the values of SignalParams
    function [varargout] = getParamVals(signalparams)
        for i = 1:nargout
            varargout{i}=signalparams(i,:);
        end
    end
end

