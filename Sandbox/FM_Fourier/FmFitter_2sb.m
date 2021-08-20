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
%try using Hilbert for a better initial guess
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
        Afm1 = abs(Af); phif1 = angle(Af);
%       Afm1 = abs(Af1); phif1 = angle(Af1);
        Afm2 = abs(Af); phif2 = angle(Af);       
%       Afm2 = abs(Af2); phif2 = angle(Af2);
       
       % model matrix
%        H = [cos(wF*tn + Afm*sin(wm*tn + phif))', ...
%            sin(wF*tn + Afm*sin(wm*tn + phif))', ...
%            cos((wF-wm)*tn + Afm1*sin(wm*tn + phif1))', ...
%            sin((wF-wm)*tn + Afm1*sin(wm*tn + phif1))', ...
%            cos((wF+wm)*tn + Afm1*sin(wm*tn + phif1))', ...
%            sin((wF+wm)*tn + Afm1*sin(wm*tn + phif1))', ...
%            ones(1,NSamples)'];
       
       % model matrix
       H = [cos(wF*tn + Afm*sin(wm*tn + phif))', ...
           sin(wF*tn + Afm*sin(wm*tn + phif))', ...
           cos((wF-wm)*tn + Afm1*sin(wm*tn + phif1))', ...
           sin((wF-wm)*tn + Afm1*sin(wm*tn + phif1))', ...
           cos((wF-(2*wm))*tn + Afm2*sin(wm*tn+phif2))',...
           sin((wF-(2*wm))*tn + Afm2*sin(wm*tn+phif2))',... 
           cos((wF+wm)*tn + Afm1*sin(wm*tn + phif1))', ...
           sin((wF+wm)*tn + Afm1*sin(wm*tn + phif1))', ...
           cos((wF+(2*wm))*tn + Afm2*sin(wm*tn + phif2))', ...
           sin((wF+(2*wm))*tn + Afm2*sin(wm*tn + phif2))', ...           
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

       %pause
       %---------------DEBUG---------------------------------------
       
       AF = sqrt(S(1)^2 + S(2)^2); phiF = atan2(-S(2),S(1));
       AL1 = sqrt(S(3)^2 + S(4)^2); phiL1 = atan2(-S(4),S(3));
       AL2 = sqrt(S(5)^2 + S(6)^2); phiL2 = atan2(-S(6),S(5));
       %AL2 = 0; phiL2 = 0;
       AU1 = sqrt(S(7)^2 + S(8)^2); phiU1 = atan2(-S(8),S(7));
       %AU1 = sqrt(S(5)^2 + S(6)^2); phiU1 = atan2(-S(6),S(5));
       AU2 = sqrt(S(9)^2 + S(10)^2); phiU2 = atan2(-S(10),S(9));      
       %AU2 = 0; phiU2 = 0;
       
       %LV reverse engineering - ModPhaseFit.vi
       mi = (....
             (AU1*cos(phiU1 - phiF) + AU2*cos(phiU2 - phiF))...
           - (AL1*cos(phiL1-phiF) + AL2*cos(phiL2 - phiF))...
            )...
           + 1i*(...
                 (AU1*sin(phiU1 - phiF)+AU2*sin(phiU2 - phiF))...
                 - (AL1*sin(phiL1-phiF)+AL2*sin(phiL2-phiF))...
                 );
       su = (...
             (AU1*cos(phiU1 - phiF) + AU2*cos(phiU2 - phiF))...
             + (AL1*cos(phiL1-phiF) + AL2*cos(phiL2 - phiF))...
            )...
            + 1i*(...
                  (AU1*sin(phiU1 - phiF) + AU2*sin(phiU2 - phiF))...
                  + (AL1*sin(phiL1-phiF))+(AL2*sin(phiL2 - phiF))...
            );
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
    fprintf('Iter = %d, dFm = %1.4e\n',k,dFm);




%--------------------------------------------------------------------------
% local function to get the values of SignalParams
    function [varargout] = getParamVals(signalparams)
        for i = 1:nargout
            varargout{i}=signalparams(i,:);
        end
    end
end

