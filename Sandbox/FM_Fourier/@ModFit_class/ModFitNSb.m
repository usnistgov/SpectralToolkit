function [Synx,Freqs,ROCOFs, iterations] = ModFitNSb(Fin,Fm,Km,Samples,dt,MagCorr,DelayCorr)
% Fitter for modulated signals with N pairs of sidebands, up to 10 pairs.  The number of sidebands 
% is determined by the modulation index
% inputs:
%   Fin: Fundamental frequency
%   Fm: Modulation frequency
%   Km: Modulation Index
%   Samples: windows of sample data in rows, columns of phases
%   dt = Sampling period (1/Sample Rate) 
% 
% outputs:
% Synx: complex phasors at the center of the windows
% Freqs: instantanious frequencies at the center of the windows
% ROCOFs: instantanious rate of change of frequency at the center of the window
% iterations: number of iterations to converge (or MaxIter if did not converge)

[NSamples,NPhases] = size(Samples);

% angular rotation
wF = 2*pi*Fin;  % fundamental frequency
wm = 2*pi*Fm;   % modulation frequency

% time vector
tn = linspace(-NSamples/2,NSamples/2-1,NSamples)*dt;

MaxIter = 40;
FitCrit = 1e-7;   %dFm min

% pre allocate results
Freqs = zeros(1,NPhases);
ROCOFs = zeros(1,NPhases);
Ain = zeros(1,NPhases);
Theta = zeros(1,NPhases);
iterations = zeros(1,NPhases);

% nSb is the number of PAIRS of sidebands
% This decides how many sideband pairs to use
if(Km<0.5), nSb = 1;
elseif (Km<1),nSb = 2;
elseif(Km<=3),nSb = 3;
elseif(Km<5),nSb = 4;
elseif(Km<7),nSb= 5;
elseif(Km<9),nSb=8;
else,nSb=10;
end

% loop fits 2 sideband (1 lower and 1 upper) model
for p = 1:NPhases
    
    % preallocations
    Af = complex(Km,0);
    AF = zeros(1,nSb*2);
    phiF = zeros(1,nSb*2);
    H = zeros(length(tn),4*(nSb)+3);  % 4 columns per sideband pair plus 2 columns for fundamental 
    
    % NOTE: here, the modulation phase phif is assumed to be zero. However
    % it has been shown that if the initial guess for phif is not within 1
    % Rad of the actual, the fit will diverge.  We need to add a method of
    % estimating a better initial value for the modulaton phase
    phif = zeros(nSb+1,1);
    
    for k = 1:MaxIter
        Afm = abs(Af); phim = angle(Af);
        j = 1;
        sb = 1;
        % H is our model (H for Hypothesis)
        H(:,j) = cos(wF(sb)*tn + Afm*sin(wm*tn + phim))'; j=j+1;
        H(:,j) = sin(wF(sb)*tn + Afm*sin(wm*tn + phim))'; j=j+1;
        sb = sb+1;
        while sb <= (2*nSb)+1
            H(:,j) = cos((wF-(sb-1)*wm)*tn + Afm*sin(wm*tn + phim))'; j=j+1; % lower sideband
            H(:,j) = sin((wF-(sb-1)*wm)*tn + Afm*sin(wm*tn + phim))'; j=j+1;
            sb = sb+1;
            H(:,j) = cos((wF+(sb-2)*wm)*tn + Afm*sin(wm*tn + phim))'; j=j+1; % upper sideband
            H(:,j) = sin((wF+(sb-2)*wm)*tn + Afm*sin(wm*tn + phim))'; j=j+1;
            sb = sb+1;
        end
        H(:,j) = ones(length(tn),1);
        
        S = (H'*H)\(H'*Samples);
        
        j = 1;
        for sb = 1:2*nSb+1
            cF = complex(S(j),S(j+1)); 
            AF(sb) = abs(cF); 
            phiF(sb) = -angle(cF);
            j = j+2;        
        end
    
        mi = 0; 
        su = 0;
        a = 2;
        while a <= length(AF)
            mi = mi - AF(a)*cos(phiF(a)-phiF(1)) - 1i*(AF(a)*sin(phiF(a)-phiF(1))); % lower sideband
            su = su + AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1)));
            a = a+1;
            mi = mi + AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1))); % upper sideband
            su = su + AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1)));
            a = a+1;
        end
        
       fcos = abs(mi)*cos(angle(mi))/AF(1);
       fsin = abs(su)*sin(angle(su))/AF(1);
       acos = abs(su)*cos(angle(su))/AF(1);
       asin = abs(mi)*sin(angle(mi))/AF(1);
       phim = atan2(fsin,fcos);
       dFm = sqrt(fcos^2 + fsin^2);
       Af = (Af + dFm*exp(1i*phif(1)));
       if abs(dFm) < FitCrit
           break
       end
       Afm(:) = abs(Af); phif(:) = angle(Af);
                
    end

    ma = sqrt(asin^2 + acos^2);
    phia = atan2(asin,acos);
    AF = AF*(1+ma*cos(phia));
    phiF = phiF + abs(Af)*sin(angle(Af));
    Freqs(p) = Fin(p) + abs(Af)*Fm(p)*cos(angle(Af));
    ROCOFs(p) = -wm(p)^2*abs(Af)*sin(angle(Af))/(2*pi);
    Ain(p) = abs(AF(1))*MagCorr(p);
    Theta(p) = phiF(1) + DelayCorr(p)*1e-9*wF(p);
    iterations(p) = k;   
end
Synx = (Ain/sqrt(2).*exp(1i.*Theta)).';
end