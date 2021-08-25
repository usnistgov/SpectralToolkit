function [] = FmFitter_nSb( ...
    SignalParams,...
    DelayCorr,...
    MagCorr,...
    F0,...
    AnalysisCycles,...
    SampleRate,...
    Samples ...
    )
%
% FM fitter creates a number of side bands appropriate for the size of the K = delta f / fm term
%
% Citation: G. N. Stenbakken, "Calculating combined amplitude and phase
% modulated power signal parameters," 2011 IEEE Power and Energy Society
% General Meeting, Detroit, MI, USA, 2011, pp. 1-7, doi: 10.1109/PES.2011.6039442.
%

[~,Fin,~,~,~,~,Fa,Ka,Fx,Kx,~,~,~,~,~] = getParamVals(SignalParams);

% ------------------- DEBUG ----------------------------------
       % plot an FFT of the original and of the fit
       ts = timeseries(Samples,'Name','Original');
       ts = setuniformtime(ts,'Interval',1/SampleRate);
       FT_orig = FourierSeries_class('Timeseries',ts,'Name','Original'); 
       
       figure(1)
       subplot(2,1,1)
       FT_orig.plot('SingleSided')
       xlim([0,100])          
       subplot(2,1,2)
       FT_orig.plot('Phase')
       xlim([0,100]) 
       
       % Plot a 3D polar view of the FFT
       figure(2)
       FT_orig.plot('Polar3')
       xlim([0,100])
       ylim([-1,1])
       zlim([-1,1])
       view([-57,27.8])
             
%--------------------------------------------------------------

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

% method uses fundamental and nSb sidebands
% This decides how many sideband pairs to use
%switch(Km)
    if(Km<0.5), nSb = 1;
    elseif (Km<1),nSb = 2;
    elseif(Km<=3),nSb = 3;
    elseif(Km<5),nSb = 4;
    elseif(Km<7),nSb= 5;
    elseif(Km<9),nSb=8;
    else,nSb=10;
    end
                
%end
% this will only calculate one phase loop here for more phases
    Freqs = 0; ROCOFs = 0;
    
    % preAllocate the H matrix
    H = zeros(length(tn),4*(nSb)+3);  % 4 columns per sideband pair plus 2 columns for fundamental 
    Afm = ones((nSb)+1,1)*Km;
    %phif = zeros((nSb)+1,1);    % A better guess at phif is needed
    phif = ones((nSb)+1,1)*DelayCorr;
    AF = zeros(nSb+1,1);
    phiF = zeros(nSb+1,1);
    Af = ones((nSb)+1,1)*Km;
    
    for k = 1:MaxIter
               
        j = 1;
        sb = 1;
        H(:,j) = cos(wF(sb)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
        H(:,j) = sin(wF(sb)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
        sb = sb+1;
        while sb <= nSb+1
            H(:,j) = cos((wF-(sb-1)*wm)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
            H(:,j) = sin((wF-(sb-1)*wm)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
            H(:,j) = cos((wF+(sb-1)*wm)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
            H(:,j) = sin((wF+(sb-1)*wm)*tn + Afm(sb)*sin(wm*tn + phif(sb)))'; j=j+1;
            sb = sb+1;
        end
            H(:,j) = ones(length(tn),1);
       
       S = (H'*H)\(H'*Samples);
  
%% -------------- DEBUG -----------------------------------
       % plot the model, original and the residual
       fit = sum(S.*H',1);
       figure(3)
       subplot(2,1,1)
       plot(tn,Samples,'b',tn,fit','r');

       subplot(2,1,2)
       plot(tn,(Samples-fit')) 
       
       %--------------
       ts = timeseries(fit','Name','Fitted Model');
       ts = setuniformtime(ts,'Interval',1/SampleRate);
       FT_fit = FourierSeries_class('Timeseries',ts);
       figure(4)       
       subplot(2,1,1)       
       FT_fit.plot('SingleSided')
       xlim([0,100])
       subplot(2,1,2)
       FT_fit.plot('Phase')
       xlim([0,100])        
       
       %pause
       %pause
%% ---------------DEBUG---------------------------------------
        j = 1;
        for sb = 1:2*nSb+1
            AF(sb) = sqrt(S(j)^2+S(j+1)^2);
            phiF(sb) = atan2(-S(j+1),S(j));
            j = j+2;
        end
        %-------------------- experiment, different AFm and phif for each band-----------------
%         dfm_rms = 0;        
%         a = 1;          
   
%         mi = - AF(a)*cos(-phiF(a)) - 1i*(AF(a)*sin(-phiF(a))); % fundamental
%         su = AF(a)*cos(-phiF(a)) + 1i*(AF(a)*sin(-phiF(a)));
%         fcos = abs(mi)*cos(angle(mi))/AF(1);
%         fsin = abs(su)*sin(angle(su))/AF(1);
%         acos = abs(su)*cos(angle(su))/AF(1);
%         asin = abs(mi)*sin(angle(mi))/AF(1);
%         dFm = sqrt(fcos^2 + fsin^2);
%         dFm_rms = sqrt((dfm_rms^2+dFm^2)/2);
%         phia = atan2(fsin,fcos);
%         Af(a) = (Af(a) + dFm*exp(1i*phia)) ;
%         Afm(a) = abs(Af(a)); phif(a)=angle(AF(a));
%         a = a+1;
%         j = a;
%         while a <= length(AF)    % for each sideband pair (not including the fundamental)
%             mi = - AF(a)*cos(phiF(a)-phiF(1)) - 1i*(AF(a)*sin(phiF(a)-phiF(1))); % lower sideband
%             su = AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1)));
%             a = a+1;
%             mi = mi + AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1))); % upper sideband
%             su = su + AF(a)*cos(phiF(a)-phiF(1)) + 1i*(AF(a)*sin(phiF(a)-phiF(1)));
%             a = a+1;
%             
%             fcos = abs(mi)*cos(angle(mi))/AF(1);
%             fsin = abs(su)*sin(angle(su))/AF(1);
%             acos = abs(su)*cos(angle(su))/AF(1);
%             asin = abs(mi)*sin(angle(mi))/AF(1);
%             dFm = sqrt(fcos^2 + fsin^2);
%             phia = atan2(fsin,fcos);
%             Af(j) = (Af(j) + dFm*exp(1i*phia)) ;
%             Afm(j) = abs(Af(j)); phif(j)=angle(AF(j));            
%             j = j+1;
%         end
%         if dFm_rms < FitCrit
%             break
%         end
%         fprintf('dFm_rms = %f\n',dFm_rms);
        %--------------- end of experiment ---------------------------------
        mi = 0; su = 0;
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
       dFm = sqrt(fcos^2 + fsin^2);
       %ma = sqrt(asin^2 + acos^2);
       %phia = atan2(asin,acos);
       phia = atan2(fsin,fcos);
       Af = (Af + dFm*exp(1i*phia));
       if abs(dFm) < FitCrit
           break
       end
       Afm(:) = abs(Af); phif(:) = angle(Af);
       
    end

    fprintf('Iter = %d, dFm = %1.4e\n',k,dFm);
    fprintf('Afm = %f, phif = %f\n',Afm(1),phif(1)*180/pi)
       
%     figure(3)
%     FT_orig.plot('Polar3')
%     xlim([0,100])
%     ylim([-1,1])
%     zlim([-1,1])
%     view([-57,27.8])

%--------------------------------------------------------------------------
% local function to get the values of SignalParams
    function [varargout] = getParamVals(signalparams)
        for i = 1:nargout
            varargout{i}=signalparams(i,:);
        end
    end
end

