% an experiment to perform a th order fit to a one-nominal-cycle unwrappe phase of a modulates signal

% Instantiate an AnalyticTS_class default instance
TS = AnalyticTS_class();
[~,Fin,Ps,~,~,~,Fa,Ka,~,~,~,~,~,~,~,Kn] = TS.getParamIndex();
SignalParams = TS.SignalParams;
SignalParams(Ps,:) = -90;    % creates a sine wave
SignalParams(Fa,:) = 2;      % 2 Hz modulation
SignalParams(Ka,:) = 2.5;    % +/- 5 Hz depth
SignalParams(Kn,:) = 0.00;

TS = AnalyticTS_class('SignalParams',SignalParams);

F0 = SignalParams(Fin,1);   % center frequency of the modulation
SampleRate = TS.SampleRate;
AnalysisCycles = 1;   %Analyse one nominal cycle
%AnalysisCycles = 3;   %Analyse three nominal cycles
%sps = SampleRate/F0;


% analyse the first nominal cycle
f = [];g = [];exp = [];
for i = 1:50
    Window = TS.getWindow(i-1,AnalysisCycles,'even');
    HHT = HilbertHuang_class('TimeSeries',Window,'EMD',true,'Hilbert',true,'Window','none');
    plot(HHT,'Hilbert')
    
    % for our purposes, the fundamental will be the IMF with the highest amplitude
    maxAmpl = 0; idxFund = 0;
    for i = 1:size(HHT.Hilbert,2)
        ampl = rms(abs(HHT.Hilbert{1,i}));
        if ampl > maxAmpl
            maxAmpl = ampl;
            idxFund = i;
        end
    end
    
    % to elimiate end effects, we will look at only the middle of the three analysis cycles
   %rH = real(HHT.Hilbert{1,i}(sps+1:(2*sps)));
   %iH = imag(HHT.Hilbert{1,i}(sps+1:(2*sps))); 
   
   % only one nominal cycle
   rH = real(HHT.Hilbert{1,i});
   iH = imag(HHT.Hilbert{1,i});   
   theta = unwrap(atan2(iH,rH));
   
   %theta = unwrap(HHT.Hilbert{3,i}(sps+1:2*sps));   % AnalysisCycles = 3 
   %theta = unwrap(HHT.Hilbert{3,i});   % AnalysisCycles = 1
    
    
    N = length(theta);
    %t = Window.Time(sps+1:2*sps)*SampleRate;   % AnalysisCycles = 3 
    t = Window.Time*SampleRate;    % AnalysisCycles = 1

    
    % model for the phase fitting used to remove noise
    %H = [ ones(length(t),1),t,t.^2,t.^3,t.^4,t.^5]; % fifth order model
    H = [ ones(length(t),1),t,t.^2,t.^3,t.^4];    % fourth order
    %H = [ ones(length(t),1),t,t.^2,t.^3];   % third order
    S = (H'*H)\(H'*theta);
    
    %F = S(1)+S(2)*t+S(3)*t.^2+S(4)*t.^3+S(5)*t.^4+S(6)*t.^5;  % fifth order model
    F = S(1)+S(2)*t+S(3)*t.^2+S(4)*t.^3+S(5)*t.^4;   % fourth order
    %F = S(1)+S(2)*t+S(3)*t.^2+S(4)*t.^3;  % third order
        
    % visualize the fit
    figure(100)
    subplot(3,1,1)
    plot(t/SampleRate,theta,'b',t/SampleRate,F,'r')
    title('hilbert and fitted theta')
    subplot(3,1,2)
    plot(t/SampleRate,theta-F,'r')
    title('residual')
    subplot(3,1,3)
    plot(t/SampleRate,gradient(F)*SampleRate/(2*pi))
    title('frequency (Hz)')
    refresh
    pause(1/60)
    
    f = vertcat(f,gradient(F)*SampleRate/(2*pi));
    g = vertcat(g,complex(rH,iH));
    % exp = vertcat(exp,Window.Data(sps+1:(2*sps)));
    exp = vertcat(exp,Window.Data);
    %g = vertcat(g,Window.Data);
    figure(101)
    plot(f)
    
    
    figure(102)
        plot(real(g),'b')
    hold on
    plot(imag(g),'r')
    plot(real(exp),'k')
    plot(imag(exp),'g')
    hold off
   
   end
