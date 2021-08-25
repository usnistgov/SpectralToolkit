% A script to perform a fourier analysis on a synthesized FM signal

% Instantiate an AnalyticTS_class default instance
TS = AnalyticTS_class();
[~,~,Ps,~,~,~,Fa,Ka] = TS.getParamIndex();
SignalParams = TS.SignalParams;
SignalParams(Ps,:) = -90;    % creates a sine wave
SignalParams(Fa,:) = 2;      % 2 Hz modulation
SignalParams(Ka,:) = 2;    % +/- 5 Hz depth
TS = AnalyticTS_class('SignalParams',SignalParams);

% Instantiate a FourierSeries_class object
%FS = FourierSeries_class('TimeSeries',TS.Ts);

%Set up for FmFitter
DelayCorr = [0];
%DelayCorr = 2.5656;
%DelayCorr = -1.3090
MagCorr = [1];
F0 = TS.F0;
AnalysisCycles = 25;
SampleRate = TS.SampleRate;

i = 1;
%for i = 1:18
    fprintf('DelayCorr: %f\n',wrapTo180(DelayCorr*180/pi));
    Samples = real(TS.getWindow(i,AnalysisCycles));
    
    FmFitter_nSb( ...
        SignalParams,...
        DelayCorr,...
        MagCorr,...
        F0,...
        AnalysisCycles,...
        SampleRate,...
        Samples ...
        );
    DelayCorr = DelayCorr+0.2513;
%end



