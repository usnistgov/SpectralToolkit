% A script to perform a fourier analysis on a synthesized FM signal

% Instantiate an AnalyticTS_class default instance
TS = AnalyticTS_class();
[~,~,Ps,~,~,~,Fa,Ka] = TS.getParamIndex();
SignalParams = TS.SignalParams;
SignalParams(Ps,:) = -90;    % creates a sine wave
SignalParams(Fa,:) = 2;      % 2 Hz modulation
SignalParams(Ka,:) = .1;    % +/- 5 Hz depth
TS = AnalyticTS_class('SignalParams',SignalParams);

% Instantiate a FourierSeries_class object
%FS = FourierSeries_class('TimeSeries',TS.Ts);

%Set up for FmFitter
DelayCorr = [0];
MagCorr = [1];
F0 = TS.F0;
AnalysisCycles = 6;
SampleRate = TS.SampleRate;

for i = 1:18
    Samples = real(TS.getWindow(i,AnalysisCycles));
    
    FmFitter_2sb( ...
        SignalParams,...
        DelayCorr,...
        MagCorr,...
        F0,...
        AnalysisCycles,...
        SampleRate,...
        Samples ...
        );
end



