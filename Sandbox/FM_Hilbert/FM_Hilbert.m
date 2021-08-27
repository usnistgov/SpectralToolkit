% A script to perform a hilbert alanysis on a synthesized FM signal
%
% Instantiate an AnalyticTS_class default instance
TS = AnalyticTS_class();
[~,Fin,Ps,~,~,~,Fa,Ka] = TS.getParamIndex();
SignalParams = TS.SignalParams;
SignalParams(Ps,:) = -90;    % creates a sine wave
SignalParams(Fa,:) = 2;      % 2 Hz modulation
SignalParams(Ka,:) = 2;    % +/- 5 Hz depth
TS = AnalyticTS_class('SignalParams',SignalParams);

F0 = SignalParams(Fin,1);   % center frequency of the modulation
SampleRate = TS.SampleRate;
AnalysisCycles = ceil(F0/SignalParams(Ka,1));   %Analyse one modulation cycle

%set up for FmHilbert
%i = 1;
for i = 1:25    % loop for every nominal cycle in the modulation period

    ts = timeseries(real(TS.getWindow(i,AnalysisCycles,'even')),'Name',sprintf('Mod #%d',i));
    ts = setuniformtime(ts,'Interval',1/SampleRate,'StartTime',-((length(ts.Data)-1)/(2*SampleRate)));
    %HT = HilbertHuang_class('TimeSeries',ts,'EMD',true,'Hilbert',true);
    %HT.plot('IMFs')
    HT = HilbertHuang_class('IMFs',ts,'Hilbert',true);
    HT.plot('Hilbert')
    
end