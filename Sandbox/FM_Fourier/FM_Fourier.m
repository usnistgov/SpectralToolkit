% A script to perform a fourier analysis on a synthesized FM signal

% % Instantiate an AnalyticTS_class default instance
% TS = AnalyticTS_class();
% [~,~,Ps,~,~,~,Fa,Ka] = TS.getParamIndex();
% SignalParams = TS.SignalParams;
% SignalParams(Ps,:) = -90;    % creates a sine wave
% SignalParams(Fa,:) = 0.9;      % 2 Hz modulation
% SignalParams(Ka,:) = .1;    % +/- 5 Hz depth
% TS = AnalyticTS_class('SignalParams',SignalParams);

% Instantiate a ModFit class
ModFit = ModFit_class('Duration',2);

% This first experiment will show the low modulation index using the 2 sideband method
% This experiment will loop through for 1 second of the waveform.  At the
% end it will show the errors
[~,~,Ps,~,~,~,Fa,Ka] = ModFit.TS.getParamIndex();
ModFit.SignalParams(Ps,:) = -90;    % creates a sine wave (rather than the default cosine wave)
ModFit.SignalParams(Fa,:) = 3.6;      % FM modulation frequency
ModFit.SignalParams(Ka,:) = 0.3;    % Modulation Index
ModFit.getTimeSeries();             % get the time series
ModFit.TS.Ts.Name = 'test50f0_0m0_0a3';


ModFit.runMod1Second(true);
ModFit.fig = ModFit.fig+1;          % increment the figure number for the next experiment

% plot the fourier of the time series
ts = timeseries(real(ModFit.TS.Ts.Data),'Name',ModFit.TS.Ts.Name);
ts = setuniformtime(ts,'Interval',1/ModFit.Fs);
ModFit.dispFourier(ts,ModFit.fig)
ModFit.fig = ModFit.fig+2;          % increment the figure number for the next experiment





