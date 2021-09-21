function [] = FmFitter_Theory( ...
    SignalParams,...
    DelayCorr,...
    MagCorr,...
    F0,...
    AnalysisCycles,...
    SampleRate,...
    Samples ...
    )
% This is an attempt to use the model directly from the Stenbakken Paper

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

[~,Fin,~,~,~,~,Fa,Ka,Fx,Kx,~,~,~,~,~] = getParamVals(SignalParams);

w0 = 2*pi*Fin;
wm = 2*pi*Fa;

NSamples = length(Samples);
t = (-(NSamples/2-(1/2)):NSamples/2-(1/2))*dt;
MaxIter = 40;

% First off, try a simple 1 sideband pair model for Km < 0.2
H = [...
     cos(w0*t);sin(w0*t);...                    % Fundamental
     cos((w0-w)*t+phi);sin((w0-w)*t+phi);...    % 1st lower sideband
     cos((w0+w)*t+phi);sin(w0+w)*t+phi);...
    ]