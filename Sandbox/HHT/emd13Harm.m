% our test signal will have 13 harmoncs
Fs = 4800;
F0 = 60;
AnalysisCycles = 6;
SignalParams = zeros(4+(12*3)+1,1);
SignalParams(1,:) = 1;   % Xm
SignalParams(2,:) = F0;   % Fin
SignalParams(3,:) = 0;    % Ps
SignalParams(4,:) = -1;   % delimiter
mags = [2.0,5.0,1.0,6.0,0.5,5.0,0.5,1.5,0.5,3.5,0.5,3.0];
for i = 1:12
    SignalParams(5+((i-1)*3),:) = F0*(i+1);
    SignalParams(6+((i-1)*3),:) = 0;
    SignalParams(7+((i-1)*3),:) = mags(i)/100;
end
SignalParams(4+(12*3)+1,:) = -1;  % delimiter

% instantiate the AnalyticTS class with these parameters
testCase.TS = AnalyticTS_class(...
    'SignalParams',SignalParams,...
    'SampleRate', Fs,...
    'F0', F0...
    );
Y = testCase.TS.getWindow(0,AnalysisCycles);
Y = Y.setuniformtime('Interval',1/Fs);

% x = real(Y.data);
% [imf,residual] = Emd(x);

HHT = HilbertHuang_class('TimeSeries',Y,'EMD',true);
