close all
fig = 1; % control the figure numbering
%% ------------------------------------------------------------------------
%Set up the Time Series class
%(comment out this section if the TS has already been set up)
fileName = 'F:\Projects\NetZero\Data\132401_1000.csv';
TimeInfo.Units = 'seconds';
TimeInfo.Increment = 1/300000;
TimeInfo.Start = 0;
TS = TimeSeries_class('File',fileName,'TimeInfo',TimeInfo);
% plot(TS)
%--------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% % Example of interactive switching to a different column of data
% timeVals=TS.Ts.Time;
% TS.loadTs(timeVals)
% plot(TS)

% % Example of non interactively switching to a different column of data
% colName='EntubeZ_L1_N_Meter_Calc_V';
% timeVals=TS.Ts.Time;
% TS.loadTs(timeVals,colName)
% plot(TS)


%% -------------------------------------------------------------------------
% Conduct Fourier analysis and create a Fourier Analysis Class
colName='EntubeSE_L1_N_Cond_Calc_V';
timeVals=TS.Ts.Time;
TS.loadTs(timeVals,colName)
figure(fig); fig = fig+1;
plot(TS)

FS_V = FourierSeries_class('TimeSeries',TS.Ts,'Window','Rect');
FS_V = FourierSeries_class('TimeSeries',TS.Ts);
figure(fig); fig = fig+1;
plot(FS_V,'SingleSided','yscale','log','ylim',[1e-4,1e3])


colName='P41600_L1_N_Meter_Calc_A';
timeVals=TS.Ts.Time;
TS.loadTs(timeVals,colName)
figure(fig); fig = fig+1;
plot(TS)

FS_A = FourierSeries_class('TimeSeries',TS.Ts,'Window','Rect','TsIndex',[1,60000]);
figure(fig); fig = fig+1;
plot(FS_A,'SingleSided','yscale','log','xlim',[0, 3000],'ylim',[1e-3,1e2])
%plot(FS_A,'SingleSided','yscale','dB','xlim',[0, 3000],'ylim',[-60,40])
% 
