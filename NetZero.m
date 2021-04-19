fileName = 'F:\Projects\NetZero\Data\132401_1000.csv';
TimeInfo.Units = 'seconds';
TimeInfo.Increment = 1/300000;
TimeInfo.Start = 0;
TS = TimeSeries_class('File',fileName,'TimeInfo',TimeInfo);
plot(TS.Ts)