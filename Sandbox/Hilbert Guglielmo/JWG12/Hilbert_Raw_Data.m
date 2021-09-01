clear 
close all
clc

%% Load acquired data
filename = 'SavedModWindows.mat'; % name of the mat file
foldername = 'Downloads\'; % folder with the mat file
path = 'C:\Users\gugli\'; % path to the folder
load([path, foldername, filename])

%% Raw data setup
% the raw data are organized in a struct P where the i-th entry has the 
% following fields:
% SignalParams: 6 columns (3 phase voltage, 3 phase current) with the
% nominal amplitude, frequency, phase, modulation parameters, etc
% DelayCorr: correction for the delay error introduced by the DAQ board
% MagCorr: correction for the magnitude error introduced by the DAQ board
% F0: nominal frequency (i.e. = reporting rate)
% AnalysisCycle: window length in nominal cycles
% SampleRate: sampling rate in Hz
% Samples: each row is a channel - first 3 voltage, second 3 current

num_win = length(P); % number of acquired windows
Fs = P(1).SampleRate; % sampling rate, Hz
Ts = 1/Fs; % sampling period, s
Ns = size(P(1).Samples, 2); % sample length
Tw = Ns*Ts; % window length, s
t_axis = 0:Ts:(Ts*(Ns - 1)); % time axis, s

%% Modulation parameters
f_start = P(1).SignalParams(2, 1); % initial frequency, Hz
mod_k = P(1).SignalParams(8, 1); % modulation depth, rad
mod_f = P(1).SignalParams(7, 1); % modulation frequency, Hz

%% Reporting rate parameters
Fr = P.F0; % reporting rate, Hz
Tr = 1/Fr; % reporting period, s
r_axis = (Tw/2):Tr:(Tw/2 + Tr*(num_win - 1)); % reporting time axis, s
f_true = f_start + mod_k*mod_f*sin(2*pi*mod_f*r_axis - pi); % true freq, Hz

%% Hilbert analysis
load low_pass_filter.mat

for i = 1:num_win
    x = P(i).Samples(1,:);
    xf = filtfilt(Hlp.Numerator, 1, x);
    z = hilbert(xf - mean(xf));
    instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
    f_est(i) = instfreq(round(Ns/2));
end

