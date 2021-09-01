clear 
close all
clc

%% Load acquired data
path = fileparts(mfilename('fullpath'));   % path to this script.
foldername = '\..\..\..\Test\Data\';
filename = '50f0_2m0_2a5.mat'; % name of the mat file
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
f_true = f_start - mod_k*mod_f*sin(2*pi*mod_f*r_axis + 0.1247); % true freq, Hz

%% Hilbert analysis
load low_pass_filter.mat
% windows definition
win_hann = transpose(hann(Ns));
win_bhar = transpose(blackmanharris(Ns));
win_kais = transpose(kaiser(Ns, 2.5));

figure()

% - Hanning
for i = 1:num_win
    x = P(i).Samples(1,:);
    xf = filtfilt(Hlp.Numerator, 1, x);
    xw = xf.*win_hann;
    z = hilbert(xw - mean(xw));
    instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
    %plot(instfreq),ylim([40,60]),title('Hanning'),drawnow,pause(0.1)
    f_est_hann(i) = instfreq(round(Ns/2));
end

% - blackmann-harris
for i = 1:num_win
    x = P(i).Samples(1,:);
    xf = filtfilt(Hlp.Numerator, 1, x);
    xw = xf.*win_bhar;
    z = hilbert(xw - mean(xw));
    instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
    %plot(instfreq),ylim([40,60]),title('Blackman-Harris'),drawnow,pause(0.1)
    f_est_bhar(i) = instfreq(round(Ns/2));
end

% - kaiser
for i = 1:num_win
    x = P(i).Samples(1,:);
    xf = filtfilt(Hlp.Numerator, 1, x);
    xw = xf.*win_kais;
    z = hilbert(xw - mean(xw));
    instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
    %plot(instfreq),ylim([40,60]),title('Kaiser'),drawnow,pause(0.1)
    f_est_kais(i) = instfreq(round(Ns/2));
end

%% Comparison plot
plot(r_axis(1:end-1), f_true(2:end), '-o')
hold on
plot(r_axis, f_est_hann, '-o')
plot(r_axis, f_est_bhar, '-o')
plot(r_axis, f_est_kais, '-o')
hold off
legend('true','hann','bl.-har.','kaiser')
xlabel('Time (s)'), ylabel('Frequency (Hz)')