clear 
close all
clc

%% Generation parameters
Fs = 100e3; % sampling rate, Hz
Ts = 1/Fs; % sampling period, s

%% Signal parameters
A = 1; % amplitude, pu
f = 50; % frequency, Hz
ph = 0; % phase, rad
% - modulation parameters
fm = 1; % modulation frequency, Hz
km = 25; % modulation depth

%% Observation interval
Tm = 1/fm; % at least one period of the modulation
Nm = round(Tm/Ts); % sample length
Fr = 50; % reporting rate, Hz (typ. 50 or 100 Hz)
Tr = 1/Fr; % reporting period, s
Nr = round(Tr/Ts); % sample length

%% Experiment duration
Tt = 10*Tm; % total time, s
Nt = round(Tt/Ts); % sample length
t_axis = 0:Ts:(Ts*(Nt - 1)); % time axis, s
i_axis = 1:Nr:(Nt - Nm); % axis of the starting indices of intervals
r_axis = Tm/2 + t_axis(i_axis); % reporting time axis, s

%% Signal and reference values
x = A*cos(2*pi*f*t_axis + km*cos(2*pi*fm*t_axis - pi)); % signal

f_ref = f - km*fm*sin(2*pi*fm*r_axis - pi); % frequency reference, Hz
Rf_ref = -km*(2*pi*fm^2)*cos(2*pi*fm*r_axis - pi); % ROCOF reference, Hz/s

figure, plotyy(r_axis, f_ref, r_axis, Rf_ref)

%% Hilbert analysis
for i = 1:length(i_axis)
    z = hilbert(x(i_axis(i):(i_axis(i) + Nm - 1)));
    instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
    f_est(i) = instfreq(round(Nm/2));
end

%% Error analysis
figure, subplot(211)
plot(r_axis, f_ref), hold on, plot(r_axis, f_est, '--r')
legend('reference','estimate'), title('(a)'), grid on
xlabel('Time (s)'), ylabel('Frequency (Hz)')
subplot(212)
stem(r_axis, (f_est - f_ref)*1e3)
legend('error'), title('(b)'), grid on
xlabel('Time (s)'), ylabel('Frequency Error (mHz)')