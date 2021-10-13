
% Hooke_Jeeves_Student3.m
%
% This is a Matlab port of C code written by Gregory Kyriazis (InMetro, Br).

% Citation:
%Kyriazis G. A., “Estimating parameters of complex modulated signals from
%prior information about their arbitrary waveform components,” IEEE Trans.
%Instrum. Meas., v. 62, no. 6, pp. 1681-1686, June 2013.

% Citation: Kyriazis G. A., “A Cartesian method to improve the results and
% save computation time in Bayesian signal analysis,” in Advanced
% Mathematical and Computational Tools in Metrology and Testing X (AMCTM
% X), Series on Advances in Mathematics for Applied Sciences, vol. 86, F.
% Pavese; W. Bremser; A.G. Chunovkina; N. Fischer; A.B. Forbes (eds.),
% World Scientific, 2015, pp. 229-240.


% The method appears to begin by performing an NLS fit of the modulation signal 
% based on prior information about the modulation:

HJS = FM_HJS_Class();   % Instantiates the class with default properties
HJS.configure()         
HJS.mod_Freq_NLS()
HJS.mod_Amp_NLS() 
HJS.plot('NLS')
HJS.Freq_BSA()
HJS.Ampl_BSA()
HJS.plot('BSA')