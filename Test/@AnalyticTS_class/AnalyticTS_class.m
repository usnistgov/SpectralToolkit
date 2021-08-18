classdef AnalyticTS_class
    %ANALYTICTS_CLASS creates a time series object containing one or more
    %channels of analytic (complex) signals.
    
    %The parameters of these signals match the PMU perfoemance standard.
    % This is different from the ArtificialTS class.  That class creates
    % real waveforms which are fourier combinations of pre sinusoidal
    % components. This class can create dynamic waveforms cush as AM or FM
    % modulated signals, frequency ramps, steps in amplitude or phase or
    % combinations of only two sinusoidal components one harmonic or one
    % interharmonic
    
% % signalparams  (Note that the labeling convention comes mostly from the standard)
%     Xm = signalparams(1,:)*sqrt(2);     % phase amplitude (given by the user in RMS)
%     Fin = signalparams(2,:);    % frequency (must be the same for all 6 channels or an error will be thrown
%     Ps = signalparams(3,:);     % phase 
%     Fh = signalparams(4,:);     % Frequency of the interfering signal
%     Ph = signalparams(5,:);     % Phase of the interfering signal
%     Kh = signalparams(6,:);     % index of the interfering signal    
%     Fa = signalparams(7,:);     % phase (angle) moduation frequency
%     Ka = signalparams(8,:);     % phase (angle) moduation index
%     Fx = signalparams(9,:);     % amplitude moduation frequency
%     Kx = signalparams(10,:);     % amplitude moduation index
%     Rf = signalparams(11,:);     % ROCOF
%     KaS = signalparams(12,:);   % phase (angle) step index
%     KxS = signalparams(13,:);   % magnitude step index
%     KfS = signalparams(14,:);   % frequency step index
%     KrS = signalparams(15,:);   % ROCOF step index (another way to create a frequency ramp)
%        
    
    properties
        Ts = timeseries()   % a single timeseries instance (see MATLAB timeseries object)
        SignalParams        % double array, rows of parameters and collumns of channels
        SampleRate
        T0
        Duration
        SettlingTime
        
    end
    
    %% =========================================================================
    %constructor
    methods
        function obj = AnalyticTS_class(varargin)
            %ANALYTICTS_CLASS Construct an instance of this class
            %   Useage: obj = AnalyticTS_class(<name,value>,<name,value>,...)
            % all arguments are optional
            
            defaultName = 'Time series';
            %defaultDescription = '';
            defaultT0 = 0;             % start time in seconds
            defaultSampleRate = 4800;  % samples per second
            defaultDuration = 1;       % duration of the data (in Ts.TimeInfo.Units)
            defaultUnits = 'seconds';
            defaultSettlingTime = 0;   % time of static sinusoid before any dynamics begin, also added to the end of the dynamics.
            defaultSignalParams = [1,50,0,0,0,0,0,0,0,0,0,0,0,0,0]';
            
            p = inputParser;
            
            % validation
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validScalarPosInt = @(x) isnumeric(x) && isscalar(x) && isinteger(x) && (x > 0);
            validateSignalParams = @(x) validateattributes(x,{'double'},{'nrows',15});
            
            addParameter(p,'Name',defaultName,@ischar);
            %addParameter(p,'Description',defaultDescription,@ischar);
            addParameter(p,'T0',defaultT0,validScalar);
            addParameter(p,'SampleRate',defaultSampleRate,validScalarPosNum);
            addParameter(p,'Duration',defaultDuration,validScalarPosNum);
            addParameter(p,'Units',defaultUnits,@ischar);
            addParameter(p,'SettlingTime',defaultSettlingTime,validScalarPosNum);
            addParameter(p,'SignalParams',defaultSignalParams,validateSignalParams);
            
            parse(p,varargin{:})
            obj.Ts.Name = p.Results.Name;
            %obj.Ts.Description = p.Results.Description;
            obj.T0 = p.Results.T0;
            obj.SampleRate = p.Results.SampleRate;
            obj.Duration = p.Results.Duration;
            obj.Ts.TimeInfo.Units = p.Results.Units;
            obj.SettlingTime = p.Results.SettlingTime;
            obj.SignalParams = p.Results.SignalParams;
            
            %obj = obj.makeUniformTime();
            obj = obj.AnalyticWaveforms();
            
        end
    end
    
    %% =========================================================================
    % public methods found in a class method file
    methods (Access = public)
        obj = AnalyticWaveforms(obj)
    end
    %% =========================================================================
    % static methods found in a class method file
    methods (Static)
        Size = SizeLcmPeriods(Freqs, Fs)
    end    
    
    %% =========================================================================
    % public methods
    methods (Access = public)
        
        function obj = makeUniformTime(obj)
            % Created the Ts.Time vector.
            totalTime = obj.Duration + 2*obj.SettlingTime;
            obj.Ts.Time = (0:(totalTime*obj.SampleRate)-1)/obj.SampleRate;
            obj.Ts = setuniformtime(obj.Ts,'StartTime',obj.T0,'Interval',1/obj.SampleRate);
        end
        
    end
    
    %% =========================================================================
    % static methods
    methods(Static)
        function [varargout] = getParamIndex(signalparams)
            for i = 1:nargout
                varargout{i}=signalparams(i,:);
            end
        end
        
%         %--------------------------------------------------------------------
%         % size = SizeLcmPeriod (Freqs, Fs)
%         % Determines the sample size of the least common multiple of cycles of a mix of
%         % multiple frequencies:
%         %
%         % Input:
%         %   [Freqs]: 1xN array of mixed frequencies
%         %   Fs: Sample rate
%         %
%         function Size = SizeLcmPeriods(Freqs, Fs)
%             
%             %--------------- ARG: The below requires the symbolic toolbox-----------
%             % % first determine the LCM of the periods.
%             % P = 1./sym(Freqs);      % creates a symbolic object of the periods
%             % l = lcm(P);
%             % ----------
%             % % We can work around needing the symbolic toolbox by converting to rational numbers
%             P = 1./Freqs;   % periods
%             P = rats(P);    % rational number of the periods
%             P = strsplit(P);    % separate the rational numbers into a cell array
%             P = P(2:end-1);     % first and last elements are whitespace
%             
%             % the get the numerators and denominators
%             n = zeros(numel(P),1); d = n;
%             for i = 1:numel(P)
%                 r = strsplit(P{i},'/');
%                 n(i) = str2double(r{1,1});
%                 d(i) = 1;
%                 if numel(r)==2; d(i) = str2double(r{1,2}); end;
%             end
%             nLcm = lcms(n);     % added an lcm set function to matlab path
%             dHcf = hcfs(d);     % added and hcf (gcd) set function to matlab path
%             l = nLcm/dHcf;      % lcm is the lcm of the numerators divided by the gcd of the denominators
%             %-----------------------------------------------------------------------
%             
%             
%             % This give us a "combined" frequency.  We now need to determine how many
%             % samples are needed to realize an integer number of cycles in that
%             % combined frequency.
%             
%             % First, check for fractional frequency and multiply by 10 until there is
%             % no fraction, if the number has too many decimal places, round it.
%             F = 1/double(l);    % combined "frequency"
%             %F = double(l);
%             Fmult = F;
%             for i = 1:4
%                 if rem(Fmult,1)==0;break;end
%                 Fmult = Fmult*10;
%             end
%             Fmult = round(Fmult);
%             
%             prime = factor(Fmult);
%             
%             i = 1;
%             cyc = 1;
%             while (i <= length(prime))
%                 s = cyc/F;
%                 n = Fs*s;
%                 if rem(n,1)==0;break;end;
%                 cyc = cyc*prime(i);
%                 i = i+1;
%                 
%             end
%             
%             s = cyc/F;
%             Size = Fs*s;
%             %end
%             
%             %--------------------------------------------------------------------------
%             % Josh (2021). Least Common Multiple Set
%             %(https://www.mathworks.com/matlabcentral/fileexchange/24670-least-common-multiple-set),
%             % MATLAB Central File Exchange. Retrieved January 9, 2021.
%             function output = lcms(numberArray)
%                 numberArray = reshape(numberArray, 1, []);
%                 %% prime factorization array
%                 for i = 1:size(numberArray,2)
%                     temp = factor(numberArray(i));
%                     
%                     for j = 1:size(temp,2)
%                         output(i,j) = temp(1,j);
%                     end
%                 end
%                 %% generate prime number list
%                 p = primes(max(max(output)));
%                 %% prepare list of occurences of each prime number
%                 q = zeros(size(p));
%                 %% generate the list of the maximum occurences of each prime number
%                 for i = 1:size(p,2)
%                     for j = 1:size(output,1)
%                         temp = length(find(output(j,:) == p(i)));
%                         if(temp > q(1,i))
%                             q(1,i) = temp;
%                         end
%                     end
%                 end
%                 %% the algorithm
%                 z = p.^q;
%                 output = 1;
%                 for i = 1:size(z,2)
%                     output = output*z(1,i);
%                 end
%                 
%                 
%             end
%             %--------------------------------------------------------------          
%             function out=hcfs(n1,varargin)
%                 % Computes Highest Common Factor of two or more numbers
%                 % Syntax:
%                 % [output variable]=hcf(number 1 [,number 2] [,..., number n])
%                 % For example:
%                 % hcf(6, 9)
%                 % hcf([6 9])
%                 % hcf(15, 30, 45)
%                 % combines multiple inputs into single variable
%                 %
%                 % Nitin (2021). hcf (https://www.mathworks.com/matlabcentral/fileexchange/53133-hcf),
%                 % MATLAB Central File Exchange. Retrieved January 9, 2021.
%                 
%                 if nargin > 1
%                     n1 = [n1 cell2mat(varargin(:))'];
%                 end
%                 n1=abs(n1);
%                 n1=round(n1);
%                 % Take care of zeros by eliminating them as hcf(0,a) = a
%                 n = n1(n1~=0);
%                 out=1;
%                 div=2;
%                 while div <= min(n)
%                     if sum(rem(n,div))==0
%                         out=out*div;
%                         n=n/div;
%                         div=2;
%                     else
%                         div=div+1;
%                     end
%                 end
%             end
%         end
    end
    %% =========================================================================
end

