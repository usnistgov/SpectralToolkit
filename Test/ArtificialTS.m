classdef ArtificialTS
    % Class for synthesizing a time series from sinusoids and noise
    % Name-Value Property Pairs  (default value)
    %   Name ()
    %   Description ()
    %   T0  (0) start time
    %   Extent (0) duration of TS
    %   nSamples (0) number of samples (must be integer type)
    %   Freqs () array of frequencies
    %   Amps () array of amplitudes
    %   Phases () array of phases
    %   NoiseUniformLow (0) scalor, uniform noise low value
    %   NoiseUniformHi (0) scalor, uniform noise high value
    %   NoiseGaussMean (0) scalor, normal noise mean value
    %   NoiseGaussSD (0) scalor, normal noise standard distribution
    %   rngState (-1) random number generator state.

    properties
        Name            % string name of the Time Series
        Description     % string describing the TS
        
        T0              % double starting time
        Extent          % double TS duration in seconds
        nSamples        % number of samples in the time series (Sample Rate = nSamples/Extent)
        
        Valid           % [1..nSinusoids] array of booleans (invalid if caller enters bad values)
        Freqs           % [1..nSinusoids] array of sinusoid frequencies
        Amps            % [1..nSinusoids] array of sinusoid amplitudes
        Phases          % [1..nSinusoids] array of phase angles (in radians)
        
        NoiseUniformLow % double Uniform distribution noise lowest value
        NoiseUniformHi  % double Uniform distribution noise highest value
        NoiseGaussMean  % double noise mean value
        NoiseGaussSD    % double noise standard deviation
        RngState = -1   % state of the random number generator.  -1 means no initial state
        noiseTs         % a time series containing only the additive noise.
        
        Ts              % Time Series
        time            % Time vector
    end
    
    methods
        %Constructor
        function ArtTS = ArtificialTS(varargin)
            
            
            defaultName = 'Time series';
            defaultDescription = '';
            defaultT0 = 0;
            defaultExtent = 5;
            defaultnSamples = int32(300);
            defaultFreqs = [1];
            defaultAmps = [1];
            defaultPhases = [0];
            defaultNoiseUniformLow = 0;
            defaultNoiseUniformHi = 0;
            defaultNoiseGaussMean = 0;
            defaultNoiseGaussSD = 0;
            defaultRngState = -1;
            
            p = inputParser;
            
            % validation
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validScalarPosInt = @(x) isnumeric(x) && isscalar(x) && isinteger(x) && (x > 0);
            
            
            addParameter(p,'Name',defaultName,@ischar)
            addParameter(p,'Description',defaultDescription,@ischar)
            addOptional(p,'T0',defaultT0,validScalar);
            addOptional(p,'Extent',defaultExtent,validScalarPosNum);
            addOptional(p,'nSamples',defaultnSamples,validScalarPosInt);
            addOptional(p,'Freqs',defaultFreqs,@isvector);
            addOptional(p,'Amps',defaultAmps,@isvector);
            addOptional(p,'Phases',defaultPhases,@isvector);
            addOptional(p,'NoiseUniformLow',defaultNoiseUniformLow,validScalar);
            addOptional(p,'NoiseUniformHi',defaultNoiseUniformHi,validScalar);
            addOptional(p,'NoiseGaussMean',defaultNoiseGaussMean,validScalar);
            addOptional(p,'NoiseGaussSD',defaultNoiseGaussSD,validScalar);
            addOptional(p,'RngState',defaultRngState);
            
            
            parse(p,varargin{:})
            ArtTS.Name = p.Results.Name;
            ArtTS.Description = p.Results.Description;
            ArtTS.T0 = p.Results.T0;
            ArtTS.Extent = p.Results.Extent;
            ArtTS.nSamples = p.Results.nSamples;
            ArtTS.Freqs = p.Results.Freqs;
            ArtTS.Amps = p.Results.Amps;
            ArtTS.Phases = p.Results.Phases;
            ArtTS.NoiseUniformLow = p.Results.NoiseUniformLow;
            ArtTS.NoiseUniformHi = p.Results.NoiseUniformHi;
            ArtTS.NoiseGaussMean = p.Results.NoiseGaussMean;
            ArtTS.NoiseGaussSD = p.Results.NoiseGaussSD;
            ArtTS.RngState = p.Results.RngState;
            
            ArtTS.Valid = ArtTS.checkValid;
            ArtTS.time = ArtTS.makeTime.time;    % create the time vector
            ArtTS.Ts = ArtTS.makeTS.Ts;      % create the time series
        end
        
       
        function Valid = checkValid(ArtTS)
            Valid = false;
            if length(ArtTS.Freqs) < 1
                error ('ArtTS.Freqs must have 1 or more frequencies')
            end
%             if length(ArtTS.Freqs) ~= length(ArtTS.Amps)
%                 error ('ArtTS.Freqs and ArtTS.Amps must be the same size')
%             end
%             if length(ArtTS.Freqs) ~= length(ArtTS.Phases)
%                 error ('ArtTS.Freqs and ArtTS.Phases must be the same size')
%             end
            if ArtTS.nSamples < 1
                error ('ArtTS.nSamples must be grater than 0')
            end
            if ArtTS.Extent <= 0
                error ('ArtTS.Extent must be greater than 0')
            end
            Valid = true;
        end
        
        function ArtTS = makeTime(ArtTS) % Creates a time vector from T0 to Extent+ 1 period with Extent/nSamples period
 
            ArtTS.time = (ArtTS.T0 : (ArtTS.Extent)/(ArtTS.nSamples):ArtTS.Extent-ArtTS.Extent/ArtTS.nSamples);
        end
                        
        function [ArtTS] = makeTS(ArtTS) % Creates an artificial time series
            ArtTS.Valid = ArtTS.checkValid;
            if length(ArtTS.time) < 1
                ArtTS.time = ArtTS.makeTime;
            end
            x = zeros(1,length(ArtTS.time));
            for i = 1:length(ArtTS.Freqs)
                x = x + ArtTS.Amps(i)*cos(2*pi*ArtTS.Freqs(i)*ArtTS.time - ArtTS.Phases(i));
            end
            
            % add noise
            bNoiseOn = false;
            NoiseUniformRange = ArtTS.NoiseUniformHi - ArtTS.NoiseUniformLow;
            if NoiseUniformRange ~= 0, bNoiseOn = true, end
            if ArtTS.NoiseGaussSD > 0
                bNoiseOn = true;
            end
            if bNoiseOn
                if isstruct(ArtTS.RngState), rng(ArtTS.RngState); end
                ArtTS.noiseTs = ...
                    ArtTS.NoiseGaussSD * randn(size(x))...
                    + ArtTS.NoiseGaussMean ...
                    + NoiseUniformRange * randn(size(x)) ...
                    + ArtTS.NoiseUniformLow;
                x = x + ArtTS.noiseTs;
            end
               
            ArtTS.Ts = x;                
        end
  

    % set methods
    function ArtTS = set.Name(ArtTS,str)
        if ischar(str)
            ArtTS.Name = str;
        else
            error('ArtTS.Name must be a string')
        end
    end
    
    function ArtTS = set.Description(ArtTS,str)
        if ischar(str)
            ArtTS.Description = str;
        else
            error('ArtTS.Description must be a string')
        end
    end
    
    function ArtTS = set.T0(ArtTS,dbl)
        if isscalar(dbl) && isreal(dbl)
            ArtTS.T0 = dbl;
        else
            error('ArtTS.T0 must be a real scalar value');
        end
    end
    
    function ArtTS = set.Extent(ArtTS,dbl)
        if isscalar(dbl) && isreal(dbl)
            ArtTS.Extent = dbl;
        else
            error('ArtTS.Extent must be a real scalar value');
        end
    end
    
    function ArtTS = set.nSamples(ArtTS,int)
        if isscalar(int) && isinteger(int) && int > 0
            ArtTS.nSamples = double(int);
        else
            error('ArtTS.nSamples must be a positive integer scalar value');
        end
    end
    
    function ArtTS = set.Freqs(ArtTS,dbl)
        if  isreal(dbl)
            ArtTS.Freqs = dbl;
        else
            error('ArtTS.Freqs must be an array of real values');
        end
    end
    
    function ArtTS = set.Amps(ArtTS,dbl)
        if  isreal(dbl)
            ArtTS.Amps = dbl;
        else
            error('ArtTS.Amps must be an array of real values');
        end
    end
    
    function ArtTS = set.Phases(ArtTS,dbl)
        if  isreal(dbl)
            ArtTS.Phases = dbl;
        else
            error('ArtTS.Phases must be an array of real values');
        end
    end
 
    end
    
%    methods  (Access=private)
%         
%         function time = makeTime(ArtTS)
%             time = (ArtTS.T0:ArtTS.Extent/ArtTS.nSamples:ArtTS.Extent);
%         end
%     end
        
    
end
                             