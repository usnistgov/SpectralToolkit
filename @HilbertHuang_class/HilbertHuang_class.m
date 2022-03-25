classdef HilbertHuang_class
    % For one time series, creates and contains a set of intrinsic
    % functions and Hilbert transform for each.  Optional will perform an EMD on a real timeseries of 
    % can receive a set of IMFs for analysis and visualization.
    %   Input <Name,Value> pairs:
    %       "Name" - char array name for this analysis
    %       "Timeseries" - timeseries object to be analysed
    %       "EMD" - true/false logical, if true will perform an EMD on the first column of data in the timeseries
    %       "IMFs" - timeseries object of IMFs
    %       "Hilbert" - true/false logical to perform EMD followed by IMF on a real timeseries
    %       "Window" - char array Window type to use on input data:  
    %       "EmdOpts" - EMD options structure (see the documentation of the emdOpts property).
    %
    %   Output (see the IMFs and Hilbert Properties
    %
    %
    %
    
    properties
        Ts_In           % Timeseries input
        IMFs            % Array of IMFs: (1 IMF per column ) 
                        % IMF{1,:) = signal
                        % IMF(1,:) = 
        Hilbert     % Cell Array of hilbert transform (complex) for each Intrinsic function 
                        % Hilbert{1,:} = Analytic Signal (complex)
                        % Hilbert{2,:} = time series of Energy of each IMF (abs(signal).^2
                        % Hilbert{3,:) = time series of phaseAngles for each IMF
                        % Hilbert{4,:} = time series of omega (angular frequency) values
                        %                each value falls between each time step in the time vector 
        Residual    % vector of input minus recombined hilbert transforms
        Window      % if set to a recognized window type, will window the IMF before Hilbert, defaults to 'none'
        
        
        % option structure for EMD:
            %'MaxNumIMFs' (default 10)
            %'MaxEnergyRatio' (default 20) Sifting stops whenthe ratio of the original signal energy to the IMF energy os greater than this
            %'MaxNumExtrema' (default 1) Sifting stops then the mumber of exterma is less than this
            %'SiftMaxIterations'
            %'SiftRelativeTolarance'
            %'Interpolation'            
        emdOpts = struct('MaxNumIMFs',10,...
                         'MaxEnergyRatio',20,...
                         'MaxNumExtrema',1,... 
                         'SiftMaxIterations',100,...
                         'SiftRelativeTolarance',0.2,...
                         'Interpolation','spline');
                     
        % option structure for Hilbert Transform
            % FreqLimits (default []), limits computing the Hilbert Spectrum
            % FreqResolution ((FreqLimits(2)-FreqLimits(1))/100),, discretizes the frequency limits
            % MinThreshold (defualt -inf), any hilber spectrum value less than this will be set to 0.
        hilOpts = struct('FreqLimits',[],...
                         'FreqResolution',[],...
                         'MinThreshold',-inf);
        
        fig = 1;    % a counter for the number of figures
    end
    
    
    %% --------------------------------------------------------------------
    % Constructor
    methods
        function obj = HilbertHuang_class(varargin)
            %HilbertHuang_class Construct an instance of this class
            % useage HilbertHuang_class(<name, value>)
            % optional arguments are name/value pairs:
            % name , (name)
            % TimeSeries, timeseries input to EMD. Must be a matlab timeseries type), Only the first column wil be analysed
            % EMD, (true or false): Perform an Empirical Mode Decomposition but not an HHT.
            % IMFs, timeseries containing intermediate mode functions, needed input if EMD is not performed
            % Hilbert, (true or false): Perform a Hilbert Transform on the timeseries (EMD followed by HHT)
            % 
            
            % checks for data validity
            checkTsType = @(x) isa(x,'timeseries');
                       
            % parse the incoming arguments
            p = inputParser;
            addParameter(p,'TimeSeries',timeseries(),checkTsType)
            addParameter(p,'EMD',false,@islogical)
            addParameter(p,'IMFs',timeseries(),checkTsType)
            addParameter(p,'Hilbert',false,@islogical)
            addParameter(p,'Window','none',@ischar)
            addParameter(p,'EmdOpts',obj.emdOpts,@isstruct)
            addParameter(p,'HilOpts',obj.hilOpts,@isstruct)
            
            parse(p,varargin{:})
            obj.Ts_In = p.Results.TimeSeries;
            obj.IMFs = p.Results.IMFs;
            obj.Window = p.Results.Window;
            obj.emdOpts = p.Results.EmdOpts;
            obj.hilOpts = p.Results.HilOpts;
            
            
            if p.Results.EMD
                obj = obj.Emd();
            end
            
            if p.Results.Hilbert
                obj = obj.Hht();
            end
                
            
            % is a timeseries was passed in with the arguments, then
            % calculate the HHT
%             if isempty(find(strcmp(p.UsingDefaults,'TimeSeries'),1))
%                 obj.Time = p.Results.TimeSeries.Time;
%                 obj.TimeInfo = p.Results.TimeSeries.TimeInfo;
%                 obj = obj.calcHHT(p.Results.TimeSeries);
%             end
        end
    end
 
     
%% -------------------------------------------------------------------------
%Public methods to calculate the IMFs using EMD and to calculate the Hilbert Transform
% also a public methog to create plots
methods (Access = public)
    obj = Emd(obj)
    obj = Hht(obj)
    plot(obj,varargin)
end

%% -------------------------------------------------------------------------
methods (Access = public)
   
    function obj = calcHHT(obj,timeSeries)
        % calculates the IMFs and HHTs from the first column of timeseries
        % data
        obj.Name = timeSeries.Name;
        obj.Time = timeSeries.Time;
        

        % EMD of the first column of timeseries data
        obj = obj.Emd(timeSeries.Data(:,1));
        obj = obj.Hht();        
    end
    
end

%% -------------------------------------------------------------------------
methods (Static)
    % Static methods to get default options values
    function [EmdOpts,HilOpts] = getDefaultOpts()
        EmdOpts = struct('MaxNumIMFs',10,...
            'MaxEnergyRatio',20,...
            'MaxNumExtrema',1,...
            'SiftMaxIterations',100,...
            'SiftRelativeTolarance',0.2,...
            'Interpolation','spline');
        
        HilOpts = struct('FreqLimits',[],...
            'FreqResolution',[],...
            'MinThreshold',-inf);
        
    end
    
end

end
