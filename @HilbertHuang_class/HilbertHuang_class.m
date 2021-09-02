classdef HilbertHuang_class
    % For one time series, creates and contains a set of intrinsic
    % functions
    % and Hilbert transform for each
    %   Input:
    %       Name
    %       Timeseries
    %
    %   Pro
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
            % EMD, (true or false): Perform and Empirical Mode Decomposition
            % IMFs, timeseries containing intermediate mode functions, needed input if EMD is not performed
            % Hilbert, (true or false): Perform a Hilbert Transform on the IMFs
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
            
            parse(p,varargin{:})
            obj.Ts_In = p.Results.TimeSeries;
            obj.IMFs = p.Results.IMFs;
            obj.Window = p.Results.Window;
            
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
 
     
%%-------------------------------------------------------------------------
%Public methods to calculate the IMFs using EMD and to calculate the Hilbert Transform
% also a public methog to create plots
methods (Access = public)
    obj = Emd(obj)
    obj = Hht(obj)
    plot(obj,varargin)
end

%%-------------------------------------------------------------------------
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
end

