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
        Name
        IMFs         % Cell Array of IMFs: (1 IMF per column ) 
                        % IMF{1,:) = signal
                        % IMF(1,:) = 
        Hilbert     % Cell Array of hilbert transform (complex) for each Intrinsic function 
                        % Hilbert{1,:} = Analytic Signal (complex)
                        % Hilbert{2,:} = time series of Energy of each IMF (abs(signal).^2
                        % Hilbert{3,:) = time series of phaseAngles for each IMF
                        % Hilbert{4,:} = time series of omega (angular frequency) values
                        %                each value falls between each time step in the time vector 
        Time        
        Residual    % vector of input minus recombined hilbert transforms
    end
    
    %% --------------------------------------------------------------------
    % Constructor
    methods
        function obj = HilbertHuang_class(varargin)
            %HilbertHuang_class Construct an instance of this class
            % useage HilbertHuang_class(<name, value>)
            % optional arguments are name/value pairs:
            % name , (name)
            % TimeSeries, (timeseries analysed (must be a matlab timeseries type)
            defaultName = 'unamed';
            
            % checks for data validity
            checkTsType @(x) isa(x,'timeseries')
                       
            % parse the incoming arguments
            p = inputparser;
            addParameter(p,'Name',defaultName,@ischar)
            addParameter(p,TimeSeries',timeseries(),checkTsType)
            
            parse(p,varargin{:})
            obj.Name=p.Results.Name;
            
            % is a timeseries was passed in with the arguments, then
            % calculate the HHT
            if isempty(find(strcmp(p.UsingDefaults,'TimeSeies'),1))
                obj.calcHHT(p.Results.TimeSeries)
            end
        end
    end
 
     
%%-------------------------------------------------------------------------
%Public methods to calculate the IMFs using EMD and to calculate the Hilbert Transform
methods (Access = public)
    obj = Emd(obj,x)
    obj = Hht(obj)
end

%%-------------------------------------------------------------------------
methods (Access = public)
   
    function calcHHT(obj,timeSeries)
        % calculates the IMFs and HHTs from the first column of timeseries
        % data
        obj.Name = timeSeries.Name;
        obj.Time = timeSeries.Time;
        

        % EMD of the first column of timeseries data
        obj.Emd(timeSeries.Data(:,1));
        obj.Hht();        
    end
    
end


