classdef TimeSeries_class < handle
    %Containt one time series channel, timeseries properties, and some 
    % parameters estimated from the time series
    %   Detailed explanation goes here
    
    properties
        Ts      % the Time Series ( a matlab timeseries object)
        Fs      % FourierSeries
        
        dataTable   % original data channels read in from a file
        
                
    end
    
    methods
        %% Constructor
        
        function self = TimeSeries_class(varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            % parses the incoming arguments
            validFileTypes = {'mat','csv','excel'};
            checkFileType = @(x) any(validatestring(x,validFileTypes));
            
            defaultTimeInfo = struct('Increment',1,'Start',0,'Units','seconds');
            
            p = inputParser;
            addParameter(p,'TimeInfo',defaultTimeInfo,@isstruct)
            addParameter(p,'File','',@ischar)
            addParameter(p,'FileType','csv',checkFileType)
            
            parse(p,varargin{:})
            
            if ~strcmp(p.Results.File,'')
                % Load  dataTable from the file
                if exist(p.Results.File,'file')
                    [~,~,ext] = fileparts(p.Results.File);
                    if strcmp(ext,'.csv') || strcmp(ext,'.xls')
                        disp('Loading Data Table')
                        self.dataTable = readtable(p.Results.File);
                        disp ('Data Table Load Complete')
                    else
                        warning('file %s does not exist',file)
                    end
                end
            end
            
            % if the user passed in a timeInfo struct AND the dataTable is not 
            % empty, create a time vector and load the Ts time series
            if isempty(find(strcmp(p.UsingDefaults, 'TimeInfo'),1))
                if ~isempty(self.dataTable)
                    timeVals = self.timeInfo2TimeVals(p.Results.TimeInfo);
                    loadTs(self,timeVals)                    
                end
                
            end
            
            
            
        end
            
            
         %%  
            
         %% Loads a table from "File" and asks user which column to load
         function loadTs(self,timeVals)
          % Query's the user for a channel name and creates the timeseries TS
          % Timeseries creation requires a time vector of the same length as the data
          headers = {self.dataTable.Properties.VariableNames{:}};          
          disp('Data Header Names')
          disp(headers(:))
          colName = input('Time series to import: ');
          dataVals = self.dataTable.(colName);
                  
          self.Ts = timeseries(dataVals,timeVals,'Name',colName);
         end
         
         %% Creates a time vector from TimeInfo structure data
         function timeVals = timeInfo2TimeVals(self,timeInfo)
             timeVals = (timeInfo.Start:size(self.dataTable,1)+timeInfo.Start-1)*timeInfo.Increment;
         end
         
 
    end
end

