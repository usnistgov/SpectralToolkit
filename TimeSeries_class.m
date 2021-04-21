classdef TimeSeries_class < handle
    % Contains one time series channel, timeseries properties, and some 
    
    properties
        Ts      % the Time Series ( a matlab timeseries object)        
        dataTable   % original data channels read in from a file
    end
    
    methods
        %% Constructor
        
        function self = TimeSeries_class(varargin)
            %TimeSeries_class Construct an instance of this class
            
            % parses the incoming arguments
            validFileTypes = {'mat','csv','excel'};
            checkFileType = @(x) any(validatestring(x,validFileTypes));
            
            defaultTimeInfo = struct('Increment',1,'Start',0,'Units','seconds');
            
            p = inputParser;
            addParameter(p,'TimeInfo',defaultTimeInfo,@isstruct)
            addParameter(p,'File','',@ischar)
            
            parse(p,varargin{:})
            
            if ~strcmp(p.Results.File,'')
                % Load  dataTable from the file
                if exist(p.Results.File,'file')
                    [~,~,ext] = fileparts(p.Results.File);
                    if strcmp(ext,'.csv') || strcmp(ext,'.xls')
                        disp('Loading Data Table')
                        opts = detectImportOptions(p.Results.File);
                        opts.PreserveVariableNames = true;  % allows underscores
                        self.dataTable = readtable(p.Results.File,opts);
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
            
            
         %% Plot override methods
         function plot(self)
             plot(self.Ts)             
         end
            
         %% Loads a table from "File" and asks user which column to load
         function loadTs(self,timeVals,colName)
          % Query's the user for a channel name and creates the timeseries TS
          % Timeseries creation requires a time vector of the same length as the data
          
          % colName is optional, if not included, display the columns and query the user
          if nargin == 2
              headers = {self.dataTable.Properties.VariableNames{:}};
              disp('Data Header Names')
              disp(headers(:))
              colName = input('Time series to import: ');
          end
          
          dataVals = self.dataTable.(colName);
                  
          self.Ts = timeseries(dataVals,timeVals,'Name',colName);
          
          % NB:  TimeInfo.Increment is supposed to auto-update.  If the timeVals
          % are not uniformly spaced, the value will be "NAN".  
          % However I found that with a timeVals standard deviation of only 5 x 10^-16 (floating point errors)
          % the value would be "NAN".  So here we will check for ourselves
          if isnan(self.Ts.TimeInfo.Increment)
            if std(diff(timeVals))<1e-12
                %self.Ts.TimeInfo.Increment = mean(diff(timeVals));     % being depricated by mathworks
                self.Ts = setuniformtime(self.Ts,'Interval',mean(diff(timeVals)));
            end
          end
          
         end
         
         %% Creates a time vector from TimeInfo structure data
         function timeVals = timeInfo2TimeVals(self,timeInfo)
             timeVals = (timeInfo.Start:size(self.dataTable,1)+timeInfo.Start-1)*timeInfo.Increment;
         end
         
         %% Switch to another time series from the data table
         function switchTs(self)
             % switches the Ts to a new column of the data table.  
             % This uses the same time vector as in the existing TS
             timeVals = self.Ts.Time;
             self.loadTs(timeVals)
         end
         
 
    end
end

