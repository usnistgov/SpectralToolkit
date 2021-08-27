function window = getWindow(obj,offset,analysisCycles,varargin)
% return a subset of the timeseries data of length analysisCycles*Fs/F0
% offset is the number of nominal cycles to the start index, if the end of the data
% is reached, the data wraps and data from the beginning of the ts is
% fetched. Think of the timeseries data as a circular buffer, the offset
% can be multiples of the data that is in the ts.
    
    nWindow = analysisCycles*obj.SampleRate/obj.F0;
    
    if nargin > 3
        if strcmp(varargin{1},'odd')
            nWindow = nWindow+1;
        end
    end
    
    nData = size(obj.Ts.Data,1);
    startIndex = mod(offset*obj.SampleRate/obj.F0,nData)+1; 
    
    % If the time series contains enough data:
    if (nData-(startIndex-1)) >= nWindow
        window = obj.Ts.Data(startIndex:nWindow+(startIndex-1),:);
    else
        % not enough data in the Ts. wrap and concatinate
        W = obj.Ts.Data(startIndex:end,:);
        while length(W) < nWindow
            nRemaining = nWindow - length(W);
            if nRemaining >= nData
               W = vertcat(W,obj.Ts.Data);
            else
               W = vertcat(W,obj.Ts.Data(1:nRemaining,:));
            end
            window = W;
        end
        %plot(real(window))
    end


end