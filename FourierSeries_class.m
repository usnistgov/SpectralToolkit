classdef FourierSeries_class < handle
    % Contains one Fourier Series, 
    % note, the series is double-sided:  DC is at index 1, nyquist is at N/2 and 
    % indices N/2+1 to end is a mirror image of data from 1 to N/2.  
    
    properties
        Name             % Fourier Series Name
        Data             % Fourier Series Data vector (complex)
        Freq             % Fourier Series Frequency vector (must be same length as Data)
        DataInfo         % structure of information about the Fourier Data vector
        FreqInfo         % structure of information about the Fourier Frequency vector
        
        Residual         % a timeseries of the residual after inverse transform
    end
    
    methods
        
        %% Constructor
        function self = FourierSeries_class(varargin)
            %FourierSeries_class Construct an instance of this class
            
            % default property values
            defaultName = 'unnamed';
            defaultDataInfo = struct(...
                'Length',-1,...         % Length of the Data, -1 means that it has not yet been calculated
                'TsIndex',[1,-1],...   % the start and end indices of the transformed time series, end = -1 means the entire time series will be used                
                'Window','Rect'...     % The window type used, defaults to rectangular (no window)
                );
            defaultWindow = defaultDataInfo.Window;
            defaultTsIndex = defaultDataInfo.TsIndex;
            
            defaultFreqInfo = struct(...
                'Length',-1, ...
                'BinWidth',-1 ...
            );
            
            % check data validity
            checkTsType = @(x) isa(x,'timeseries');
                
            
            % parses the incoming arguments
            p = inputParser;
            addParameter(p,'DataInfo',defaultDataInfo,@isstruct)
            addParameter(p,'FreqInfo',defaultFreqInfo,@isstruct)
            addParameter(p,'TimeSeries',timeseries(),checkTsType)
            addParameter(p,'Name',defaultName,@ischar)
            addParameter(p,'Window',defaultWindow,@ischar)
            addParameter(p,'TsIndex',defaultTsIndex,@isnumeric)
            
            parse(p,varargin{:})
            
            self.Name = p.Results.Name;
            self.DataInfo = p.Results.DataInfo;
            self.DataInfo.Window = p.Results.Window;
            self.DataInfo.TsIndex = p.Results.TsIndex;
            self.FreqInfo = p.Results.DataInfo;
            
            % if a timeseries was passed in, then calculate the Fourier Series
            if isempty(find(strcmp(p.UsingDefaults, 'TimeSeries'),1))
                self.calcFs(p.Results.TimeSeries)                
            end
            
        end
        
        %%-----------------------------------------------------------------
        function plot(self,varargin)
            % plot override functions
            yscale = 'linear';
            yMsg = 'Amplitude';
            r = [];
            N = self.DataInfo.Length;
            N_2 = floor(self.DataInfo.Length/2);
            
            i = 1;
            while i <= numel(varargin)
                switch varargin{i}
                    case 'DoubleSided'
                        x = self.Freq;
                        y = abs(self.Data/N);
                    case 'SingleSided'
                        x = self.Freq(1:N_2+1);
                        y = abs(self.Data/N);
                        y = y(1:N_2+1);
                        y(2:end-1) = 2*y(2:end-1);
                    case 'Shift'
                        x = [-flip(self.Freq(2:N_2+1)); self.Freq(1:N_2)];
                        y = fftshift(abs(self.Data/N));
                    case 'xlim'
                        i = i+1;                        
                        xlims = varargin{i};  
                    case 'ylim'
                        i = i+1;
                        ylims = varargin{i};
                    case {'residual','Residual'}
                        r = self.Residual;
                    case {'yscale'}
                        i = i+1;
                        switch varargin{i}
                            case 'log'
                                yMsg = ('Amplitude (log_10)');
                                yScale = 'log';
                            case 'dB'
                                yMsg = ('Amplitude (dB)');
                                yScale = 'dB';
                        end
                        
                     otherwise
                        warning('Plot argument %s not supported',varargin{i})
                end
                i = i+1;
            end
            
            
            if ~isempty(r)      % residual?
                subplot(2,1,2)
                plot(r)
                title('Residual')
                subplot(2,1,1)
            end         
            
            if strcmp(yScale,'dB')
               y = 20*log10(y); 
            end
            
            % plot function            
            h=stem(x,y);       % you should never plot Frequency series with interpolated lines
            title(strcat('Frequency Series: ',self.Name),'Interpreter','none')
            msg = sprintf('Frequency (Hz) (Bin Width %5.2f)',self.FreqInfo.BinWidth);
            xlabel(msg)               
            ylabel(yMsg)
                                    
            if strcmp(yScale,'log')
                set(gca,'yscal','log')
            end
            
            if strcmp(yScale,'dB')
                set(h,'BaseValue',-120);  % assumes Digitizer dynamic range is 120dB
            end

            
            if exist('xlims','var')
                set(gca,'xlim',xlims)
            end
            if exist('ylims','var')
               set(gca,'ylim',ylims) 
            end
        end
        
            
        
        %% ----------------------------------------------------------------
        % use with caution, these windows have not all been verified yet.
        function calcFs(self,timeSeries)
            self.Name = timeSeries.Name;
            % calculates frequency series (Fourier) data from the time series
           
            % first grab the timeSeries data 
            if self.DataInfo.TsIndex(2) == -1
                self.DataInfo.TsIndex(2) = timeSeries.Length;
            end
            tsData = timeSeries.Data(self.DataInfo.TsIndex(1):self.DataInfo.TsIndex(2));
            dt = timeSeries.TimeInfo.Increment;
            N = length(tsData);
            self.Freq = ((0:N-1)/(N*dt))';
            self.FreqInfo.Length = length(self.Freq);
            self.FreqInfo.BinWidth = 1/(N*dt);
            
            % Create a Window function
            switch self.DataInfo.Window
                case 'Rect'
                    w = ones(length(tsData),1);
                case 'Hamming'
                    w = self.Hamming(length(tsData)); 
                case 'Hann'
                    w = self.Hann(length(tsData));
                case 'Blackman'
                    w = self.Blackman(length(tsData));
                otherwise
                    warning('Window type: %s is not supported, Rect will be used',self.DataInfo.Window)
                    w = ones(length(tsData),1);
            end
            
            tsData = tsData.*w;     % window the data            
            self.Data = fft(tsData);    
            self.DataInfo.Length = length(self.Data);
            
            % calculate the inverse fourier transform and the residual
            x = ifft(self.Data);
            self.Residual = tsData - x;
            
        end
        
        
%% ============================================================================
%  % WORK IN PROGRESS, NOT READY TO USE YET
%        % Calculate (and plot?) the expectation value due to the windowing function
%         function Ws = LeakFunc(self, Window)
%             % to start withm only the rectangular function
%             M = 1000;
%             switch Window
%                 case 'Rect'
%                     w = ones(M,1);
%                 case 'Hamming'
%                     w = self.Hamming(M);
%                 case 'Hann'
%                     w = self.Hann(M);
%                 case 'Blackman'
%                     w = self.Blackman(M);
%                 otherwise
%                     warning('Window type: %s is not supported, Rect will be used',self.DataInfo.Window)
%                     w = ones(M,1);
%             end
%             plot(w)
%             
%             Wss = M*sum(w.^2);
%             s = (0:M-1)/100; %
%            W = (1/M^2)*((sin(pi*s)./sin((pi*s)/M)).^2); %Rectangular
% %             for k = 1:M
% %                 W(k) = (1/Wss)*(sum(w(k).*exp(1i*2*pi*s*k/M)))^2;
% %             end
% %                           
%             semilogy(s,W);
%             set(gca,'ylim',[1e-6; 1])
%         end
%% ============================================================================
         
    
           
        
        
 
    end
    
    methods (Static)
        
        % Window functions override the Signal Processing Toolbox in case it is not present
        function w = Hamming(L) 
            N = L-1;
            w = zeros(N,1);
            for n = -N/2:N/2
                w(n+1+N/2,1) = 0.54+0.46*cos(2*pi*n/N);
            end
        end
        
        function w = Hann(L)
            N = L-1;            
            w = zeros(L,1);
            for n=0:N
               w(n+1) = 0.5*(1-cos(2*pi*n/N));
            end
        end
        
        function w = Blackman(L)
            %N = round(L/2);
            N = L;
            w = zeros(N,1);
            for n = 0:N-1
                w(n+1) = 0.42-0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1));
            end
%             for n = 1:N-1
%                 w(N+n) = w(N-n);
%             end 
%             w(L)=w(1);
        end
    end
  
    
end

