classdef ModFit_class < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TS      % Time series structure with the full duratio of signal
        Window  % Subset of the TS data to be analysed
        SignalParams    % Parameters of the signal in the TS
        Duration        % Duration of the TS in seconds
        Fs              % SampleRate of the TS
        F0              % Nominal Frequency (50 Hz or 60 Hz)
        AnalysisCycles  % number of nominal cycles in the Window to be analysed 
        DlyCorr         % Delay Correction Factor(s)
        MagCorr         % Magnitude Correction Factors
        fig = 1;        % a counter for figure numbers
    end
    
    
    %% --------------------------------------------------------------------
    % Constructor    
    methods
        function obj = ModFit_class(varargin)
            defaultName = 'unnamed';
            defaultF0 = 50;
            defaultFs = 4800;  % samples per second
            defaultDuration = 1;       % duration of the data (in Ts.TimeInfo.Units)
            defaultAnalysisCycles = 6;
        
            defaultSignalParams = zeros(15,1);   % default to only 1 phase of signal
            defaultSignalParams(1,:) = 1;       
            defaultSignalParams(2,:) = 50;
            
            defaultMagCorr = ones(1,6);
            defaultDlyCorr = zeros(1,6);
            p = inputParser;
        
            % validation
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %validScalarPosInt = @(x) isnumeric(x) && isscalar(x) && isinteger(x) && (x > 0);
            validateSignalParams = @(x) validateattributes(x,{'double'},{'nrows',15});
            validateCorrVals = @(x) validateatributes(x,{'double'},{'nrows',6});
            
            addParameter(p,'Name',defaultName,@ischar);
            addParameter(p,'F0',defaultF0,validScalar);
            addParameter(p,'Fs',defaultFs,validScalarPosNum);
            addParameter(p,'AnalysisCycles',defaultAnalysisCycles,validScalarPosNum);            
            addParameter(p,'Duration',defaultDuration,validScalarPosNum);
            addParameter(p,'SignalParams',defaultSignalParams,validateSignalParams);
            addParameter(p,'DlyCorr',defaultDlyCorr,validateCorrVals);
            addParameter(p,'MagCorr',defaultMagCorr,validateCorrVals);
                       
            parse(p,varargin{:})
            
            obj.F0 = p.Results.F0;
            obj.Fs = p.Results.Fs;
            obj.Duration = p.Results.Duration;
            obj.SignalParams = p.Results.SignalParams;
            obj.AnalysisCycles = p.Results.AnalysisCycles;
            obj.DlyCorr = p.Results.DlyCorr;
            obj.MagCorr = p.Results.MagCorr;
            
            %Instantiate a default Time series
            obj.TS = AnalyticTS_class(...
                                        'Name',p.Results.Name,...
                                        'F0',obj.F0,...
                                        'SampleRate',obj.Fs,...
                                        'Duration',obj.Duration,...
                                        'SignalParams',obj.SignalParams);
                   
        end
        
    end
%% -----------------------------------------------------------------------
    % Public Methods from external .m files
    methods (Access = public)
        [Synx,Freq,ROCOF,iterations] = ModulationFit(obj)
        
    end
    
    % Static Methods from external .m files
    methods(Static)
         [Synx,Freqs,ROCOFs, iterations] = ModFit2Sb(obj,Fin,Fm,Km,Samples,dt,MagCorr,DelayCorr);
         [Synx,Freqs,ROCOFs, iterations] = ModFitNSb(Fin,Fm,Km,Samples,dt,MagCorr,DelayCorr)
    end
    
%% -----------------------------------------------------------------------
    % Public Methods
    methods (Access = public)
       
        %==================================================================
        function obj = getTimeSeries(obj)
            % TS is a timeseries that creates and contains one period of an
            % analytic arbitrary waveform. That period can be used as a
            % circular buffer.  It has a getWindow method which will return
            % AnalysisCycle nominal periods of that waveform with an offset
            % in nominal periods.
            obj.TS = AnalyticTS_class(...
                                        'SignalParams',obj.SignalParams,...
                                        'SampleRate',obj.Fs,...
                                        'F0',obj.F0,...
                                        'Duration',obj.Duration...
                                        );                                                                                            
        end
        
        %==================================================================
 
         function runMod1Second(obj,bDisplay)
            % run the modularon fitter for 1 second and verify the TVE, 
            % FE, and RFE with the expected result.  Optionally display the 
            % results and pause
            disp(obj.TS.Ts.Name)
                                    
            % pre allocate actual (returned fitter) data arrays and expected value arrays
            actSynx = zeros(length(obj.SignalParams(1,:))+2,obj.F0);
            expSynx = actSynx;            
            actFreq = zeros(1,obj.F0);
            expFreq = actFreq;
            actROCOF = actFreq;
            expROCOF = actROCOF;
            
            % run a loop for 1 second of simulated data, comparing the fitted SynX against the center values of the window
            for i = 1:obj.F0
                %obj.getWindow(i*(obj.Fs/obj.F0));
                
                obj.Window = obj.TS.getWindow(i,obj.AnalysisCycles,'even');
                [actSynx(:,i), actFreq(i), actROCOF(i), iter] = obj.ModulationFit();
                                               
                %expSynx(:,i) = obj.calcSymComp(obj.Window.UserData.Vals.')/sqrt(2);
                expSynx(:,i) = obj.Window.UserData.Vals/sqrt(2);
                expFreq(i) = obj.Window.UserData.Freqs;
                expROCOF(i) = obj.Window.UserData.ROCOFs;
                
                msg = sprintf('Cycle: %d, Iterations: %d',i,iter);
                disp(msg);
            end
            
            
            act = struct('Synx',actSynx,'Freq',actFreq,'ROCOF',actROCOF);
            exp = struct('Synx',expSynx,'Freq',expFreq,'ROCOF',expROCOF);
            if bDisplay == true, obj.dispErrors(act,exp,obj.fig),end
            
        end
       
        %==================================================================
        
        
    end
    
%% ------------------------------------------------------------------------    
% Static methods
    methods(Static)
        %=======================================
        function dispErrors(act,exp,fig)
            TVE = sqrt(((real(act.Synx)-real(exp.Synx)).^2+(imag(act.Synx)-imag(exp.Synx)).^2)./(real(exp.Synx).^2+imag(exp.Synx).^2))*100;
            ME =  (abs(act.Synx)-abs(exp.Synx))./ abs(exp.Synx)*100;
            PE = wrapToPi(angle(act.Synx)-angle(exp.Synx)).*(180/pi);
            %PE = wrapToPi(angle(act.Synx)-angle(exp.Synx))./(2*pi.*exp.Freq); 
            FE = act.Freq-exp.Freq;
            RFE = act.ROCOF-exp.ROCOF;
            
            figure(fig);
            subplot(5,1,1)
            plot(TVE');
            title('TVE (%)')
            subplot(5,1,2)
            plot(ME');
            title('ME (%)')
            subplot(5,1,3)
            plot(PE');
            title('PE (sec)')
            subplot(5,1,4)
            plot(FE');
            title('FE (Hz)')
            subplot(5,1,5)
            plot(RFE');
            title('RFE (Hz/s)')
                        
        end
    
        %=======================================
        function dispFourier(ts,fig)
            FT = FourierSeries_class('Timeseries',ts,'Name',ts.Name);
            figure(fig)
            subplot(2,1,1)
            FT.plot('SingleSided')
            xlim([0,100])
            subplot(2,1,2)
            FT.plot('Phase')
            xlim([0,100])
            
            % Plot a 3D polar view of the FFT
            figure(fig+1)
            FT.plot('Polar3')
            xlim([0,100])
            ylim([-1,1])
            zlim([-1,1])
            view([-57,27.8])
                        
        end

        
    end

end

