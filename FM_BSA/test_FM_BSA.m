classdef test_FM_BSA < matlab.unittest.TestCase
    % Class-based unit testing and experiments for FM_BSA_Class
    %
    %   run these tests with two command-line commands:
    %   >> testCase = test_FM_BSA();
    %   >> res = run(testCase); OR
    %   >> res = testCase.run;
    
     properties
        Name    % Test name
        T0             % start time
        F0             % nominal frequency
        Fs             % sample rate
        AnalysisCycles % number of fundamental cycles to be analysed
        SignalParams    % Parameters input to PmuTestSignals
        Duration    % maximum signal duration is seconds
        SettlingTime
        TS      % Time series structure to be analysed
        DlyCorr % Delay correction factors
        MagCorr % Magnitude correction factors
        expect     %expected values
        fig = 1 % figure numbers
        Window  % the window of data to be analysed
        even = true;
        verbose
        debug
        makeAnimation = false;
     end  
    
    % % Signal params.  Note that the labeling convention comes mostly from the
    % PMU standard
    % At the bottom of this class is a static methods to retrieve the
    % indexes to the parameters.
    %     Xm = signalparams(1,:)*sqrt(2);     % phase amplitude (given by the user in RMS
    %     Fin = signalparams(2,:);    % frequency (must be the same for all 6 channels or an error will be thrown
    %     Ps = signalparams(3,:);     % phase
    %     Fh = signalparams(4,:);     % Frequency of the interfering signal
    %     Ph = signalparams(5,:);     % Phase of the interfering signal
    %     Kh = signalparams(6,:);     % index of the interfering signal
    %     Fa = signalparams(7,:);     % phase (angle) moduation frequency
    %     Ka = signalparams(8,:);     % phase (angle) moduation index
    %     Fx = signalparams(9,:);     % amplitude moduation frequency
    %     Kx = signalparams(10,:);    % amplitude moduation index
    %     Rf = signalparams(11,:);    % ROCOF
    %     KaS = signalparams(12,:);   % phase (angle) step index
    %     KxS = signalparams(13,:);   % magnitude step index
    %     KfS = signalparams(14,:);   % frequency step index
    %     KrS = signalparams(15,:);   % ROCOF step index (redundant with Rf)
     
    %% ----------------------------------------------------------------------
    %Constructor Method
    methods (Access = public)
        
        function self = test_FM_BSA(varargin)
        % Class Constructor
        
            defaultName = 'unnamed';
            defaultT0 = 0;             % start time in seconds
            defaultF0 = 50;
            defaultFs = 4800;  % samples per second
            defaultDuration = 1;       % duration of the data (in Ts.TimeInfo.Units)
            defaultAnalysisCycles =6;
            defaultSettlingTime = 0;   % time of static sinusoid before any dynamics begin, also added to the end of the dynamics.
            
            [Xm, Fin, Ps] = self.getParamIndex();
            defaultSignalParams = zeros(15,6);
            defaultSignalParams(Xm,:) = 1;
            defaultSignalParams(Fin,:) = 50;
            
            self.SignalParams(Ps,:) = [0,-120,120,0,-120,120];            
            defaultMagCorr = ones(1,6);
            defaultDlyCorr = zeros(1,6);
            p = inputParser;
            
            % validation
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %validScalarPosInt = @(x) isnumeric(x) && isscalar(x) && isinteger(x) && (x > 0);
            validateSignalParams = @(x) validateattributes(x,{'double'},{'nrows',15});
            validateCorrVals = @(x) validateatributes(x,{'double'},{'nrows',6});
            validateLogical = @(x) islogical(x);
            
            addParameter(p,'Name',defaultName,@ischar);
            addParameter(p,'T0',defaultT0,validScalar);
            addParameter(p,'F0',defaultF0,validScalar);
            addParameter(p,'Fs',defaultFs,validScalarPosNum);
            addParameter(p,'AnalysisCycles',defaultAnalysisCycles,validScalarPosNum);            
            addParameter(p,'Duration',defaultDuration,validScalarPosNum);
            addParameter(p,'SettlingTime',defaultSettlingTime,validScalarPosNum);
            addParameter(p,'SignalParams',defaultSignalParams,validateSignalParams);
            addParameter(p,'DlyCorr',defaultDlyCorr,validateCorrVals);
            addParameter(p,'MagCorr',defaultMagCorr,validateCorrVals);
            addParameter(p,'verbose',false,validateLogical);
            addParameter(p,'debug',false,validateLogical);
                        
            parse(p,varargin{:})
            
            self.Name = p.Results.Name;
            self.T0 = p.Results.T0;
            self.F0 = p.Results.F0;
            self.Fs = p.Results.Fs;
            self.Duration = p.Results.Duration;
            self.SettlingTime = p.Results.SettlingTime;
            self.SignalParams = p.Results.SignalParams;
            self.AnalysisCycles = p.Results.AnalysisCycles;
            self.DlyCorr = p.Results.DlyCorr;
            self.MagCorr = p.Results.MagCorr;
            self.verbose = p.Results.verbose;
            self.debug = p.Results.debug;
        end
        
    end
    
   %----------------------------------------------------------------------
    %% Test Methods
    % These functions will be called on   >> "res = run(testCase);"
    methods (Test)
        function regressionTests(self)
            self.fig = 1;
            test50f0_2m0_2a5(self); self.fig=self.fig+1; % Phase Modulation, fm = 2, k = 2.5
            test50f0_5m0_5a0(self); self.fig=self.fig+1; % Phase Modulation, fm = 2, k = 2.5
        end
        
        function experiments(self)
            FcarrDfContour(self) 
            %self.debug=true;  % If you want to see all the contour plots
            GridSearchThreshold(self)
        end
    end
    
    
    %% --------------------------------------------------------------------
    % The test methods themselves
    methods (Access = public)
        
        % Regression Tests
        %------------------------------------------------------------------
        function test50f0_2m0_2a5(self)
            % this is a very high rate FM modulation (peak frequency 5 Hz, peak ROCOF 62 Hz)
            self.setTsDefaults();
            self.Duration = 2;
            [ ~, ~, ~, ~, ~, ~, Fa, Ka, ~, ~] = self.getParamIndex();
            self.SignalParams(Fa,:) = 2.0;
            self.SignalParams(Ka,:) = 2.5;
            self.AnalysisCycles = self.F0/self.SignalParams(Fa,1);  % one modulation cycle
            self.getTimeSeries();
            self.TS.Ts.Name = 'test50f0_2m0_2k5';
            self.runMod1Second(true);           
        end    
        
        function test50f0_5m0_5a0(self)
            % this is a very high rate FM modulation (peak frequency 5 Hz, peak ROCOF 62 Hz)
            self.setTsDefaults();
            self.Duration = 2;
            [ ~, ~, ~, ~, ~, ~, Fa, Ka, ~, ~] = self.getParamIndex();
            self.SignalParams(Fa,:) = 5.0;
            self.SignalParams(Ka,:) = 5.0;
            self.AnalysisCycles = self.F0/self.SignalParams(Fa,1);  % one modulation cycle
            self.getTimeSeries();
            self.TS.Ts.Name = 'test50f0_5m0_5k0';
            self.runMod1Second(true);     
        end
        
        % Experiments
        %------------------------------------------------------------------
        function FcarrDfContour(self)
            % Create a 3D contour plot of the objective function over
            % carrier frequency vs peak frequency deviation while holding
            % the modulation phase at 0
            
            % First, set up a range of carrier frequencies and
            % delta frequencies for which to plot the objective function
            % contour

            self.Name = 'Fcarr dF Contour';
            self.SignalParams = zeros(15,1);        % only one phase of data
            [Xm, Fin, Ps, ~, ~, ~, Fa, Ka] = self.getParamIndex();
            self.SignalParams(Xm,:) = 1;
            self.SignalParams(Fin,:) = 50;
            self.SignalParams(Ps,:) = 0;
            
            self.SignalParams(Fa,:) = 5.0;
            self.SignalParams(Ka,:) = 5.0;
            
            dT = 1/self.Fs;    % sampling period
            dF = 2*pi*self.SignalParams(Fa,:)*self.SignalParams(Ka,:)*dT;
            self.AnalysisCycles =  ceil(self.SignalParams(Fin,1)/self.SignalParams(Fa,1));         
            
            self.getTimeSeries();
            self.Window = self.TS.getWindow(7,self.AnalysisCycles,self.even);
            
            
            FM = FM_BSA_Class(...
                self.SignalParams(Fin,:),...
                self.SignalParams(Fa,:),...
                self.SignalParams(Ka,:),...
                dT,...
                real(self.Window.Data)...
                );
            
            wCarr = 2*pi*self.F0*dT;
            
            figure(self.fig),self.fig = self.fig+1;
            omega_F_dF = [[wCarr/2,wCarr*2];[0,0];[0,2*dF]];
            FM.phase = 1;
            FM.fcontour3(omega_F_dF,30,@FM.objFun);
            
            figure(self.fig),self.fig = self.fig+1;
            omega_F_phi = [[wCarr/2,wCarr*2];[-pi,pi];[dF,dF];];
            FM.fcontour3(omega_F_phi,30,@FM.objFun);
            
            figure(self.fig),self.fig = self.fig+1;
            omega_phi_dF = [[wCarr,wCarr];[-pi,pi];[0,2*dF]];
            FM.fcontour3(omega_phi_dF,30,@FM.objFun);                       
        end
        
        function GridSearchThreshold(self)
            % For a range of modulation frequencies and indices, plot the
            % minimum value of the objective function.  Next use these
            % minimum calues to create a 3D plot versus the Fm and Km.
            %
            % the values of this plot, will then be used to find a function
            % which can predict good grid search stopping threshold values.

            self.Name = 'Grid Search Threshold Experiment';
            disp(self.Name)
            
            %self.Fs = 2400;
            self.Duration = 2;

            self.SignalParams = zeros(15,1);        % only one phase of data
            [Xm, Fin, Ps, ~, ~, ~, Fa, Ka] = self.getParamIndex();
            self.SignalParams(Xm,:) = 1;
            self.SignalParams(Fin,:) = 50;
            self.SignalParams(Ps,:) = 0;

            
            FmStart = 0.5;
            FmEnd = 5;
            FmIncr = 0.5;
            
            KmStart = 0.5;
            KmEnd = 5;
            KmIncr = 0.5;
            
            grid = 20;
            res = 30;
            
            dT = 1/self.Fs;
            wCarr = 2*pi*self.F0*dT;
            
            % Set up a table for the results
            nRows = (round((FmEnd-FmStart)/FmIncr)+1)*(round((KmEnd-KmStart)/KmIncr)+1);
            varNames = {'Km','Fm','zMin'};
            varTypes = {'double','double','double'};
            T = table('Size',[nRows,length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
            
            %---------------Movie Making----------------
            if self.makeAnimation
                vidfile = VideoWriter('GridSearchContour.mp4','MPEG-4');
                open(vidfile);
            end
            %-------------------------------------------            
            
            p = 1; % counter for the table
            if self.debug,figure(self.fig),self.fig=self.fig;end
            for Km = KmStart:KmIncr:KmEnd
                for Fm = FmStart:FmIncr:FmEnd
                    self.SignalParams(Fa,:) =Fm;
                    self.SignalParams(Ka,:) = Km;
                    dF = 2*pi*Km*Fm*dT;
                    self.AnalysisCycles =  ceil(self.SignalParams(Fin,1)/self.SignalParams(Fa,1));
                    
                    self.getTimeSeries();
                    self.Window = self.TS.getWindow(0,self.AnalysisCycles,self.even);
                    
                    FM = FM_BSA_Class(...
                        self.SignalParams(Fin,:),...
                        self.SignalParams(Fa,:),...
                        self.SignalParams(Ka,:),...
                        dT,...
                        real(self.Window.Data)...
                        );
                    
                    FM.phase = 1;
                    startPt1 = 2*pi*self.SignalParams(Fin,1)*dT;
                    OMEGA2 = linspace(-pi,pi,grid);
                    OMEGA3 = linspace(0,2*dF,grid);
                    
                    if self.debug
                        OMEGA = [startPt1,startPt1;-pi,pi;0,2*dF];
                        FM.fcontour3(OMEGA,res,@FM.objFun)
                        view([70,30])
                        zlim([-3e4,0])
                        ylim([0,0.04])
                        xlim([0,pi])
                        title(sprintf('Fm = %1.1f,     Km = %1.1f',Fm,Km),'FontSize',18)
                        hold on
                    end
                    
                    z = zeros(grid,grid);
                    for k = 1:grid
                        for j = 1:grid
                            z(j,k) = FM.objFun([startPt1,OMEGA2(j),OMEGA3(k)]);
                            if self.debug
                               plot3(OMEGA2(j),OMEGA3(k),z(j,k),'*') 
                            end
                        end                        
                    end
                    
                    if self.debug
                        hold off
                        refresh
                        %----------------- Write Video --------------------
                        if self.makeAnimation
                            frame = getframe(gcf);
                            for v = 1:8
                                writeVideo(vidfile,frame)
                            end
                        end
                        %--------------------------------------------------
                    end
                    
                    % Store the parameters and the minimum function value
                    T(p,:) = {Km,Fm,min(min(z))};
                    p = p+1;                                        
                    
                end
            end
            
            % Now perform a curve fit
            x = T{:,'Km'};
            y = T{:,'Fm'};
            zMin = T{:,'zMin'};
            
            % plot the surface
            %nX = length(unique(x)); nY=length(unique(y));
            [X,Y] = meshgrid(KmStart:KmIncr:KmEnd,FmStart:FmIncr:FmEnd);            
            zMat=reshape (zMin,length(FmStart:FmIncr:FmEnd),length(KmStart:KmIncr:KmEnd));
            figure(self.fig);self.fig=self.fig+1;
            surf(X,Y,zMat)
            colormap('parula')
            set(gca,'ZScale','log')
            xlabel('Km')
            ylabel('Fm')
            zlabel('Minimum objective function value')
            title('Minimum objective value contour')
            %---------- Movie making --------------------------------
            if (self.debug && self.makeAnimation)
                frame = getframe(gcf);
                for v = 1:24
                    writeVideo(vidfile,frame)
                end
            end
            %--------------------------------------------------------
            
            % Now we are going to fit the log function values using a CFIT object
            % ** WARNING**  This next bit of code requires the Fitting Toolbox.
            %fitFun = ('p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.*2');  %Second order polynomial surface
            SF = fit([x,y],real(-log(zMin)),'poly22');  % fit is our cfit object
            figure(self.fig);self.fig=self.fig+1;
            plot(SF,[x,y],real(-log(zMin)))
            colormap(flipud(parula));
            xlabel('Km');
            ylabel('Fm');
            zlabel('-log( objective function value)')
                      
         end
    end
    
    %% --------------------------------------------------------------------    
    % private methods
    methods (Access = private)
        
        %===================================
        function self=setTsDefaults(self)
            % Clear the SignalParams and set default values
            self.Name = 'Default';
            self.SignalParams = zeros(15,3);
            [Xm, Fin, Ps] = self.getParamIndex();
            self.SignalParams(Xm,:) = 1;
            self.SignalParams(Fin,:) = 50;
            self.SignalParams(Ps,:) = [0,-120,120];
        end
        
        
        %======================================
        function self = getTimeSeries(self)
            % TS is a timeseries that creates and contains one period of an
            % analytic arbitrary waveform. That period can be used as a
            % circular buffer.  It has a getWindow method which will return
            % AnalysisCycle nominal periods of that waveform with an offset
            % in nominal periods.
            self.TS = AnalyticTS_class(...
                                        'SignalParams',self.SignalParams,...
                                        'SampleRate',self.Fs,...
                                        'F0',self.F0,...
                                        'T0',self.T0,...
                                        'Duration',self.Duration,...
                                        'SettlingTime',self.SettlingTime...
                                        );                                                                                
        end
        
        %=======================================

        function runMod1Second(self,bDisplay)
            % run the modularon fitter for 1 second and verify the TVE, 
            % FE, and RFE with the expected result.  Optionally display the 
            % results and pause
            disp(self.TS.Ts.Name)
            [~,Fin,~,~,~,~,Fa,Ka] = self.getParamIndex();
                                    
            % pre allocate actual (returned fitter) data arrays and expected value arrays
            nSymComp = floor(size(self.SignalParams,2)/3);
            actSynx = zeros(length(self.SignalParams(1,:))+nSymComp,self.F0);
            expSynx = actSynx;            
            actFreq = zeros(1,self.F0);
            expFreq = actFreq;
            actROCOF = actFreq;
            expROCOF = actROCOF;
            
            % run a loop for 1 second of simulated data, comparing the fitted SynX against the center values of the window
            for i = 1:self.F0                
                self.Window = self.TS.getWindow(i,self.AnalysisCycles,self.even);
                
                FM = FM_BSA_Class(...
                    self.SignalParams(Fin,:),...
                    self.SignalParams(Fa,:),...
                    self.SignalParams(Ka,:),...
                    1/self.Fs,...
                    real(self.Window.Data),...
                    'verbose',self.verbose,...
                    'debug', self.debug....
                    );
                 
                Synx = zeros(1,FM.nPhases);
                Freq = Synx;
                ROCOF = Synx;
                FM.fig = 1;
                for phase = 1:FM.nPhases
                    FM.phase = phase;
                    if FM.debug,figure(FM.fig);FM.fig=FM.fig+1;end
                    [startPt] = FM.GridSearch(); 
                    if FM.debug,figure(FM.fig);FM.fig=FM.fig+1;end
                    [endpt_BSA] = FM.BSA_Est(startPt);
                    if FM.debug,figure(FM.fig);FM.fig=FM.fig+1;end
                    [estParams] = FM.Param_Est(endpt_BSA);
                    [Synx(phase),Freq(phase),ROCOF(phase)] = FM.Synx_Calc(estParams);
                end
                actFreq(i) = mean(Freq); actROCOF(i)=mean(ROCOF);
                actSynx(:,i) = self.calcSymComp(Synx.');
                
                expSynx(:,i) = self.calcSymComp(self.Window.UserData.Vals.')/sqrt(2);
                expFreq(i) = mean(self.Window.UserData.Freqs(1:3));
                expROCOF(i) = mean(self.Window.UserData.ROCOFs(1:3));
                %disp([i, iter]);
            end
            
            act = struct('Synx',actSynx,'Freq',actFreq,'ROCOF',actROCOF);
            exp = struct('Synx',expSynx,'Freq',expFreq,'ROCOF',expROCOF);
            if bDisplay == true, self.dispErrors(act,exp,self.fig),end
            
            self.verifyEqual(actSynx,expSynx,'AbsTol',1e-6)
            self.verifyEqual(actFreq,expFreq,'RelTol',1e-5)
            self.verifyEqual(actROCOF,expROCOF,'RelTol',1e-5)
            
        end
        
        
    end
    
    
    
    %% --------------------------------------------------------------------   
    % Static method to get indexes into the signalParams matrix
    methods(Static)
        %        returns Xm=1;Fin=2;Ps=3;Fh=4;Ph=5;Kh=6;Fa=7;Ka=8;Fx=9;Kx=10;Rf=11;KaS=12;KxS=13;KfS=14;KrS=15
        %         where variables after Ps are variable output arguments
        function [varargout] = getParamIndex()
            
            for i = 1:nargout
                varargout{i}=i;
            end
        end
        
       %=======================================
        
        function symComp = calcSymComp(x)
            nPhases = size(x,1);
            
            if nPhases == 1
                symComp = x;
                return
            else
                %Calculating symmetrical components
                alfa = exp(2*pi*1i/3);
                Ai = (1/3)*[1 1 1; 1 alfa alfa^2; 1 alfa^2 alfa];
            
                Vabc = x(1:3,:);
                Vzpn = Ai*Vabc; %voltage: zero, positive and negative sequence
                symComp = horzcat(Vabc.',Vzpn(2));
                
                if nPhases < 4, return, end                                       
                Iabc = x(4:6,:);
                Izpn = Ai*Iabc; %curren: zero, positive and negative sequence
            
                symComp = horzcat(symComp,Iabc.',Izpn(2));
            end
        end  
        
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
        
    end
    
    

end