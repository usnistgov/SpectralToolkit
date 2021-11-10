classdef FM_NMS_Class < handle
    %FM_NMS_CLASS A port of the Hooke-Jeeves Student3 C program by
    % Gregory Kyriazis of InMetro but using Nelder-Mead direct search instread of Hooke-Jeeves.
    %
    % Citation:
    %Kyriazis G. A., “Estimating parameters of complex modulated signals from
    %prior information about their arbitrary waveform components,” IEEE Trans.
    %Instrum. Meas., v. 62, no. 6, pp. 1681-1686, June 2013.  
    %
    % Citation: 
    % Kyriazis G. A., “A Cartesian method to improve the results and
    % save computation time in Bayesian signal analysis,” in Advanced
    % Mathematical and Computational Tools in Metrology and Testing X (AMCTM
    % X), Series on Advances in Mathematics for Applied Sciences, vol. 86, F.
    % Pavese; W. Bremser; A.G. Chunovkina; N. Fischer; A.B. Forbes (eds.),
    % World Scientific, 2015, pp. 229-240.
    %
    %  Ported by Allen Goldstein, NIST
    %
    %   Major differences between this and the original CVIWindows code :
    %       This has no user interface but is designed to be driven from
    %   scripts, this allowing regression testing 
    %       Also this will not only simulate waveforms (both modulated and
    %   modulating), but also allows for sampled signals to be input
    %
    %
    %   Steps to run a single test case with default parameters and simulated
    %   modulating signal analysed via NLS and modulating signals analysed via BSA    
    %         HJS = FM_HJS_Class()  % Instantiates the class with default properties
    %         HJS.configure()       %        
    %         HJS.mod_Freq_NLS()    %
    %         HJS.mod_Amp_NLS()     %
    %         HJS.plot('NLS')       %
    %         HJS.Freq_BSA()        %
    %         HJS.Ampl_BSA()        %
    %         HJS.plot('BSA')       %      

    
    
    properties 
        
% % SignalParams  (Note that the labeling convention comes mostly from the standard)
%     Xm = signalparams(1,:)*sqrt(2);     % phase amplitude (given by the user in RMS)
%     Fin = signalparams(2,:);    % frequency (must be the same for all 6 channels or an error will be thrown
%     Ps = signalparams(3,:);     % phase 
%     Fh = signalparams(4,:);     % Frequency of the interfering signal
%     Ph = signalparams(5,:);     % Phase of the interfering signal
%     Kh = signalparams(6,:);     % index of the interfering signal    
%     Fa = signalparams(7,:);     % phase (angle) moduation frequency
%     Ka = signalparams(8,:);     % phase (angle) moduation index
%     Fx = signalparams(9,:);     % amplitude moduation frequency
%     Kx = signalparams(10,:);     % amplitude moduation index
%     Rf = signalparams(11,:);     % ROCOF
%     KaS = signalparams(12,:);   % phase (angle) step index
%     KxS = signalparams(13,:);   % magnitude step index
%     KfS = signalparams(14,:);   % frequency step index
%     KrS = signalparams(15,:);   % ROCOF step index (another way to create a frequency ramp)
%        
        SignalParams         % an array of signal parameters as documented above
   
% % NoiseParams A structure as follows:
% struct NoiseParams('NoiseUniformLow',0,'NoiseUniformHi',0,'NoiseGaussMean',0,'NoiseGaussSD',0)
%         NoiseUniformLow % double Uniform distribution noise lowest value
%         NoiseUniformHi  % double Uniform distribution noise highest value
%         NoiseGaussMean  % double noise mean value
%         NoiseGaussSD    % double noise standard deviation
%         RngState = -1   % state of the random number generator.  -1 means no initial state                       
        AmplNoiseParams      % a structure of amplitude noise parameters 
        PhaseNoiseParams     % a structure of phase noise parameters as documented above
        
        % needed info for any signal, simulated or not
        Duration        % signal duration in seconds
        SampleRate      % Samples per Second
        NumHarm
        
        % The modulating signal.  Can be input using the 
        % GenModData constructor argument.  If not, then a simulated
        % signal will be generated during obj.configure().  This may also be overwritten later.        
        data            % The modulating signal (simulated or uploadded   
           
        % The modulated signal.  Can be input using the GenData constructor argument.        
        % if not, then a simulated signal will be generated during obj.configure().  This may be overwritten later.
        data1           
                        
        phaseNoise      % Phase Noise to be added to the phase of the modulating signal, modulated signal, or both 
        amplNoise       % Amplitude Noise to be added
        modNoise        % logical value "True" if noise is to be applied to the modulating signal
        sigNoise        % logical value "True" if noise is to be applied to the modulated signal
        
        
       % Modulating Signal parameters 
       Ampl_Gen_Mod     % Normally set by SignalParams, but can be overwritten after obj.configure() runs
       Freq_Gen_Mod     % Normally set by SignalParams, but can be overwritten after obj.configure() runs
       Phi_Gen_Mod=0      % Defaults to 0 but can be overwritten after obj.configure() runs
       dc_Gen_Mod       % Defaults to 0 but can be overwritten after obj.configure() runs
       THD_NLS
       
       % Modulated Signal parameters
       Ampl_Gen
       Freq_Gen
       Phi_Gen
       dc_Gen
       DeltaF_Gen       % peak frequency deviation
       SSE_BSA       % Sum of Squared Errors of the BSA
       R_Sq_BSA      % Coefficient of determination
       funevals
       
                              
        Name            % Name of the class instance           
        verbose         % logical to show status messages if "true"
        debug             % logical true if in debugging mode
        fig = 1         % keeps track of plot figure numbers
        
       
       Hooke = false;    % use Hooke-Jeeves rather than Nelder-Mead for BSA

        
        %Debugging properties
        hookeContour = {} % cell array, columns: funevals;w0;ka;phi
        
        
    end
    
   properties (Access = private)
       x_NLS
       Tsamp    % sample period
       ino
       ifun     % number of harmonics to fit       
       noiseTs         % a time series containing only the additive noise.
       Gaussrnd
       
       % NLS of the modulating signal
       Freq_NLS
       Mod_NLS
       phi_NLS 
       Result_NLS       
       Residue_NLS
       iterMax_NLS
       epsilon_NLS
       rho_NLS       
       end_bi_NLS
       
       % BSA of the modulated signal
       iterMax_BSA
       epsilon_BSA
       rho_BSA
       endpt_BSA    % the omega result for the BSA       
       sigma
       % stloge
       Phi_BSA
       Modulo_BSA
       Result_BSA
       Residue_BSA
       %zloge        % I do not know what this does, it is always = 0
       %st           % this value never changes either
       %nst
       %phat         % phat was calculated but never used
       
       % BSA model properties
       V_norm           % normalized eiganvector of the model
       invD_norm
       hi       
             
   end   
   
    %% =========================================================================
   % Constructor
    methods
        function obj = FM_NMS_Class(varargin)
            % The constructor accepts name,value pair arguments.  If the
            % argument is not included in the constructor call the default
            % value will ne used.  the arguments and their default values
            % are shown here:
            %
            % Example: FM_NMS_Class(Name1,Value1,Name2,Value2,...NameN,ValueN,)
            %
            % Argument Name    , type  , default value           % comment
            % 'Name'           , char  , 'default FM_NMS'        % 
            % 'SampleRate'     , double,  1/.00006520            %
            % 'Duration'       , double,  8000/defaultSampleRate % signal duration in seconds
            % 'SignalParams    , 15 x 1 array of doubles, [1,49.9876543210,0,0,0,0,4.9876543210,5,0,0,0,0,0,0,0]' % See the SignalParams property description
            % 'PhaseNoiseParams', struct, struct('NoiseUniformLow',0,'NoiseUniformHi',0,'NoiseGaussMean',0,'NoiseGaussSD',0.000001) %
            % 'AmplNoiseParams' , struct, struct('NoiseUniformLow',0,'NoiseUniformHi',0,'NoiseGaussMean',0,'NoiseGaussSD',0.000001) %
            % 'ModNoise'        , logical, 'true'                 % if true, phase noise will be added to the simulated modulating signal
            % 'SigNoise'        , logical, 'true'                 % if true, Amplitude noise will be added to the simulated modulated signal
            % 'NumHarm'         , 3                               % Number of harmonics to use during NLS and BSA analysis
            % 'GenModData'      , []                              % Uploaded modulating signal.  If empty, a simulated signal will be created
            % 'GenData'         , []                              % Uploaded modulated signal.  If empty, a simulated signal will be created
            % 'Verbose'         , 'true                           % If true, computed values will be displayed in the console
                        
            defaultName = 'default FM_NMS';
            defaultSignalParams = [1,49.9876543210,0,0,0,0,4.9876543210,5,0,0,0,0,0,0,0]';
            defaultSampleRate = 1/.00006520;  % samples per second
            defaultDuration = 8000/defaultSampleRate;
            defaultNumHarm = 3;
            defaultGenModData = [];
            defaultGenData = [];
            defaultVerbose = true;
            defaultDebug = false;
            defaultNoiseParams = struct('NoiseUniformLow',0,'NoiseUniformHi',0,'NoiseGaussMean',0,'NoiseGaussSD',0.000001);
            defaultModNoise = true;
            defaultSigNoise = true;
            
            
            p = inputParser;
            
            validBool = @(x) islogical(x);
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);            
            validateSignalParams = @(x) validateattributes(x,{'double'},{'nrows',15});
            validateStruct = @(x) isstruct(x);

            addParameter(p,'Name',defaultName,@ischar)
            addParameter(p,'SampleRate',defaultSampleRate,validScalarPosNum)
            addParameter(p,'Duration',defaultDuration,validScalarPosNum)
            addParameter(p,'SignalParams',defaultSignalParams,validateSignalParams)
            addParameter(p,'PhaseNoiseParams',defaultNoiseParams,validateStruct)
            addParameter(p,'AmplNoiseParams',defaultNoiseParams,validateStruct)
            addParameter(p,'ModNoise',defaultModNoise,validBool)
            addParameter(p,'SigNoise',defaultSigNoise,validBool)            
            addParameter(p,'NumHarm',defaultNumHarm,validScalarPosNum)
            addParameter(p,'GenModData',defaultGenModData)
            addParameter(p,'GenData',defaultGenData)            
            addParameter(p,'Verbose',defaultVerbose,validBool)
            addParameter(p,'Debug', defaultDebug,validBool)
            
            
            parse(p,varargin{:})
            
            obj.Name = p.Results.Name;
            obj.SampleRate = p.Results.SampleRate;
            obj.Duration = p.Results.Duration;
            obj.SignalParams = p.Results.SignalParams;
            obj.PhaseNoiseParams = p.Results.PhaseNoiseParams;
            obj.AmplNoiseParams = p.Results.AmplNoiseParams;
            obj.modNoise = p.Results.ModNoise;
            obj.sigNoise = p.Results.SigNoise;            
            %obj.Sigma = p.Results.Sigma;
            obj.NumHarm = p.Results.NumHarm;
            obj.data = p.Results.GenModData;
            obj.data1 = p.Results.GenData;
            obj.verbose = p.Results.Verbose;
            obj.debug = p.Results.Debug;
            
            % set up some global variables
            obj.ino = int32(obj.Duration*obj.SampleRate);          
            
         end
        
    end
    
    %% =========================================================================
    % Public Methods in external files
    methods (Access = public)
        plot(obj, varargin)
    end

    %% =========================================================================
    % Public Methods
    methods (Access = public)
        function obj = configure(obj) 
            % Configures the sample period and number of harmonics to analyse.
            % Generates random noise properties for later use
            % Configures modulating signal parameters and generates the modulating signal if it was not previously uploaded
            % Configures NLS properties for analysis of the modulating signal
            % Configures modulated signal parameters
            % Generates the modulated signal if it was not previously uploaded
            % Configures BSA properties for analysis of the modulated signal
            
            obj.Tsamp = 1/obj.SampleRate;
            obj.ifun = 2*obj.NumHarm+1;
            [Xm,Fin,Ps,~,~,~,Fa,Ka] = obj.getParamVals(obj.SignalParams);
            
            % generate Phase Noise and Amplitude Noise
            obj.phaseNoise = obj.genNoise(obj.PhaseNoiseParams,obj.ino);
            obj.amplNoise = obj.genNoise(obj.AmplNoiseParams,obj.ino);
            
            % Modualting Signal
            %obj.Ampl_Gen_Mod = 1;
            obj.Ampl_Gen_Mod = Ka;         % why is the ampl not the index?
            obj.Freq_Gen_Mod = Fa;
            obj.dc_Gen_Mod = 0;             % the modulating signal DC offset is not configurable in SignalParams but may be changed here
            %obj.Phi_Gen_Mod = 0;            % the modulating signal initial phase is not configurable in SignalParams but may be changed here
            if isempty(obj.data),obj.genwave();end  
            
            obj.iterMax_NLS = 40;
            obj.epsilon_NLS = 1e-8;
            obj.rho_NLS = 1;
            
            % Modulated Signal
            obj.Ampl_Gen = Xm;
            obj.Freq_Gen = Fin;
            obj.Phi_Gen = Ps*pi/180;
            obj.dc_Gen = 0;
            obj.DeltaF_Gen = Ka*Fa;         % peak frequency deviation = Index * Fmod
            if isempty(obj.data1),obj.gendata();end
            
            % BSA
            obj.iterMax_BSA = 5000;
            obj.epsilon_BSA = 1e-8;
            %obj.rho_BSA = 0.85;              
            %obj.rho_BSA = [0.975,0.85,0.95];  
            obj.rho_BSA = [0.85,0.85,0.85];
            
        end
                
        function genwave(obj)
            % generates a simulated modulation signal
            % Note, may be expanded in the future to include triangular and square wave modulations
            
            
            %obj.Gaussrnd = obj.Sigma*randn(obj.ino,1);            
            t = linspace(0,double(obj.ino)*obj.Tsamp-obj.Tsamp,obj.ino)';
            
             obj.data = obj.Ampl_Gen_Mod*cos(2*pi*obj.Freq_Gen_Mod*t-pi/2);
             if obj.modNoise, obj.data = obj.data+obj.phaseNoise; end   
        end
        
        function obj = mod_Freq_NLS(obj)
            % Least-Squared fit of the modulation signal
            
            new_bi = zeros(obj.ifun+1,1);
            startpt_NLS = 2*pi*obj.Freq_Gen_Mod*obj.Tsamp;
            omega2 = startpt_NLS;
            GIJ_NLS = obj.f1(omega2);
            newbi = obj.ampli_est(GIJ_NLS,obj.data);
            new_bi(1)=0;
            new_bi(2:end) = newbi;
            
            % fitting to the modulation frequency (this is a typical 4-parameter sine fit)
            iters = 0;
            stepLength_NLS = 1;
            while (iters < obj.iterMax_NLS && stepLength_NLS > obj.epsilon_NLS)
                GIJC_NLS = obj.f2(omega2,new_bi);
                end_bi = obj.ampli_est(GIJC_NLS,obj.data);
                omega2 = omega2 + obj.rho_NLS*end_bi(1);
                GIJ_NLS = obj.f1(omega2);
                newbi = obj.ampli_est(GIJ_NLS,obj.data);
                new_bi(1)=end_bi(1);
                new_bi(2:end)=newbi;
                iters = iters+1;
                stepLength_NLS = abs(end_bi(1));
            end
            if iters == obj.iterMax_NLS
                warning('mod_Freq_NLS did not converge after %d iterations',iters)
            end
            
            obj.end_bi_NLS=new_bi;
            
            obj.Freq_NLS = omega2;
            
            % display the frequencies if verbose
            if obj.verbose
                nfun = (obj.ifun-1)/2;
                %fList = zeros(nFun,1);
                fList = '';
                for i=1:nfun+1
                    fItem = sprintf('%d: %f \n',(i-1),(i-1)*omega2);
                    fList = sprintf('%s %s',fList,fItem);
                end
                fprintf('Frequency (rad/unit):\n%s',fList);
            end
            
            
        end
        
        function obj = mod_Amp_NLS(obj)   
            % Using the results calculated by mod_Freq_NLS(), 
            % calculates the amplitude and phase of the modulating signal
            
            % intialize the coeffs with the results from mod_Freq_NLS
            a_NLS = zeros((obj.ifun+1)/2,1);
            b_NLS = a_NLS;
            
            a_NLS(1) = obj.end_bi_NLS(2);
            b_NLS(1) = 0;
            k=1;
            while (k<obj.ifun+1)
                if (k<(obj.ifun+1)/2)
                    a_NLS(k+1)=obj.end_bi_NLS(k+2);
                else
                    b_NLS(k-(obj.ifun-1)/2)=obj.end_bi_NLS(k+1);
                end
                k=k+1;
            end
            
            c_NLS = complex(a_NLS,-b_NLS);
            obj.Mod_NLS = abs(c_NLS);
            obj.phi_NLS = angle(c_NLS); obj.phi_NLS(1)=0;           
            phi_corr_NLS = zeros((obj.ifun+1)/2,1);
            
            % RMS
            a_NLS = a_NLS/sqrt(2);
            b_NLS = b_NLS/sqrt(2);
            a_corr_NLS = a_NLS.*cos(phi_corr_NLS)+b_NLS.*sin(phi_corr_NLS);
            b_corr_NLS = b_NLS.*cos(phi_corr_NLS)-a_NLS.*sin(phi_corr_NLS);
            
            c_corr_NLS = complex(a_corr_NLS,b_corr_NLS);  % complex
            Acrms_NLS = sum(abs(c_corr_NLS).^2);              % sum of squares
            Mod_Fundamental_NLS = abs(c_corr_NLS(2));     % Fundamental amplitude
            Mod_corr_NLS = 100*abs(c_corr_NLS)/Mod_Fundamental_NLS;
            phi_ini = obj.phi_NLS(2)*180/pi;                  % phase in degrees
            phi_corr_NLS = obj.phi_NLS * 180/pi + (1:(obj.ifun+1)/2)'*(-90-phi_ini);
            phi_corr_NLS(1)=0;  % DC has no phase
            
            Acrms_NLS = sqrt(Acrms_NLS);
            if Acrms_NLS^2 > Mod_Fundamental_NLS^2
                obj.THD_NLS = 100*sqrt(Acrms_NLS^2-Mod_Fundamental_NLS^2)/Mod_Fundamental_NLS; % THD in %
            else
                obj.THD_NLS = 0;
            end
            

            i = (0:double(obj.ino-1))';
            obj.Result_NLS = ones(obj.ino,1).*obj.Mod_NLS(1);    % DC part
            for j=1:(obj.ifun-1)/2
                obj.Result_NLS = obj.Result_NLS + obj.Mod_NLS(j+1)*cos(j*obj.Freq_NLS*i+obj.phi_NLS(j+1));
            end
            obj.Residue_NLS = obj.Result_NLS-obj.data;
            % if verbose is true, display amplitudes, phases 
            if obj.verbose
                nfun = (obj.ifun-1)/2;
                fList = '';
                for i=1:nfun+1
                    fItem = sprintf('%d: %f \n',i-1,Mod_corr_NLS(i));
                    fList = sprintf('%s %s',fList,fItem);
                end
                fprintf('NLS Amplitudes:\n%s',fList);
                
                fList = '';
                for i=1:nfun+1
                    fItem = sprintf('%d: %f \n',i-1,phi_corr_NLS(i));
                    fList = sprintf('%s %s',fList,fItem);
                end
                fprintf('NLS Phases (rad/unit):\n%s',fList);
                
            end
        end
        
        function gendata(obj) 
            % Generates the modulated signal using the results from the NLS of the modulating signal
                                    
            i = double(0:obj.ino-1)';
            
            % generate the modulating signal (a very convoluted method)
            Onda_Seno = (obj.Ampl_Gen_Mod/obj.Freq_Gen_Mod)*sin(2*pi*obj.Freq_Gen_Mod*i*obj.Tsamp+obj.Phi_Gen_Mod);
            if obj.modNoise, Onda_Seno = Onda_Seno + obj.phaseNoise; end
            
            obj.data1 = sqrt(2)*obj.Ampl_Gen*cos((2*pi*obj.Freq_Gen + 2*pi*(obj.DeltaF_Gen/obj.Ampl_Gen_Mod)*obj.dc_Gen_Mod)*i*obj.Tsamp...
                                                 + obj.Phi_Gen...
                                                 + (obj.DeltaF_Gen/obj.Ampl_Gen_Mod)*Onda_Seno);                                             
             if obj.sigNoise, obj.data1 = obj.data1 + obj.amplNoise; end
            
        end
        
        function Freq_BSA(obj)
            % BSA of the modulated signal using results of the NLS of the modulated signal
            %
            % we need to find a good startpoint:
            %   startpt(1) = best guess of the carrier frequency
            %   startpt(2) = best guess of the modulation phase (Phi_Gen_Mod)
            %   startpt(3) = best guess of Delta Frequency (DeltaF_Gen)
            % for the modulated signal, we usually have a pretty good idea
            % of the carrier frequency, however, for the modulated signal, the
            % final modulating frequency phase and delta frequency, are not the
            % same values as those used to generate.
            
            % Observation of the objective function contour shows us that a
            % grid of 20 points across 0 to 2 pi and from 0 to 4DF
            % will have at least one or two good initial guesses up to DF of about 50.  
            % This would be 400 function evals if allowed to search all points.
            % But we also know that there are no local minima, so once we
            % reach a threshold of the objective function, we know our
            % start value will be good enough.  Typically we get less than 36 fevals
            
            % look for a difference between zBest and zWorst to stop searching
            
            obj.funevals = 0;
            startpt(1) = 2*pi*obj.Freq_Gen*obj.Tsamp;     % Carrier frequency is assumed to be good enough
            grid = 20;      % This is the grid search resolution, number of points per axis.
            thresh = 3000;   % this is the grid search threshold:  difference between worst and best before stopping the search
            OMEGA2 = linspace(0,2*pi,grid);
            OMEGA3 = linspace(0,4*2*pi*obj.DeltaF_Gen*obj.Tsamp,grid);
            z = zeros(12,10);
            zWorst = obj.f([startpt(1),OMEGA2(1),OMEGA3(1)]);
            zBest = zWorst;
            for i = 1:grid
                for j = 1:grid
                    z(i,j) = obj.f([startpt(1),OMEGA2(i),OMEGA3(j)]);
                    if z(i,j) < zBest, zBest = z(i,j); end
                    if z(i,j) > zWorst, zWorst = z(i,j); end 
                    if abs(zWorst-zBest) > thresh,break,end                                            
                end
                if abs(zWorst-zBest) > thresh,break,end 
            end
            startpt(2) = OMEGA2(i);
            startpt(3) = OMEGA3(j);
            
            %  vebose display ----------------------------
            if obj.verbose
                fprintf('Grid Search: Phase start = %e, DF start = %e, funevals = %d\n',startpt(2),startpt(3),obj.funevals)                
            end
            % ---------------------------------------------
            % debugging display of contour plots -----------
            if obj.debug
                figure(obj.fig);obj.fig=obj.fig+1;
                dF = 2*pi*obj.DeltaF_Gen*obj.Tsamp;
                obj.plot('fcontour3',[startpt(1),startpt(1);0,2*pi;0,2*dF],30)
                hold on
                for k = 1:i
                    for l = 1:j
                        plot3(OMEGA2(k),OMEGA3(l),z(k,l),'.')
                    end
                end
                hold off
            end
            % ----------------------------------------------------
                                                             
            nvars = length(startpt);
            obj.funevals = 0;
            
            
            if ~obj.Hooke
                % ------Replacing Hooke-Jeeves with Nelder Mead(NM)
                %[iters,obj.endpt_BSA] = obj.hooke(nvars,startpt, obj.rho_BSA, obj.epsilon_BSA, obj.iterMax_BSA);
                
                % Using fminsearch(problem)where problem is a structure.  see MATLAB "doc fminsearch"
                opts = optimset(@fminsearch);    % configures NM options with default values
                opts.TolX = 1e-8;
                opts.TolFun = 1e8;
                PROBLEM = struct('objective',[],'X0',[],'options',opts,'solver','fminsearch');
                X0 = startpt;
                PROBLEM.X0 = X0;
                PROBLEM.objective = @(X0) obj.f(X0);
                %PROBLEM.objective = @(X0) obj.fopt1(omega1,X0);         % handle to the objective function
                
                if obj.debug
                    %figure(obj.fig);obj.fig=obj.fig+1;                    
                    %PROBLEM.options.PlotFcns = 'optimplotfval';  % can plot val or x but not both
                    %PROBLEM.options.PlotFcns = 'optimplotx';     % bar chart of x but cannot also plot val
                    %PROBLEM.options.Display = 'iter'             % display info on a per-iteration basis
                    PROBLEM.options.OutputFcn = @obj.out;         % plots the contour map and points
                    [obj.endpt_BSA,fval,exitflag,outpoint] = fminsearch(PROBLEM);
                    hold off
                    %pause
                else
                    [obj.endpt_BSA] = fminsearch(PROBLEM);
                end
                % --------------------------------------------------------
            else
                [iters,obj.endpt_BSA] = obj.hooke(nvars,startpt, obj.rho_BSA, obj.epsilon_BSA, obj.iterMax_BSA);
                % display the result
                if obj.verbose
                    fprintf('\n==========\nTotal iterations: %d\n',iters)
                    nfun = (obj.ifun-1)/2;
                    fList = '';
                    for k = 1:nfun
                        fItem = sprintf('%d: %f \n',k,obj.endpt_BSA(k));
                        fList = sprintf('%s %s',fList,fItem);
                    end
                    fprintf('BSA Frequencies (rad/unit):\n%s',fList);
                end
                
            end
           
            
        end
        
        function Ampl_BSA(obj)
            % Calculates the Amplitude and Phase of the modulated signal after Freq_BSA() has been run
            % Generated the "best fit" modulated signal and the residual.
            %
            % if the modulated signal is uploaded and no modulating signal
            % analysis is performed, the following properties must first be
            % directly written by the calling script:
            % obj.Mod_NLS, obj.Freq_NLS, obj.phi_NLS
            
            
            % preallocate vectors a and b
            nfun = (obj.ifun-1)/2;
            a = zeros((nfun+1)/2,1);
            b = zeros((nfun+1)/2,1);
            
            B = obj.V_norm*obj.invD_norm;
            bi = B*obj.hi;
            a(1) = bi(1);
            %b(1) = 0;   % not needed due to preallocation  
            k = 2;
            while k <= nfun
                if k < (nfun+1)/2+1
                    a(k) = bi(k);
                else
                    b(k-(nfun-1)/2) = bi(k); 
                end
                k = k+1;
            end  
            cBSA = complex(a,b);        % complex
            obj.Modulo_BSA = abs(cBSA);
            obj.Phi_BSA = angle(cBSA);
            
            %a_corr = zeros((nfun+1)/2,1);    % preallocation not needed
            %b_corr = zeros((nfun+1)/2,1);
            % Modulo_corr = zeros((nfun+1)/2,1)
            %Phi_corr = zeros((nfun+1)/2,1);
            Acrms = 0;
            a = a/sqrt(2);
            b = b/sqrt(2);
            a_corr = a.*cos(obj.Phi_BSA)+b.*sin(obj.Phi_BSA);
            b_corr = b.*cos(obj.Phi_BSA)-a.*sin(obj.Phi_BSA);
            cCorr = complex(a_corr,b_corr);
            Acrms = Acrms + sum(abs(cCorr));
            Mod_Fundamental = abs(cCorr(2));
            
            if obj.verbose
                %Modulo_corr(1)=100*obj.Modulo_BSA(1)/Mod_Fundamental;
                %Modulo_corr(2)=100*abs(cCorr(2:end)/Mod_Fundamental);
                Modulo_corr = 100*abs(cCorr)/Mod_Fundamental;
                Phi_corr = angle(cCorr)*180/pi;
                fList = '';
                for k = 1:nfun-1
                    fItem = sprintf('%d: %f \n',k,Modulo_corr(k));
                    fList = sprintf('%s %s',fList,fItem);                    
                end
                fprintf('BSA Amplitudes:\n%s',fList);
                fList = '';
               for k = 1:nfun-1
                    fItem = sprintf('%d: %f \n',k,Phi_corr(k));
                    fList = sprintf('%s %s',fList,fItem);                    
                end
                fprintf('BSA Phases (deg):\n%s',fList);
                
            end
            
            % Modulating signal
            Onda_Result = zeros(obj.ino,1);
            i = (0:double(obj.ino)-1)';
            k = 1;
            while k < (obj.ifun+1)/2
                Onda_Result = Onda_Result + (obj.Mod_NLS(k+1)/(k*obj.Freq_NLS))*sin(k*((obj.Freq_NLS*i)+obj.endpt_BSA(2))+obj.phi_NLS(k+1));
                k = k+1;
            end
            
%             % A convoluted method of finding the residue after finding the modulating signal but before finding the modulated signal
%             % (This was shown to yield the same results as siply subtracting obj.data1 from obj.Result.BSA so is replaced below)
%             obj.Residue_BSA = obj.Modulo_BSA(1)+obj.Modulo_BSA(2)...
%                                *cos(...
%                                     (obj.endpt_BSA(1)...
%                                      +(obj.endpt_BSA(3)/obj.Mod_NLS(2))*obj.Mod_NLS(1))*i...
%                                    +((obj.endpt_BSA(3)/obj.Mod_NLS(2))*Onda_Result+obj.Phi_BSA(2)))...
%                               - obj.data1; 
                         
           % Now we find the modulated signal           
           Onda_Result = Onda_Result * obj.endpt_BSA(3)/obj.Mod_NLS(2);
           Onda_Result = Onda_Result + (obj.endpt_BSA(1) +  (obj.endpt_BSA(3)/obj.Mod_NLS(2))*obj.Mod_NLS(1))*i + obj.Phi_BSA(2);
           %Onda_Result_Temp = cos(Onda_Result);     % why not calculate this in place?
           %Onda_Result = Onda_Result_Temp;          % not needed, just calculate in place!
           Onda_Result = cos(Onda_Result);      % The above, only calculated in place
           Onda_Result = Onda_Result * obj.Modulo_BSA(2);   % first harmonic amplitude 
           Onda_Result = Onda_Result + obj.Modulo_BSA(1);   % DC offset
           obj.Result_BSA = Onda_Result;       
           
           % A simpler ,ethod of findint the residue
           obj.Residue_BSA = obj.Result_BSA - obj.data1;
           
           % some goodness of fit parameters
           obj.SSE_BSA = sum(obj.Residue_BSA.^2);
           obj.R_Sq_BSA = 1 - sum((obj.data1 - obj.Result_BSA).^2)/sum((obj.data1 - mean(obj.data1)).^2);
           
           
        end
                   
    end
    

    %% =========================================================================
    methods (Access = private)
                
        function [GIJ_NLS] = f1(obj,x_NLS)                        
            % f1 (initial estimates for 4-parameter sine fit)
            
            if x_NLS == 0, x_NLS=1e-6; end
            GIJ_NLS = zeros(obj.ino,obj.ifun);
            GIJ_NLS(:,1) = 1;       % DC
            i = double(0:obj.ino-1);
            nfun = (obj.ifun-1)/2;
            for j=1:nfun
                GIJ_NLS(:,j+1)=cos(j*x_NLS*i);
                GIJ_NLS(:,nfun+j+1) = sin(j*x_NLS*i);                
            end                        
        end 
        
        function [GIJC_NLS] = f2(obj,x_NLS,bi_NLS)
            % function to find the frequency of the modulating signal
            % (This is basically the iterative part of a 4-parameter sine fit)
            a_NLS = zeros((obj.ifun-1)/2,1);
            b_NLS = a_NLS;
            
            if x_NLS == 0, x_NLS=1e-6; end
            a_NLS(1) = bi_NLS(2);       % DC
            b_NLS(1) = 0;
            
            k = 2;
            while (k<obj.ifun)
                if (k<(obj.ifun+1)/2)
                    a_NLS(k)=bi_NLS(k+1);  % cosines
                else
                    b_NLS(k-(obj.ifun-1)/2)=bi_NLS(k+1);  % sines
                end
                k=k+1;
            end
            
            % the GIJC model (hypothesis) matrix
            GIJC_NLS = zeros(obj.ino,(obj.ifun));
            GIJC_NLS(:,2)=1;
            i = double(0:obj.ino-1); 
            nfun = (obj.ifun-1)/2;            
            for j = 1:nfun
                GIJC_NLS(:,1)= GIJC_NLS(:,1) + (-a_NLS(j)*j*i.*sin(j*x_NLS*i)...
                               + b_NLS(j)*j*i.*cos(j*x_NLS*i))';
            end
            for j = 2:nfun+1
                GIJC_NLS(:,j+1)=cos((j-1)*x_NLS*i);
                GIJC_NLS(:,nfun+j+1) = sin((j-1)*x_NLS*i);
            end                
        end
        
        
        
        function [iters,endpt] = hooke(obj, nvars, startpt, rho, epsilon, itermax)
            % A Hooke-Jeeves pattern search function
            %
            [newx, xbefore] = deal(startpt);
            delta = abs(startpt .* rho);
            % replace any 0's with rho
            %delta(delta==0) = rho;
            idx = find(delta==0);       % experiment to use a array of rho values
            delta(idx)=rho(idx);
            
            iadj = 0;
            %steplength;
            steplength = rho(2);
            newf = 0;
            iters = 0;
            [fbefore] = obj.f(newx);      % first call to the objective function
            newf = fbefore;               % newf will be changed by best_nearby
            while ((iters < itermax) && (steplength > epsilon))
                iters = iters+1;
                iadj = iadj+1;
               
                % print the intermediate values
                for i=1:nvars
                    if obj.verbose
                        fList = '';
                        for k = 1:nvars
                            fItem = sprintf('%d: %f \n',k,xbefore(k));
                            fList = sprintf('%s %s',fList,fItem);                            
                        end
                        fprintf("\nAfter %d funevals, f(x)=%1.4e at\n%s",obj.funevals,fbefore,fList);  
                        %fprintf("steplength = %e, df(x) = %1.10e\n",steplength, abs(newf-fbefore) )                      
                    end
                    
                    % save the intermediate values for later plotting
                    if obj.debug
                        element = num2cell(horzcat(obj.funevals,fbefore,xbefore));
                        obj.hookeContour=[obj.hookeContour;element];                        
                    end
                    
                    % find the best new point, one coordinate at a time
                    newx = xbefore;
                    [newf, newx] = obj.best_nearby(delta, newx, fbefore, nvars);
                    
                    % if we made an improvement, persue that direction
                    keep = 1;
                    while newf < fbefore && keep == 1
                        iadj = 0;
                        for i = 1 : nvars
                            % determine the (sign) direction of detla
                            if newx(i) <= xbefore(i)
                                delta(i) = 0 - abs(delta(i));
                            else
                                delta(i) = abs(delta(i));
                            end
                            % move in this direction
                            tmp = xbefore(i);
                            xbefore(i) = newx(i);
                            newx(i) = newx(i)*2-tmp;
                        end
                        fbefore = newf;
                        [newf, newx] = obj.best_nearby(delta, newx, fbefore, nvars);
                        % if we got worse, end the search (may be able to
                        % later check for local minima with trial
                        % checks in both directions, this of course would take time)
                        if newf >= fbefore
                            break
                        end
                        % make sure that the differences between the new
                        % and old points are due to actual displacements;
                        % bewre of roundoff errors that might cause newf <
                        % fbefore
                        keep = 0;
                        for i = 1 : nvars
                            keep = 1;
                            if abs(newx(i) - xbefore(i)) > 0.5 * abs(delta(i))
                                break
                            else
                                keep = 0;
                            end
                        end
                    end
                     if steplength >= epsilon && newf >= fbefore
                         steplength = steplength * rho(2);
                         %steplength = abs(newf-fbefore);
                         delta = delta .* rho;
                     end

                end
                 if obj.debug
                    fprintf("steplength = %e, df(x) = %1.10e\n",steplength, abs(newf-fbefore) )
                 end
            end
            endpt = xbefore;
        end
         
         function [minf,point] = best_nearby(obj, delta, point, prebest, nvars)
             % function called by Hooke-jeeves to find the best point
             f = @(z) obj.f(z);         % function handle
             minf = prebest;
             z = point;
             for i = 1: nvars
                z(i) = point(i) + delta(i);
                ftmp = f(z);        % function called by it's handle
                if ftmp < minf
                    minf = ftmp;
                else
                    delta(i) = 0 - delta(i);
                    z(i) = point(i) + delta(i);
                    ftmp = f(z);
                    if ftmp < minf
                       minf = ftmp;
                    else
                        z(i) = point(i);
                    end
                end
             end
             point = z;
         end
         
         function [y] = fopt1(obj,omega1,X0)
             x = [omega1,X0];
             y = obj.f(x);
         end
% This can be deleted         
%          function [y] = f_PD(obj,w1,w2,w3)
%              x =[w1,w2,w3];
%              y = obj.f(x);
%          end             
                  
         function [y] =f(obj,x)
             obj.funevals = obj.funevals+1;
             omega = x;
             % obj.zloge = 0;        % Commented out because it never changes
             % nfun = (obj.ifun+1)/2;
             Onda = zeros(obj.ino,1);
             i = double(0:obj.ino-1)';
             
             % set GIJ
             k=1;
             while k<(obj.ifun+1)/2
                 Onda = Onda + (obj.Mod_NLS(k+1)/(k*obj.Freq_NLS))*sin(k*(obj.Freq_NLS*i+omega(2))+obj.phi_NLS(k+1));
                k = k+1;
             end
             GIJ = zeros(obj.ino,3);
             GIJ(:,1) = 1;      % DC
             GIJ(:,2) = cos((omega(1)+(omega(3)/obj.Mod_NLS(2))*obj.Mod_NLS(1))*i...
                            +(omega(3)/obj.Mod_NLS(2))*Onda);   
             GIJ(:,3) = sin((omega(1)+(omega(3)/obj.Mod_NLS(2))*obj.Mod_NLS(1))*i...
                            +(omega(3)/obj.Mod_NLS(2)).*Onda); 
             stloge = obj.prob(GIJ);
             y = -stloge;
             %obj.nst = -obj.st;    % Commented out because it does nothng
                                   
         end
         
         function stloge = prob(obj,GIJ)
            iNo = double(obj.ino);      % makes the below more readable
            HIJ = obj.ortho(GIJ);
             nfun = size(HIJ,2);
             obj.hi = zeros(nfun,1);
             h2=0;
             for j=1:nfun
                 h1 = sum(obj.data1.*HIJ(:,j));
                 obj.hi(j)=h1;
                 %h2 = h2 + h1.*h1;
             end
             h2 = sum(obj.hi.^2);
             h2bar = h2/nfun;
             
             y2 = sum(obj.data1.^2);
             y2bar = y2/iNo;
             
             qq = 1-h2/y2;
             if qq<=0
                qq = 1e-16; 
             end

             stloge = log(qq)*(nfun-iNo)/2;
             dif = y2bar-nfun*h2bar/iNo;
             obj.sigma = sqrt(abs(dif)*iNo/(iNo-nfun-2));

% The commented code below does nothing because obj.zloge is always 0             
%             ahold = stloge-obj.zloge;
             
             % this code is useless, obj.zloge is always 0
%              if obj.zloge == 0
%                  obj.st = 0;
%              else
%                  obj.st = exp(ahold);
%              end
             
%             obj.phat = nfun*(obj.sigma*obj.sigma+h2bar)*obj.st;
             
         end
         
         function [HIJ] = ortho(obj,GIJ)
             % Orthogonalizes the hypothesis function
             M = GIJ'*GIJ;
             M = (M+M.')/2; % The matlab eig function must recognize the matrix as symetrical
             [V,D_vector] = eig(M,'vector');         % eiganvalues and eiganvectors
             SqrSumCol = sum(V.^2);
             norm = sqrt(SqrSumCol);
             obj.V_norm = V./norm;
             D_vec = sqrt(abs(D_vector));
             D_norm = diag(D_vec);
             obj.invD_norm = inv(D_norm);
             A = GIJ*obj.V_norm;
             HIJ = A*obj.invD_norm; 
             %HIJ = A\D_norm;         % may not need to calculate the inverse of D_norm, needs to be tested
         end                                            
        
         
         
    end
    
    %% =========================================================================
    % Debugging tools
   % draw the contour of the optimazation problems 
    methods (Access = public)
        
        function stop = out(obj,x, optimValue, state)
            % output function for ploting the trajectory during NM minimization
            stop = false;
            switch state
                case 'init'
                    hold off
                    f = 2*pi*50*obj.Tsamp;
                    f1 = 2*pi*45*obj.Tsamp;
                    f2 = 2*pi*55*obj.Tsamp;
                    dF = 2*pi*obj.DeltaF_Gen*obj.Tsamp;
                    figure(obj.fig);obj.fig=obj.fig+1;
                    obj.plot('fcontour3',[f,f;0,2*pi;0,2*dF],30)                    
                    %obj.plot('fcontour3',[f1,f2;-pi,pi;dF,dF],100)
                    hold on;
                case 'iter'
                    %plot3(x(1),x(2),optimValue.fval,'.','MarkerSize',optimValue.iteration + 1);
                    %plot3(x(2),x(3),optimValue.fval,'.','MarkerSize',optimValue.iteration + 1);
                    plot3(abs(x(2)),abs(x(3)),optimValue.fval,'.');
                    drawnow
            end
        end
% replaced by obj.plot('fcontour3')            
%             function fContourPD(obj)
%             % Objective function contour for carrier phase and delta-freq
%             
%             omega1 = 2.047812652508917e-02;
%             omega2 = 1.570796319854123e0;
%             omega3 = 1.021630416912134e-02; 
% %             x = [omega1,omega2,omega3];
% %             y = obj.f(x);
%             
%             iter = 100;
%             y = zeros(iter,iter);
%             dOmega1 = linspace(0,omega1*2,iter);
%             dOmega2 = linspace(0,omega2*2,iter);
%             dOmega3 = linspace(0,omega3*2,iter);
%             count=0;
%             wb=waitbar(count,'wait');
%             for i = 1:iter
%                 for j = 1:iter
%                     y(i,j) = obj.f([omega1,dOmega2(i),dOmega3(j)]);
%                     count = count+1;
%                     waitbar(count/iter^2);
%                 end
%             end
%             close(wb)
%             contour3(dOmega2,dOmega3,y,50)
%             colormap(hsv)           
%             xlabel('Carrier Phase')
%             ylabel('Delta Freq')
%             colorbar
%             
%        end
    end
    
    %% =========================================================================
   
    % Static Methods
    methods(Static)
        % gets the parameter values contained in SignalParams
        function [varargout] = getParamVals(signalparams)
            varargout = cell(nargout,1);
            for i = 1:nargout
                varargout{i}=signalparams(i,:);
            end
        end
            
        function [varargout] = getParamIndex()
            % gets the indexes into SignalParams to later set or get the values
            varargout = cell(nargout,1);
            for i = 1:nargout
                varargout{i}=i;
            end            
        end
                
        function [S] = ampli_est(G,samples)
            % Linear least-squares function using MATLAB mldivide operator
             S = (G'*G)\(G'*samples);  
        end
        
        function [noise] = genNoise(NoiseParams,length)
           % Generates random noise based on the structure elements.  Generally either uniform
           % or Gaussian noise will be generated by setting some parameters to 0, but they can be combiled
           %
           % NoiseParams.NoiseUniformHi - NoiseParams.NoiseUniformLow is a range of uniform noise
           % NoiseParams.NoiseGaussSD is the standard deviation for Gaussian Noise
           % NoiseParams.NoiseGaussMean offsets the noise
           
           noise = zeros(length,1);
           noise = noise + (NoiseParams.NoiseUniformHi - NoiseParams.NoiseUniformLow)*randn(length,1);
           noise = noise + NoiseParams.NoiseGaussSD*randn(length,1);
           noise = noise + NoiseParams.NoiseGaussMean;
        end
 
    end
    %% =========================================================================
    
    
end   
