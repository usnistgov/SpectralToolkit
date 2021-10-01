classdef FM_HJS_Class < handle
    %FM_HJS_CLASS A port of the Hooke-Jeeves Student3 C program by
    % Gregory Kyriazis of InMetro
    %
    %  Allen Goldstein, NIST
    %
    %   Major differences between this and the original CVIWindows code :
    %       This has no user interface but is designed to be driven from
    %   scripts, this allowing regression testing 
    %       Also this will not only simulate waveforms (both modulated and
    %   modulating), but also allows for sampled signals to be input
    %
    % Citation:
    %Kyriazis G. A., “Estimating parameters of complex modulated signals from
    %prior information about their arbitrary waveform components,” IEEE Trans.
    %Instrum. Meas., v. 62, no. 6, pp. 1681-1686, June 2013.
    
    % Citation: 
    % Kyriazis G. A., “A Cartesian method to improve the results and
    % save computation time in Bayesian signal analysis,” in Advanced
    % Mathematical and Computational Tools in Metrology and Testing X (AMCTM
    % X), Series on Advances in Mathematics for Applied Sciences, vol. 86, F.
    % Pavese; W. Bremser; A.G. Chunovkina; N. Fischer; A.B. Forbes (eds.),
    % World Scientific, 2015, pp. 229-240.
    %
    %   Detailed explanation goes here
    
    properties 
        Name
        SignalParams
        
        Duration        % signal duration in seconds
        SampleRate
        Sigma
        NumHarm
        
        data            % The modulating signal.  Can be input using the 
                        % GenModData constructor argument.  If not, then a simulated
                        % signal will be generated during configure.  This may also be overwritten later
 
        verbose         % boolean to show status messages
        
    end
    
   properties (Access = private)
       x_NLS
       Tsamp    % sample period
       ino
       ifun     % number of harmonics to fit
       data1
       Gaussrnd
       residue
       
       % Modulating Signal properties
       Ampl_Gen_Mod
       Freq_Gen_Mod
       
       % NLS of the modulating signal
       Freq_NLS
       Mod_NLS
       phi_NLS 
       Result_NLS
       THD
       Residue_NLS
       iterMax_NLS
       epsilon_NLS
       rho_NLS
       
       end_bi_NLS

             
   end
   
% % signalparams  (Note that the labeling convention comes mostly from the standard)
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
   
   
    %% =========================================================================
   % Constructor
    methods
        function obj = FM_HJS_Class(varargin)
            defaultName = 'default FM_HJS';
            defaultSignalParams = [1,50,0,0,0,0,5,5,0,0,0,0,0,0,0]';
            defaultSampleRate = 1/.00006520;  % samples per second
            defaultDuration = 8000/defaultSampleRate;
            defaultSigma = 0.000001;
            defaultNumHarm = 3;
            defaultGenModData = [];
            defaultVerbose = true;
            
            p = inputParser;
            
            validBool = @(x) islogical(x);
            validScalar = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);            
            validateSignalParams = @(x) validateattributes(x,{'double'},{'nrows',15});

            addParameter(p,'Name',defaultName,@ischar);
            addParameter(p,'SampleRate',defaultSampleRate,validScalarPosNum);
            addParameter(p,'Duration',defaultDuration,validScalarPosNum);
            addParameter(p,'SignalParams',defaultSignalParams,validateSignalParams);
            addParameter(p,'Sigma',defaultSigma,validScalar);
            addParameter(p,'NumHarm',defaultNumHarm,validScalarPosNum)
            addParameter(p,'GenModData',defaultGenModData);
            addParameter(p,'Verbose',defaultVerbose,validBool);
            
            parse(p,varargin{:})
            
            obj.Name = p.Results.Name;
            obj.SampleRate = p.Results.SampleRate;
            obj.Duration = p.Results.Duration;
            obj.SignalParams = p.Results.SignalParams;
            obj.Sigma = p.Results.Sigma;
            obj.NumHarm = p.Results.NumHarm;
            obj.data = p.Results.GenModData;
            obj.verbose = p.Results.Verbose;
            
            % set up some global variables
            obj.ino = int32(obj.Duration*obj.SampleRate);          
            
         end
        
    end
    
    %% =========================================================================
    % Public Methods
    methods (Access = public)
        function configure(obj)   
            obj.Tsamp = 1/obj.SampleRate;
            obj.ifun = 2*obj.NumHarm+1;
            [~,~,~,~,~,~,Fa,Ka] = obj.getParamVals(obj.SignalParams);
            obj.Ampl_Gen_Mod = 1;
            %obj.Ampl_Gen_Mod = Ka;         % why is the ampl not the index?
            obj.Freq_Gen_Mod = Fa;
                       
            if isempty(obj.data),obj.genwave();end  
            
            obj.iterMax_NLS = 40;
            obj.epsilon_NLS = 1e-8;
            obj.rho_NLS = 1;
        end
                
        function genwave(obj)
            % generates a simulated modulation signal
            % Note, may be expanded in the future to include triangular and square wave modulations
            obj.Gaussrnd = obj.Sigma*randn(obj.ino,1);            
            t = linspace(0,double(obj.ino)*obj.Tsamp-obj.Tsamp,obj.ino)';
            
             obj.data = obj.Ampl_Gen_Mod*cos(2*pi*obj.Freq_Gen_Mod*t-pi/2);
             obj.data = obj.data+obj.Gaussrnd;   
        end
        
        function mod_Freq_NLS(obj)
            % Least-Squared fit of the modulation signal
            new_bi = zeros(obj.ifun+1,1);
            startpt_NLS = 2*pi*obj.Freq_Gen_Mod*obj.Tsamp;
            omega2 = startpt_NLS;
            GIJ_NLS = obj.f1(omega2);
            new_bi(2:end) = obj.ampli_est(GIJ_NLS,obj.data);
            
            % fitting to the modulation frequency (this is a typical 4-parameter sine fit)
            iters = 0;
            stepLength_NLS = 1;
            while (iters < obj.iterMax_NLS && stepLength_NLS > obj.epsilon_NLS)
                GIJC_NLS = obj.f2(omega2,new_bi);
                end_bi = obj.ampli_est(GIJC_NLS,obj.data);
                omega2 = omega2 + obj.rho_NLS*end_bi(1);
                new_bi(1)=end_bi(1);
                GIJ_NLS = obj.f1(omega2);
                new_bi(2:end) = obj.ampli_est(GIJ_NLS,obj.data);
                iters = iters+1;
                stepLength_NLS = abs(end_bi(1));
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
        
        function mod_Amp_NLS(obj)            
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
                obj.THD = 100*sqrt(Acrms_NLS^2-Mod_Fundamental_NLS^2)/Mod_Fundamental_NLS; % THD in %
            else
                obj.THD = 0;
            end
            

            i = (0:double(obj.ino-1))';
            obj.Result_NLS = ones(obj.ino,1).*obj.Mod_NLS(1);    % DC part
            for j=1:(obj.ifun-1)/2
                obj.Result_NLS = obj.Result_NLS + obj.Mod_NLS(j+1)*cos(j*obj.Freq_NLS*i+obj.phi_NLS(j+1));
            end
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
            a_NLS = zeros((obj.ifun-1)/2,1);
            b_NLS = a_NLS;
            
            if x_NLS == 0, x_NLS=1e-6; end
            a_NLS(1) = bi_NLS(2);
            b_NLS(1) = 0;
            
            k = 2;
            while (k<obj.ifun)
                if (k<(obj.ifun+1)/2)
                    a_NLS(k)=bi_NLS(k+1);
                else
                    b_NLS(k-(obj.ifun-1)/2)=bi_NLS(k+1);
                end
                k=k+1;
            end
            
            % the GIJC model (hypothesis) matrix
            GIJC_NLS = zeros(obj.ino,(obj.ifun));
            GIJC_NLS(:,2)=1;
            i = double(0:obj.ino-1); 
            nfun = (obj.ifun-1)/2;            
            for j = 1:nfun
                GIJC_NLS(:,1)= GIJC_NLS(:,1) + (-a_NLS(j)*i.*sin(j*x_NLS*i)...
                               + b_NLS(j)*i.*cos(j*x_NLS*i))';
            end
            for j = 2:nfun+1
                GIJC_NLS(:,j+1)=cos((j-1)*x_NLS*i);
                GIJC_NLS(:,nfun+j+1) = sin((j-1)*x_NLS*i);
            end
                
                
        end
        
    end
    
    
    %% =========================================================================
   
    % Static Methods
    methods(Static)
        function [varargout] = getParamVals(signalparams)
            varargout = cell(nargout,1);
            for i = 1:nargout
                varargout{i}=signalparams(i,:);
            end
        end
            
        function [varargout] = getParamIndex()
            varargout = cell(nargout,1);
            for i = 1:nargout
                varargout{i}=i;
            end            
        end
                
        function [S] = ampli_est(G,samples)
            % Linear least-squares
            S = (G'*G)\(G'*samples);            
        end

    end
    %% =========================================================================
    
    
    end