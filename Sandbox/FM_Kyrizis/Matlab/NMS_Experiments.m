classdef NMS_Experiments < matlab.unittest.TestCase
    % Class-based unit testing for Spectral Toolkit
%
%   run these tests with two command-line commands:
%   - testCase = NMS_Experiments;
%   - res = testCase.run OR res = run(testCase);
%
    
    properties
        NMS     % a single instance of a FM_HJS class

    end
    
    %% Constructor ========================================================
    methods 
        function obj = NMS_Experiments(obj)
            obj.NMS = FM_NMS_Class();
        end
    end
    
   %% Test methods ======================================================== 
    methods (Test)
        % These functions will be automatically run at the command line command:
        % res = testCase.run OR res = run(testCase);
        function Tests(obj)
            %testDefaultNMS(obj);
            testModIndexRange(obj)
            %troubleshootLocalMax(obj)
            %plotKContours(obj)
            %plotPContours(obj)
        end
        

    end
    
    %% ====================================================================
    methods (Access = private)
        % First test will check the default NMS settings
        function obj = testDefaultNMS(obj)
            obj.NMS = FM_NMS_Class();
            obj.runOne()
        end
        
        function obj = testModIndexRange(obj)
            % ranges across modulation indexes from 0.1 to 10 in 0.5 increments.  With noise
            % As of this update, the Freq_BSA is not finding correct answers for Ka = 2.6,3.1,3.6,6.1,7.1,8.1 
           %obj.NMS.modNoise = false;    % remove amplitude noise
           %obj.NMS.phaseNoise = false;   % remove phase noise

            obj.NMS.verbose = false;
            obj.NMS.debug = true;
            Kincr = 0.5;
            Kstart = 0.1;
            Kend = 10;
            [~,Fin,~,~,~,~,Fa,Ka] = obj.NMS.getParamIndex();
            obj.NMS.SignalParams(Fin) = 50.0;
            obj.NMS.SignalParams(Fa) = 5.0;
            k = Kstart;
            while k <= Kend
                obj.NMS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
                %obj.NMS.fig = 1;
                obj.NMS.SignalParams(Ka) = k;
                obj.runOne()  
                obj.NMS.plot('BSA')
                
                k = k + Kincr;
                %pause
            end
        end
        
        function obj = troubleshootLocalMax(obj)
            %obj.NMS = FM_NMS_Class();
            % 50f0_5Fa0_5Ka0 with no additive noise
            [~,Fin,~,~,~,~,Fa,Ka] = obj.NMS.getParamIndex();
            obj.NMS.SignalParams(Fin) = 50.0;
            obj.NMS.SignalParams(Fa) = 5.0;
            

           %obj.NMS.fig = 1;             % reset the figure number 
           %obj.NMS.hookeContour = {}; % clear the contour from last run
           obj.NMS.verbose = false;   % for now verbose is not working
           obj.NMS.debug = true;      % set the debug mode
           obj.NMS.modNoise = false;    % remove amplitude noise
           obj.NMS.phaseNoise = false;   % remove phase noise
           
           %k = 4.1;
           k = 10.0;
           obj.NMS.SignalParams(Ka) = k;
           obj.NMS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
           obj.runOne()
           obj.NMS.plot('BSA')           
           %obj.NMS.plot('hookeContour');
        end
        
        function obj = plotKContours(obj)
            [~,Fin,Ps,~,~,~,Fa,Ka] = obj.NMS.getParamIndex();
            obj.NMS.SignalParams(Fin) = 50.0;
            obj.NMS.SignalParams(Fa) = 5.0;
            fCarr = 2*pi*50/obj.NMS.SampleRate;
            dF = 2*pi*25/obj.NMS.SampleRate;
            
            kStart = .1;
            kEnd = 5;
            kIncr = 0.5;
            k = kStart;
            while k <= kEnd
                obj.NMS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
                %obj.NMS.fig = 1;
                obj.NMS.SignalParams(Ka) = k;
                
                % run one iteraton
                obj.NMS.data = [];      % clear the modulating signal data from the preceeding test
                obj.NMS.data1 = [];      % clear the modulated signal data from the preceeding test
                obj.NMS.configure();
                obj.NMS.mod_Freq_NLS()
                obj.NMS.mod_Amp_NLS()
                
                obj.NMS.plot('fcontour3',[fCarr,fCarr;-pi,pi;0,2*dF],30)
                
                
                k = k+kIncr;
                %pause
            end
            
                
        end
        
        function plotPContours(obj)
            [~,Fin,Ps,~,~,~,Fa,Ka] = obj.NMS.getParamIndex();
            obj.NMS.SignalParams(Fin) = 50.0;
            obj.NMS.SignalParams(Fa) = 5.0;
            fCarr = 2*pi*50/obj.NMS.SampleRate;
            dF = 2*pi*25/obj.NMS.SampleRate;
            pStart = -pi;
            pEnd = pi;
            pIncr = pi/10;
            p = pStart;
           while p <= pEnd
                %obj.NMS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
                %obj.NMS.fig = 1;
                %obj.NMS.SignalParams(Ps) = p*180/pi;
                %disp(obj.NMS.SignalParams(Ps))
                obj.NMS.Phi_Gen_Mod = p;
                
                % run one iteraton
                obj.NMS.data = [];      % clear the modulating signal data from the preceeding test
                obj.NMS.data1 = [];      % clear the modulated signal data from the preceeding test
                obj.NMS.configure();
                obj.NMS.mod_Freq_NLS()
                obj.NMS.mod_Amp_NLS()
                
                figure(1)
                obj.NMS.plot('fcontour3',[fCarr,fCarr;-pi,pi;0,2*dF],30)
                
%                 figure(2)
%                 plot(obj.NMS.data1)
                
                
                p = p+pIncr;
                disp(p)
                %pause
           end    
        end
            
            
       
        
        function runOne(obj)
            % runs one simulation of NMS
            obj.NMS.data = [];      % clear the modulating signal data from the preceeding test
            obj.NMS.data1 = [];      % clear the modulated signal data from the preceeding test
            obj.NMS.configure();
            obj.NMS.mod_Freq_NLS()
            obj.NMS.mod_Amp_NLS()
            %obj.NMS.plot('NLS')            
            obj.NMS.Freq_BSA()
            obj.NMS.Ampl_BSA()
            %obj.NMS.plot('BSA')
            fprintf('%s THD_NLS = %1.3e, R-sq_BSA = %1.3f, funevals = %d\n',obj.NMS.Name,obj.NMS.THD_NLS,obj.NMS.R_Sq_BSA,obj.NMS.funevals)
        end
        
    end
    
end

