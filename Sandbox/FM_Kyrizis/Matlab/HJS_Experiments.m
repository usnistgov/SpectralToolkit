classdef HJS_Experiments < matlab.unittest.TestCase
    % Class-based unit testing for Spectral Toolkit
%
%   run these tests with two command-line commands:
%   - testCase = HJS_Experiments;
%   - res = testCase.run OR res = run(testCase);
%
    
    properties
        HJS     % a single instance of a FM_HJS class

    end
    
    %% Constructor ========================================================
    methods 
        function obj = HJS_Experiments(obj)
            obj.HJS = FM_HJS_Class();
        end
    end
    
   %% Test methods ======================================================== 
    methods (Test)
        % These functions will be automatically run at the command line command:
        % res = testCase.run OR res = run(testCase);
        function Tests(obj)
            %testDefaultHJS(obj);
            %testModIndexRange(obj)
            troubleshootLocalMax(obj)
        end
        

    end
    
    %% ====================================================================
    methods (Access = private)
        % First test will check the default HJS settings
        function obj = testDefaultHJS(obj)
           obj.HJS = FM_HJS_Class(); 
           obj.runOne()
        end
        
        function obj = testModIndexRange(obj)
            % ranges across modulation indexes from 0.1 to 10 in 0.5 increments.  With noise
            % As of this update, the Freq_BSA is not finding correct answers for Ka = 2.6,3.1,3.6,6.1,7.1,8.1 
           obj.HJS.modNoise = false;    % remove amplitude noise
           obj.HJS.phaseNoise = false;   % remove phase noise

            obj.HJS.verbose = false;
            Kincr = 0.5;
            Kstart = 0.1;
            Kend = 10;
            [~,Fin,~,~,~,~,Fa,Ka] = obj.HJS.getParamIndex();
            obj.HJS.SignalParams(Fin) = 50.0;
            obj.HJS.SignalParams(Fa) = 5.0;
            k = Kstart;
            while k <= Kend
                obj.HJS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
                %obj.HJS.fig = 1;
                obj.HJS.SignalParams(Ka) = k;
                obj.runOne()  
                obj.HJS.plot('BSA')
                
                k = k + Kincr;
                %pause
            end
        end
        
        function obj = troubleshootLocalMax(obj)
            % 50f0_5Fa0_3Ka1 with no additive noise
           %obj.HJS.fig = 1;             % reset the figure number 
           obj.HJS.hookeContour = {}; % clear the contour from last run
           obj.HJS.verbose = true;
           obj.HJS.debug = true;      % set the debug mode
           obj.HJS.AmplNoiseParams.NoiseGaussSD = 0;    % remove amplitude noise
           obj.HJS.PhaseNoiseParams.NoiseGaussSD = 0;   % remove phase noise
           
           [~,~,~,~,~,~,~,Ka] = obj.HJS.getParamIndex();
           %k = 4.1;
           k = 5.0;
           obj.HJS.SignalParams(Ka) = k;
           obj.HJS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),round(abs(fix(k)-k)*10));
           obj.runOne()
           obj.HJS.plot('BSA')           
           obj.HJS.plot('hookeContour');
          

        end
        
        
        function runOne(obj)
            % runs one simulation of HJS
            obj.HJS.data = [];      % clear the modulating signal data from the preceeding test
            obj.HJS.data1 = [];      % clear the modulated signal data from the preceeding test
            obj.HJS.configure();
            obj.HJS.mod_Freq_NLS()
            obj.HJS.mod_Amp_NLS()
            %obj.HJS.plot('NLS')
            obj.HJS.Freq_BSA()
            obj.HJS.Ampl_BSA()
            %obj.HJS.plot('BSA')
            fprintf('%s THD_NLS = %1.3e, R-sq_BSA = %1.3f, funevals = %d\n',obj.HJS.Name,obj.HJS.THD_NLS,obj.HJS.R_Sq_BSA,obj.HJS.funevals)
        end
        
    end
    
end

