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
            testModIndexRange(obj)
        end
        

    end
    
    %% ====================================================================
    methods (Access = private)
        % First test will check the default HJS settings
        function obj = testDefaultHJS(obj)
           obj.HJS = FM_HJS_Class(); 
           obj.HJS.configure();         
           obj.HJS.mod_Freq_NLS()
           obj.HJS.mod_Amp_NLS() 
           obj.HJS.plot('NLS')
           obj.HJS.Freq_BSA()
           obj.HJS.Ampl_BSA()
           obj.HJS.plot('BSA')
           fprintf('%s THD_NLS = %d, R-sq_BSA = %1.3f\n',obj.HJS.Name,obj.HJS.THD_NLS,obj.HJS.R_Sq_BSA)           
        end
        
        function obj = testModIndexRange(obj)
            obj.HJS.verbose = false;
            Kincr = 0.5;
            Kstart = 0.1;
            Kend = 10.1;
            [~,~,~,~,~,~,~,Ka] = obj.HJS.getParamIndex();
            k = Kstart;
            while k <= Kend
                %obj.HJS.Name = sprintf('50F0_5Fa0_%dKa%i',floor(k),floor((k-floor(k))*10));
                obj.HJS.Name = sprintf('50F0_5Fa0_%dKa%i',fix(k),floor(abs(fix(k)-k)*10));
                obj.HJS.fig = 1;
                obj.HJS.data = [];      % clear the modulating signal data from the preceeding test
                obj.HJS.data1 = [];      % clear the modulated signal data from the preceeding test
                obj.HJS.SignalParams(Ka) = k;
                obj.HJS.configure();
                obj.HJS.mod_Freq_NLS()
                obj.HJS.mod_Amp_NLS()
                obj.HJS.plot('NLS')
                obj.HJS.Freq_BSA()
                obj.HJS.Ampl_BSA()
                obj.HJS.plot('BSA')  
                fprintf('%s THD_NLS = %d, R-sq_BSA = %1.3f\n',obj.HJS.Name,obj.HJS.THD_NLS,obj.HJS.R_Sq_BSA)
                k = k + Kincr;
                %pause
            end
        end
        
        
        
    end
    
end

