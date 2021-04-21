% Class-based unit testing for Spectral Toolkit
%
%   run these tests with two command-line commands:
%   - testCase = testMFT;
%   - res = run(testCase);
%
classdef TestSpectral < matlab.unittest.TestCase
    %TESTSPECTRAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        localFunctions
    end
    
    methods (TestClassSetup)
        % any setup function calls should go here
    end
    
    methods (Test)
        % These functions will be automatically run at the command line command:
        % res = run(testCase);
        function regressionTests (testCase)
        end
    
end

