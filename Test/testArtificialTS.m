% Class-based unit testing for ArtificialTS.m
%
%   run these tests with two command-line commands:
%   - testCase = testArtificialTS;
%   - res = run(testCase);
%
classdef testArtificialTS < matlab.unittest.TestCase
    
    properties
        TestFigure
        TS
    end
    
    methods (TestClassSetup)
        function constructTS (testcase)
            testCase.TS = ArtificialTS;
        end
    end
            
    
    
    methods (TestMethodSetup)
        function createFigure(testCase)
            testCase.TestFigure = figure;
            testCase.TS = ArtificialTS;
        end
    end
    
    methods (TestMethodTeardown)
        function closeFigure(testCase)
            close(testCase.TestFigure)
        end
    end
    
    
    methods (Test)
        
        function testSetMethods (testCase)
            testCase.TS = ArtificialTS;
            testCase.TS.Name = 'testTS';
        end
        
        function testArtTS (testCase)
            Name = 'testTS';
            Description = 'Artificial Time Series Test Object';
            
            T0 = 0;
            Extent = 2;     % SECONDS
            nSamples = 960 * Extent;
            
            Freqs = [1, 2];
            Amps = [1, 1];
            Phases = [0 0];
            
            NoiseUniformLow = .05;
            NoiseUniformHi = .1;
            NoiseGaussMean = 0;
            NoiseGaussSD = 0;
            
            testCase.TS = ArtificialTS(...
                Name, ...
                Description, ...
                T0, ...
                Extent, ...
                nSamples, ...
                Freqs, ...
                Amps, ...
                Phases, ...
                NoiseUniformLow, ...
                NoiseUniformHi, ...
                NoiseGaussMean, ...
                NoiseGaussSD ...
                );
        
        plot(testCase.TS.time,testCase.TS.Ts);
        pause;
            
        end
        
        

                
                        
            
    end
end
        
        