% Class-based unit testing for Spectral Toolkit
%
%   run these tests with two command-line commands:
%   - testCase = TestSpectral;
%   - res = testCase.run;
%
classdef TestSpectral < matlab.unittest.TestCase
    %TESTSPECTRAL Unit testing got the Spectral Toolbox
    %   Functions in the following methods all run in order automatically when the class is "run":
    %   TestClassSetup, TestMethodSetup, Test, TestMethodTeardown, TestClassTeardown
    
    properties
        TS      % Time Series class instance
        FS      % Fourier Series class instantance
        
        fig = 1;
    end
    
    methods (TestClassSetup)
        % any test class setup function calls should go here
         function TcSetup (testCase)
%             createArtTsFile (testCase)  % only used to Fourier-based so moved into regressionTests function
              testCase.TS = AnalyticTS_class();
         end
        
    end
    
       methods (TestMethodSetup)
        % any class setup function calls should go here
    end
    
     
    methods (Test)
        % These functions will be automatically run at the command line command:
        % res = run(testCase);
        function regressionTests (testCase)
            % All the listed functions will run.  Comment out any function you do not want to run
            
%             %--------------------------------------------------------------
%             %Tests of Fourier-based analysis methods
%             createArtTsFile (testCase)
%             testTimeSeriesClass (testCase)            
%             disp('Begin FourierSeries class testing')
%             testFourierSeriesClass (testCase)
%             testFourierPlots (testCase)
%             testFourierWindows (testCase)
%             % testLeakageFunctionPlots (testCase)   % not working yet
%             
            %--------------------------------------------------------------
            % Tests of Hilbert Huang based analysis
            testHhtFmAnalysis(testCase);
            
        end        
    end
    
    methods (TestMethodTeardown)
        % Any teardown of the test methods would go here
    end
    
    methods (TestClassTeardown)
        % any teardown of the test class would go here
        function TcTearDown(testCase)
            destroyArtTsFile (testCase)
        end
    end
    
    methods (Access = private)
        % These methods can only be called by other methods in the class itself
        function createArtTsFile (testCase)
            % Create a .csv file with columns of timeseries test data
            cTS = {};     % A temporary cell array to hold a collection of artTS objects
            
            % This first TS will have 60 Hz at 100 peak, 120 Hz at 50 peak and 75.25 at 20 peak
            % it will also have Uniform noise at +/- 1 (-40 dB from 60 Hz)
            Name = '60f0_75f25_120f0';
            Description = '60Hz@100 peak, 120Hz@50peak, 75.25Hz@20 peak plus 1 uniform noise';
            T0 = 0;
            Extent = 1;                 % seconds
            nSamples = uint32(9600 * Extent);   % Fs * extent
            
            Freqs = [60 120 75.25];
            Amps = [100 50 20];
            Phases = [0 0 60];
            
            NoiseUniformLow = 1;
            NoiseUniformHigh = 1;
            NoiseGaussMean = 0;
            NoiseGaussSD = 0;
            
            cTS{1} = ArtificialTS(...
                Name, ...
                Description, ...
                T0, ...
                Extent, ...
                nSamples, ...
                Freqs, ...
                Amps, ...
                Phases, ...
                NoiseUniformLow, ...
                NoiseUniformHigh, ...
                NoiseGaussMean, ...
                NoiseGaussSD ...
                );
            
            % Here, you can create more TS{n} with different parameters do not change
            % the Extent or nSamples or they wiill not be suitable to be written to the 
            % same .csv file
            % ...
            
            
            % After all the TS{n} have been created, loop through them to create a table
            % to be written to a .csv file.  This assumes the Extent and nSamples of all the time series 
            % are the same.
            for i = 1:numel(cTS)
                varName = cTS{i}.Name;
                data = cTS{i}.Ts;
                if i == 1
                    time = cTS{i}.time;
                    T = table(time',data');
                    T.Properties.VariableNames={'time',varName};
                else
                    t = table(data');
                    t.Properties.VariableNames={varName};
                    T = [T;t];
                end                
            end
            
            % Writes a .csv file to the working folder.  
            % This file will be deleted on Class destruction
            writetable(T,'ArtTS.csv')
        end
        
        %------------------------------------------------------------------
        function destroyArtTsFile (testCase)
            if exist('ArtTS.csv','file')
                delete('ArtTS.csv')
            end            
        end
        
        %------------------------------------------------------------------
        function testTimeSeriesClass (testCase)
            % the goal of this function is to exercize all the code in the TimeSeries Class
            disp('Begin TimeSeries class test')
            testCase.TS = TimeSeries_class('File','ArtTs.csv');    % instantiate the class and read in the file
            varNames = testCase.TS.dataTable.Properties.VariableNames;
            timeVals = testCase.TS.dataTable.(varNames{1});      % the first column of data should be the time vector
            
            % Test TS.(loadTs) function uging the second column of data
            testCase.TS.loadTs(timeVals,varNames{2});
            
            % test the plot function.  This just checks for errors, the user will not see a plot
            f = figure();
            plot(testCase.TS);       
            close(f)
            
            % Destroy the object
            clear TS
            
            % Set up to test the constructor             
            TimeInfo.Units = 'seconds';
            TimeInfo.Increment = mean(diff(timeVals));
            TimeInfo.Start = timeVals(1);
            
            %This constructor test is interactive and requires input from te user
            testCase.TS = TimeSeries_class('File','ArtTs.csv','TimeInfo',TimeInfo);
            f = figure();
            plot(testCase.TS)
            close(f)
            
            % This switching test is interactive and requires input from the user
            testCase.TS.switchTs()
            disp('TimeSeries class test complete')            
        end
        
        %------------------------------------------------------------------
            % the goal of the following set of functions is exercize all of the code in the
            % FourierSeries class
        
            function testFourierSeriesClass (testCase)                
                % set up a TimeSeies object using the first timeSeries in ArtTS
                testCase.TS = TimeSeries_class('File','ArtTs.csv');    % instantiate the class and read in the file
                varNames = testCase.TS.dataTable.Properties.VariableNames;
                timeVals = testCase.TS.dataTable.(varNames{1});      % the first column of data should be the time vector
                testCase.TS.loadTs(timeVals,varNames{2});
                
                % Instantiate a default FourierSeries object
                testCase.FS = FourierSeries_class('TimeSeries',testCase.TS.Ts);
            end            
            
            function testFourierPlots (testCase)
                % test all the plot options.  Looking for errors, The user will not see these plots
                f = figure();
                plot(testCase.FS,'Shift')
                %pause
                plot(testCase.FS,'DoubleSided');
                %pause
                plot(testCase.FS,'SingleSided','yScale','dB','xlim',[0, 200],'ylim',[-50, 70],'Residual')
                %pause
                close(f)                
            end
            
            function testFourierWindows (testCase)
                % test all the Windowing options (add more tests if you add more windows)
                testCase.FS.DataInfo.Window = 'Hamming';
                testCase.FS.calcFs(testCase.TS.Ts);
                plot(testCase.FS,'SingleSided','yScale','dB','xlim',[0, 200],'ylim',[-50, 70])
                %pause
                testCase.FS.DataInfo.Window = 'Hann';
                testCase.FS.calcFs(testCase.TS.Ts);
                plot(testCase.FS,'SingleSided','yScale','dB','xlim',[0, 200],'ylim',[-120, 70])
                %pause
                testCase.FS.DataInfo.Window = 'Blackman';
                testCase.FS.calcFs(testCase.TS.Ts);
                plot(testCase.FS,'SingleSided','yScale','dB','xlim',[0, 200],'ylim',[-120, 70])                
                %pause    
                close(f)
            end
            
            function testLeakageFunctionPlots (testCase)
                % the code in FourierSeries_class is not working yet.
            end
            
            %--------------------------------------------------------------
            % The goal of the following code is to ecersize all the code in the
            % AnalyticTS_class
            function testAnalyticTsClass
                
            end
            
            %--------------------------------------------------------------
            % code to test HHT analysis of windowed FM signal
            function testCase = testHhtFmAnalysis(testCase)
                
                % Get a window of FM data
                SignalParams = testCase.TS.SignalParams;
                [~,~,~,~,~,~,Fa,Ka] = testCase.TS.getParamIndex();
                SignalParams(Fa,:) = 2.0;
                SignalParams(Ka,:) = 2.5;
                % replace the default TS with a new one
                testCase.TS = AnalyticTS_class('SignalParams',SignalParams);
                Y = testCase.TS.getWindow(0,6);
                
                % Instatiate a HilbertHuang_class object.
                % Passing a timeseries object triggers the analysis
                increment = testCase.TS.Ts.TimeInfo.Increment;
                HHT = HilbertHuang_class('TimeSeries',timeseries(real(Y)));  
                
                % frequency is the derivitive of phase with respect to time
                increment = testCase.TS.Ts.TimeInfo.Increment;
                fVector = HHT.Hilbert{4,1}/(increment*2*pi);
                
                % calculate the actual frequency from the complex Y
                actAngle = angle(Y);
                actFreq = gradient(unwrap(actAngle))/(2*pi*increment);
                
                figure(testCase.fig); testCase.fig = testCase.fig + 1;
                subplot(2,1,1)
                plot(fVector,'r');hold on;plot(actFreq,'b');hold off
                subplot(2,1,2)
                plot(fVector-actFreq)
                
            end
            
        
    end
    
end

