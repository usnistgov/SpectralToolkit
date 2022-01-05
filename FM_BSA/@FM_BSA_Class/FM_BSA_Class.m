% LICENSE:
%NIST-developed software is provided by NIST as a public service. You may
%use, copy, and distribute copies of the software in any medium, provided
%that you keep intact this entire notice. You may improve, modify, and
%create derivative works of the software or any portion of the software,
%and you may copy and distribute such modifications or works. Modified
%works should carry a notice stating that you changed the software and
%should note the date and nature of any such change. Please explicitly
%acknowledge the National Institute of Standards and Technology as the
%source of the software.

%NIST-developed software is expressly provided "AS IS." NIST MAKES NO
%WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION
%OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF
%MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND
%DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF
%THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL
%BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING
%THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED
%TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

%You are solely responsible for determining the appropriateness of using
%and distributing the software and you assume all risks associated with its
%use, including but not limited to the risks and costs of program errors,
%compliance with applicable laws, damage to or loss of data, programs or
%equipment, and the unavailability or interruption of operation. This
%software is not intended to be used in any situation where a failure could
%cause risk of injury or damage to property. The software developed by NIST
%employees is not subject to copyright protection within the United States.


classdef FM_BSA_Class
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
    
    properties (Access = public)
        verbose = false;    % If true, computed values will be displayed in the console 
        debug = false;      % If true, plots showing contours and intermediate values will be displayed
        fig = 1;
                
        % Modulated signal starting parameters
        Fcarr               % carrier frequency (required input)
        Fm                  % modulation frequncy (required input)
        Km                  % index of modulation (required input)
        dT                  % sampling period (required input)
        Samples             % 1D or 2D array of samples in rows and phases in columns
        phase               % which phase are we working on?
        nPhases
        
        % correction factors (used for hardware sampling systems)
        MagCorr             % magnitude correction factor
        DlyCorr           % delay correction factor
                
    end
    
    properties (Access = private)
        contourRes = 30;    % Resolution for the objective function contour plots
        grid = 20;          % length and width of the startpoint grid search
        nSamples  
        funEvals = 1;
        
    end
    
    %% =========================================================================
    % Constructor
    methods
        function obj = FM_BSA_Class(FCarr, Fm, Km, dT, Samples, varargin)
            % Constructor
            % The constructor accepts name,value pair arguments.  If the
            % argument is not included in the constructor call the default
            % value will ne used.  the arguments and their default values
            % are shown here:
            %
            % Example: FM_BSA_Class(Fcarr, Fm, Km, dT, Samples <,Name1,Value1,Name2,Value2,...NameN,ValueN,...>)
            %
            % Argument Name    , type   ,   default value           % comment
            % 'Fcarr'          ' [double]   '  n/a (required input) % carrier frequency or [frequencies]
            % 'Fm'             ' [double]   '  n/a (required input) % modulation frequency or [frequencies]
            % 'Km'             ' [double]   '  n/a (required input) % modulation index or [indices]
            % 'dT'             ' double     '  n/a (required input) % sampling period
            % 'Samples'        ' [double]   '  n/a (required input) % 1D or 2D array of samples in rows and phases in columns
            % 'MagCorr'        ' [double]   ' [1, 1, 1,...]         % Magnitude correction factors (one per phase)
            % 'DlyCorr'        ' [double]   ' [0, 0, 0...]          % Delay correction factors (one per phase)
            
            % 'verbose'        , Logical,   'false'                 % If true, computed values will be displayed in the console
            % 'debug'          , Logical,   'false'                 % If true, plots showing intermediate values will be displayed
            
            p = inputParser;
            p.addRequired('Fcarr')
            p.addRequired('Fm')
            p.addRequired('Km')
            p.addRequired('dT')
            p.addRequired('Samples')
            [obj.nSamples,obj.nPhases]=size(Samples);
            p.addParameter('MagCorr',ones(1,obj.nPhases))
            p.addParameter('DlyCorr',zeros(1,obj.nPhases))
            p.addParameter('verbose',false)            
            p.addParameter('debug',false)
                        
            parse(p, FCarr, Fm, Km, dT, Samples, varargin{:})
            obj.Fcarr = p.Results.Fcarr;
            obj.Fm = p.Results.Fm;
            obj.Km = p.Results.Km;
            obj.dT = p.Results.dT;
            obj.Samples = p.Results.Samples;
            obj.MagCorr = p.Results.MagCorr;
            obj.DlyCorr = p.Results.DlyCorr;
            obj.verbose = p.Results.verbose;
            obj.debug = p.Results.debug;
                        
        end
    end
    
    %% =========================================================================
    % Public Methods in external files
    methods (Access = public)
        [startPoint] = GridSearch(obj);
        y = objFun(obj, x); 
        [endpt_BSA,argout] = BSA_Est(obj,startPt);
        [estParams] = Param_Est(obj,endpt_BSA);
        [Synx,Freq,ROCOF] = Synx_Calc(obj,estParams)
    end

    %% =========================================================================
    % Static Methods in external files
    methods(Static)
        fcontour3(omega,resolution,fun);
    end
    
end

