% For a range of modulation frequencies and indices, plot the
% minimum value of the objective function.  Next use these
% minimum calues to create a 3D plot versus the Fm and Km.
%
% the values of this plot, will then be used to find a function
% which can predict good grid search stopping threshold values.

% Analytic Time Series class.  This is the object that created the signal
% and contains properties and methods to manipulate it
TS = AnalyticTS_class(...
                         'Name','Grid Search Threshold Demo',...
                         'F0', 50,...
                         'SampleRate', 4800,...
                         'Duration', 2.0...
                         );
                     
[~, Fin, ~, ~, ~, ~, Fa, Ka] = TS.getParamIndex();  % indexes into SignalParams
                     

FmStart = 0.5;
FmEnd = 5;
FmIncr = 0.5;

KmStart = 0.5;
KmEnd = 5;
KmIncr = 0.5;

grid = 20;
res = 30;

dT = 1/TS.SampleRate;       % Sampling Period
wCarr = 2*pi*TS.F0*dT;      % Carrier angular frequency

% Set up a table for the results
nRows = (round((FmEnd-FmStart)/FmIncr)+1)*(round((KmEnd-KmStart)/KmIncr)+1);
varNames = {'Km','Fm','zMin'};
varTypes = {'double','double','double'};
T = table('Size',[nRows,length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

p = 1;   %counter for the table
Fig = figure(1);

% loop through the Km values for each Fm
for Fm = FmStart:FmIncr:FmEnd
    for Km = KmStart:KmIncr:KmEnd
        TS.SignalParams(Fa,:) =Fm;
        TS.SignalParams(Ka,:) = Km;
        dF = 2*pi*Km*Fm*dT;
        AnalysisCycles =  ceil(TS.F0/TS.SignalParams(Fa,1));  % Analyse one modulation cycle of data
        Window = TS.getWindow(0,AnalysisCycles,'even');  %Get a Window of data
        
        % instantiate the FM_BSA class
        FM = FM_BSA_Class(...
            TS.SignalParams(Fin,:),...
            TS.SignalParams(Fa,:),...
            TS.SignalParams(Ka,:),...
            dT,...
            real(Window.Data)...
            );
        
        FM.phase = 1;
        startPt1 = 2*pi*TS.SignalParams(Fin,1)*dT;
        OMEGA2 = linspace(-pi,pi,grid);
        OMEGA3 = linspace(0,2*dF,grid);
        
        OMEGA = [startPt1,startPt1;-pi,pi;0,2*dF];
        FM.fcontour3(OMEGA,res,@FM.objFun)
        hold on
        
        z = zeros(grid,grid);
        for k = 1:grid
            for j = 1:grid
                z(j,k) = FM.objFun([startPt1,OMEGA2(j),OMEGA3(k)]);
                plot3(OMEGA2(j),OMEGA3(k),z(j,k),'*')
                
            end
        end
        
        hold off
        
        % Store the parameters and the minimum function value
        T(p,:) = {Km,Fm,min(min(z))};
        p = p+1;
                              
    end
end
