A scale of 0.5 should be sufficient to be sure thfunction [startpt] = GridSearch(obj)
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

Fin = obj.Fcarr;
Fm = obj.Fm;
Km = obj.Km;
dT = obj.dT;
Delta_Freq = Fm.*Km;

% Grid search threshold prediction
% The threshold for the grid search can be predicted using a function of Fm and Km
p00 = -10.88;
p10 = 0.3855;
p01 = 1.038;
p20 = -0.02687;
p11 = -0.0004427;
p02 = -0.1067;

thrLog = p00 + p10.*Fm + p01.*Km + p20.*Fm.^2 + p11.*Fm.*Km + p02.*Km.^2;
ePoint = exp(-thrLog);

% The threshold is scaled by the number of samples.  The above was sampled
% at 4800 samples per second.
thresh = ePoint * (-.5/(4800*obj.dT));

% display the threshold
if obj.verbose
    fprintf('Threshhold = %f, Fm = %f, Km = %f',thresh(1),Fm(1),Km(1))
end

% perform the grid search
for phase = 1:obj.nPhases
    
    % for each phase, set up the grid search parameters
    startpt(1) = 2*pi*Fin(phase)*dT;   % carrier angular frequency normalized for samplerate
    OMEGA2 = linspace(-pi,pi,obj.grid);
    OMEGA3 = linspace(0,2*2*pi*Delta_Freq(phase)*dT,obj.grid);
    
    z = zeros(obj.grid,obj.grid);
    
    if obj.debug
        figure(obj.fig), obj.fig=obj.fig+1;
        dF = 2*pi*Delta_Freq(phase)*dT;
        obj.fcontour3([startpt(1),startpt(1);-pi,pi;0,2*dF],obj.contourRes,@obj.objFun,phase)
        hold on
    end
     
    % perform the grid search
    zWorst = obj.objFun([startpt(1),OMEGA2(1),OMEGA3(1)],phase);
    zBest = zWorst;
    idxBest = [1,1];
    for m = 1:obj.grid        % Delta-Freq in columns
        for k = 1:obj.grid    % Phi in rows
            z(k,m) =  obj.objFun([startpt(1),OMEGA2(k),OMEGA3(m)],phase);
            if obj.debug
                plot3(OMEGA2(k),OMEGA3(m),z(k,m),'.')
            end           
            if z(k,m) < zBest
                zBest = z(k,m);
                idxBest = [k,m];
            end
            if z(k,m) > zWorst, zWorst = z(k,m);end
            if zBest-zWorst < thresh(phase),break,end
        end
        if zBest-zWorst < thresh(phase),break,end        
    end

    % grid search found starting parameters
    startpt(2) = OMEGA2(idxBest(1));
    startpt(3) = OMEGA3(idxBest(2));
    
    if obj.debug,hold off,end

end