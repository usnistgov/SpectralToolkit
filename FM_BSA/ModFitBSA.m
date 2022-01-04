function [Synx,Freqs,ROCOFs, iterations] = ModFitBSA(Fin,Fm,Km,Samples,dT,MagCorr,DelayCorr)
%MODFITBSA Summary of this function goes here
%Fitter for modulated signals with modulation index greater than 0.2
%   Fits to a model with only 2 sidebands
% inputs:
%   Fin: Fundamental frequency
%   Fm: Modulation frequency
%   Km: Modulation Index
%   Samples: Rows of samples, columns of phases
%   dt = Sampling period (1/Sample Rate)
%
% Citation:
%Kyriazis G. A., "Estimating parameters of complex modulated signals from
%prior information about their arbitrary waveform components," IEEE Trans.
%Instrum. Meas., v. 62, no. 6, pp. 1681-1686, June 2013.
%
% Citation:
% Kyriazis G. A., "Cartesian method to improve the results and
% save computation time in Bayesian signal analysis,"? in Advanced
% Mathematical and Computational Tools in Metrology and Testing X (AMCTM
% X), Series on Advances in Mathematics for Applied Sciences, vol. 86, F.
% Pavese; W. Bremser; A.G. Chunovkina; N. Fischer; A.B. Forbes (eds.),
% World Scientific, 2015, pp. 229-240.

% With additional grid search by Allen Goldstein (NIST) (unpublished)


% for debugging and visualization
verbose = false;
debug = true;
fig = 1;
res = 30;
zp = zeros(res,res);

% Basic configuration
[nSamples,nPhases] = size(Samples);
nHarm = 3;          % the number of harmonics to analyse
iFun = 2*nHarm+1;   % number of orthogonal vectors in the hypothesis


% iterMax_BSA = 5000;
% epsilon_BSA = 1e-8;
% rho_BSA = [0.85, 0.85, 0.85];
grid = 20;      % row length for the initial grid search (increase for higher delta-frequency)

% A-priori knowledge about the modulating signal
Delta_Freq = Km.*Fm;
mod_Freq = 2*pi*Fm*dT;       % modulation angular frequency normalized for sample rate
mod_Ampl = Km;               % modulation amplitude
mod_Phase = zeros(1,nPhases); % modulation phase assumed 0, BSA will find this

% Function scope variables
V_norm = [];        % normalized eiganvector of the model
invD_norm = [];     % normalized inverted eiganvalues of the model
sigma = 0;
hi = zeros(nHarm,1);

% results
Modulo_BSA = zeros(2,nPhases);
Phi_BSA = zeros(2,nPhases);
Phim_BSA = zeros(1,nPhases);

% BSA of the modulated signal using results of the NLS of the modulated signal
%
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
%thresh = -0.5*ePoint;

% The threshold is scaled by the number of samples.  The above was sampled
% at 4800 samples per second.
thresh = ePoint * (-.5/(4800*dT));

if verbose
    fprintf('Threshhold = %f, Fm = %f, Km = %f',thresh(1),Fm(1),Km(1))
end

for phase = 1:nPhases
    
    % for each phase, set up the grid search parameters
    startpt(1) = 2*pi*Fin(phase)*dT;   % carrier angular frequency normalized for samplerate
    OMEGA2 = linspace(-pi,pi,grid);
    %OMEGA3 = linspace(0,4*2*pi*Delta_Freq(phase)*dT,grid);
    OMEGA3 = linspace(0,2*2*pi*Delta_Freq(phase)*dT,grid);
    
    z = zeros(grid,grid);
    funEvals = zeros(1,nPhases);               % count the function evaluations
    
    if debug
        figure(fig), fig=fig+1;
        dF = 2*pi*Delta_Freq(phase)*dT;
        fcontour3([startpt(1),startpt(1);-pi,pi;0,2*dF],res)
        hold on
    end
     
    % perform the grid search
    zWorst = f([startpt(1),OMEGA2(1),OMEGA3(1)]);
    zBest = zWorst;
    idxBest = [1,1];
    idxWorst = [1,1];
    for m = 1:grid        % Delta-Freq in columns
        %zWorst = f([startpt(1),OMEGA2(1),OMEGA3(m)]);  % reset worst and best for each row
        %zBest = zWorst;
        for k = 1:grid    % Phi in rows
            z(k,m) =  f([startpt(1),OMEGA2(k),OMEGA3(m)]);
            if debug
                plot3(OMEGA2(k),OMEGA3(m),z(k,m),'.')
            end           
            if z(k,m) < zBest
                zBest = z(k,m);
                idxBest = [k,m];
            end
            if z(k,m) > zWorst, zWorst = z(k,m);end
            %if abs(zWorst-zBest) > abs(thresh(phase)),break,end
            if zBest-zWorst < thresh(phase),break,end
        end
        %if abs(zWorst-zBest) > abs(thresh(phase)),break,end
        if zBest-zWorst < thresh(phase),break,end        
    end

    % grid search found starting parameters
    startpt(2) = OMEGA2(idxBest(1));
    startpt(3) = OMEGA3(idxBest(2));
        
    % debugging contour plots ---------------------------------------------
    if debug
        hold off
    end
    % end debug ------------------------------------------------
    
    % Using fminsearch(problem)where problem is a structure.  see MATLAB "doc fminsearch"
    funEvals(phase)=0;               % restart the function counter
    opts = optimset(@fminsearch);    % configures NM options with default values
    opts.TolX = 1e-8;
    opts.TolFun = 1e8;
    PROBLEM = struct('objective',[],'X0',[],'options',opts,'solver','fminsearch');
    X0 = startpt;
    PROBLEM.X0 = X0;
    PROBLEM.objective = @(X0) f(X0);
    
    if debug
        PROBLEM.options.OutputFcn = @out;             % plots the contour map and points
        [endpt_BSA,fval,exitflag,outpoint] = fminsearch(PROBLEM);
        hold off
    else
        [endpt_BSA] = fminsearch(PROBLEM);
    end
    
    % endpt_BSA: [normalized modulation frequency, modulaton phase, delta frequency]
    
    % verbose status ------------------------------------------------------
    if verbose
        fprintf('Total iterations: %d\n',funEvals(phase))
        nfun = (iFun-1)/2;
        fList = '';
        for m = 1:nfun
            fItem = sprintf('%d: %f \n',m,endpt_BSA(m));
            fList = sprintf('%s %s',fList,fItem);
        end
        fprintf('BSA Frequencies (rad/unit):\n%s',fList);
    end
    % end verbose ---------------------------------------------------------
    
    
    % Calculates the Amplitude and Phase of the modulated signal after Freq_BSA() has been run
    % Generated the "best fit" modulated signal and the residual.
    %
    
    % preallocate vectors a and b
    nfun = (iFun-1)/2;
    a = zeros((nfun+1)/2,1);
    b = zeros((nfun+1)/2,1);
    
    B = V_norm*invD_norm;
    bi = B*hi;
    a(1) = bi(1);
    k = 2;
    while k <= nfun
        if k < (nfun+1)/2+1
            a(k) = bi(k);
        else
            b(k-(nfun-1)/2) = bi(k);
        end
        k = k+1;
    end
    
    cBSA = complex(a,b);        % complex
    Modulo_BSA(:,phase) = abs(cBSA);
    Phi_BSA(:,phase) = -angle(cBSA);
    Fcarr_BSA(phase) = endpt_BSA(1);
    Phim_BSA(phase) = endpt_BSA(2);
    dF_BSA(phase) = endpt_BSA(3);
    
    % debug: determine the best fit and residual then plot ----------------
    if debug
        % determine the best-fit signal
        wm = 2*pi*Fm(phase);
        n = double(0:nSamples-1)'*dT;   % time vector
        Result = Modulo_BSA(2,phase).*cos(Fcarr_BSA(phase)/dT.*n + dF_BSA(phase)/mod_Freq(phase) * sin(mod_Freq(phase)/dT.*n+Phim_BSA(phase))+Phi_BSA(2,phase));
        figure(fig);fig=fig+1;
        subplot(2,1,1)
        plot(n,Samples(:,phase),n,Result)
        subplot(2,1,2)
        plot(n,Samples(:,phase)-Result)
    end
    % end debug -----------------------------------------------------------
    
    
    
end % for phase = 1:nPhases


% Calculate the results: Synx, Freqs, ROCOFS, iterations
iterations = funEvals;

% angle at the center of the window
n = ceil(nSamples/2);
%wm = 2*pi*Fm;
wf = 2*pi*Fin;
Af = (Modulo_BSA(2,:)/sqrt(2)).*MagCorr;
Theta = Fcarr_BSA.*n + (dF_BSA/mod_Freq) .* sin(mod_Freq.*n + Phim_BSA) + Phi_BSA(2,:);
Theta = Theta + DelayCorr*1e-9.*wf;
Synx = (Af.*exp(-1i*Theta))';

% Frequency at the center of the window
Freqs = (Fcarr_BSA/(2*pi*dT) - (dF_BSA/mod_Freq) .* Fm .* cos(2*pi*Fm.*n+Phim_BSA))';
ROCOFs = ((dF_BSA/mod_Freq) .* Fm.^2.* sin(2*pi*Fm.*n+Phim_BSA)*2*pi)';

%% ========================================================================
% objective function
    function y = f(x)
        funEvals(phase) = funEvals(phase)+1;
        omega = x;  % [carrier freq, modulation phase, delta_Freq;
        t = double(0:nSamples-1)';   % time vector
        
        % The hypothesis
        GIJ = zeros(nSamples,3);
        GIJ(:,1) = 1;       %DC
        GIJ(:,2) = cos(omega(1)*t+omega(3)/mod_Freq(phase)*sin(mod_Freq(phase)*t+omega(2)));
        GIJ(:,3) = sin(omega(1)*t+omega(3)/mod_Freq(phase)*sin(mod_Freq(phase)*t+omega(2)));
        
        y = -prob(GIJ);
    end

    function stloge = prob(GIJ)
        iNo = double(nSamples);
        HIJ = ortho(GIJ);    % orthogonolize the objective function
        nFun = size(HIJ,2);
        hi = zeros(nFun,1);
        h2 = 0;
        for j=1:nFun
            h1 = sum(Samples(:,phase).*HIJ(:,j));
            hi(j)=h1;
        end
        h2 = sum(hi.^2);
        h2bar = h2/nFun;
        
        y2 = sum(Samples(:,phase).^2);
        y2bar = y2/iNo;
        
        qq = 1-h2/y2;
        if qq<=0
            qq = 1e-6;
        end
        
        stloge = log(qq)*(nFun-iNo)/2;
        dif = y2bar-nFun*h2bar/iNo;
        sigma = sqrt(abs(dif)*iNo/(iNo-nFun-2));
    end

    function HIJ = ortho(GIJ)
        % Orthogonalizes the hypothesis function
        M = GIJ'*GIJ;
        M = (M+M.')/2; % The matlab eig function must recognize the matrix as symetrical
        [V,D_vector] = eig(M,'vector');         % eiganvalues and eiganvectors
        SqrSumCol = sum(V.^2);
        norm = sqrt(SqrSumCol);
        %V_norm = V./norm;                       % matlab 2015 has issues with elementwise divide
        V_norm = bsxfun(@rdivide,V,norm);       % this is the workaround
        D_vec = sqrt(abs(D_vector));
        D_norm = diag(D_vec);
        invD_norm = inv(D_norm);
        A = GIJ*V_norm;
        HIJ = A*invD_norm;
    end
%% ========================================================================
% contour plots for debugging
    function fcontour3(omega,resolution)
        % plots a 3D contour map of the objective function f(x)
        % fcontour3(omega,fun, res)
        %
        %   omega:
        % an N x 2 matrix with N dimensions and the columns being the begin and
        % end points for each dimension.  All dimensions except 2 of
        % them must be the same begin and end point.
        %
        %   res:
        % The number f points to plot in the x and y axis.  The Z axis
        % will have resolution^2 points
        %
        %   fun:
        % Handle to the objective function
        % Create vectors of omega values to check.  Also x and y will be the
        % values to be plotted
        colOMEGA = size(omega,1);
        x = []; y = [];
        OMEGA = [];
        for ip = 1 : colOMEGA
            OMEGA = horzcat(OMEGA,linspace(omega(ip,1),omega(ip,2),resolution)');
            if omega(ip,1)~=omega(ip,2)
                if isempty(x)
                    x = OMEGA(:,ip);
                    if ip==1,xLbl='Carrier Freq (Hz)';else,xLbl='Phase (rad)';end
                else
                    y = OMEGA(:,ip);
                    if ip==2,yLbl='Phase (rad)';else,yLbl='Delta Freq (Hz/s)';end
                end
            end
        end
        count = 0;
        wb=waitbar(count,'Drawing Contour');
        for ip = 1:resolution
            for jp = 1:resolution
                zp(ip,jp) = f([OMEGA(ip,1),OMEGA(jp,2),OMEGA(ip,3)]);
                count = count+1;
                waitbar(count/resolution^2)
            end
        end
        close(wb)
        contour3(x,y,zp,30,'-k')
        xlabel(xLbl),ylabel(yLbl),zlabel('f eval')
    end


    function stop = out(x, optimValue, state)
        % fminsearch output function for ploting the trajectory during NM minimization
        stop = false;
        switch state
            case 'init'
                hold off
                f = 2*pi*Fin(phase)*dT;
                figure(fig);fig=fig+1;
                fcontour3([f,f;-pi,pi;0,2*dF],30)
                hold on;
            case 'iter'
                plot3((x(2)),abs(x(3)),optimValue.fval,'.');
                
                drawnow
        end
    end

end