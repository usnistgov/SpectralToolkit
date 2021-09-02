function obj = Emd(obj)

% empirical mode decomposition practice function
%
% Input:
%   x must be a timeseries object.  inly the first data column will be analuses
% Output:
%   obj.IMFs is a cell array of time series objects, 

x = real(obj.Ts_In.Data);

MaxNumIMFs = 10;
MaxEnergyRatio = 20;
MaxNumExtrema = 1;
SiftMaxIterations = 100;
SiftRelativeTolarance = 0.2;
Interpolation = 'spline';

if isrow(x)
    x = x(:);
end

rsig = x;       % residual signal
t = (1:length(x))';
N = length(x);
IMFs = zeros(N, MaxNumIMFs, 'double');

%extract IMFs
i = 0;
outerLoopExitFlag = 0;

while(i<MaxNumIMFs)
    % convergence checking
    %sprintf ('IMF %d\n',i);
    [peaksIdx, bottomsIdx] = localFindExtremaIdx(rsig);
        
    numResidExtrema = length(peaksIdx) + length(bottomsIdx);
    energyRatio = 10*log10(norm(x,2)/norm(rsig,2));
    
    if energyRatio > MaxEnergyRatio
        outerLoopExitFlag = 1;
        break
    end
    
    if numResidExtrema < MaxNumExtrema
        outerLoopExitFlag = 2;
        break
    end
    
    % SIFTING process initialization
    rsigL = rsig;
    rtol = ones(1,'double');
    k = 0;
    SiftStopCriterionHit = 'SiftMaxIteration';
    
    % Sifting process
    while (k<SiftMaxIterations)
       % check convergence 
       if(rtol<SiftRelativeTolarance)
          SiftStopCriterionHit = 'SiftMaxRelativeTolerance'; 
          break;          
       end
       
       % store previous residal
       rsigPrev = rsigL;
       
       % find the peaks
       [peaksIdx, bottomsIdx] = localFindExtremaIdx(rsigL);
       
       if((length(peaksIdx) + length(bottomsIdx))>0)
           %compute the upper and lower envelope 
           [uLoc, uVal, bLoc, bVal] = computeSupport(t, rsigL, peaksIdx, bottomsIdx);
           upperEnvelope(:,1) = interp1(uLoc, uVal, t, Interpolation);
           lowerEnvelope(:,1) = interp1(bLoc, bVal, t, Interpolation);
           
           % subtract mean envelope from residual
           mVal = (upperEnvelope + lowerEnvelope)/2;
       else
           mVal(1:N,1) = 0;
       end
       
%        %------PLOT the rsig and the peak indexes
%        figure()
%        plot(rsigL,'-k')
%        hold on
%        %pause
%        plot(peaksIdx,rsigL(peaksIdx),'og')
%        plot(bottomsIdx,rsigL(bottomsIdx),'db')
%        %pause
%        plot(upperEnvelope,'-g')
%        plot(lowerEnvelope,'-b')
%        %pause
%        plot (mVal,'-r')
%        hold off
%        %pause
%        %------
       
       rsigL = rsigL - mVal;
       
       % residual tolerance
       rtol = (norm(rsigPrev-rsigL,2)/norm(rsigPrev,2))^2;
       k = k + 1;
       
       
       
    end  % inner while loop
    
    % extract new IMF and subtract from residual
    IMFs(:,i+1) = rsigL;
    
    %----------- debug plot -------------------
    %figure()
    % plot(IMFs(:,i+1))
    %ylabel(sprintf('IMF %d',i+1));
    %pause
    %------------------------------------------
    rsig = rsig - IMFs(:,i+1);
    i = i +1;    
end  % outer while loop

obj.IMFs = timeseries(IMFs(:,1:i),obj.Ts_In.Time,'Name',strcat(obj.Ts_In.Name,' IMF'));
obj.IMFs = setuniformtime(obj.IMFs,'StartTime',obj.Ts_In.TimeInfo.Start,'EndTime',obj.Ts_In.TimeInfo.End);
obj.Residual = rsig;

end  % function

%--------------------------------------------------------------------------
function [locMax, locMin] = localFindExtremaIdx(x)
    nx = length(x);
    nmax = 0;
    nmin = 0;
    locMax = zeros(nx,1);
    locMin = zeros(nx,1);
    if nx >=3
        dright = x(2) - x(1);
        for k = 2:nx-1
            dleft = dright;
            dright = x(k+1) - x(k);
            if dleft > 0 && dright <= 0
                nmax = nmax + 1;
                locMax(nmax) = k;
            elseif dleft < 0 && dright >=0
                nmin = nmin +1;
                locMin(nmin) = k;              
            end            
        end
    end
    locMax = locMax(1:nmax,1);
    locMin = locMin(1:nmin,1);    
end

%--------------------------------------------------------------------------
function [uLoc, uVal, bLoc, bVal] = computeSupport(t, rsigL, pksIdx, btmsIdx)
% compute support for upper and lower envelope given input signal rsigL
N = length(t);
if(isempty(pksIdx))
    pksIdx = [1; N];
end

if(isempty(btmsIdx))
    btmsIdx = [1; N];
end

pksLoc = t(pksIdx);
btmsLoc = t(btmsIdx);

% compute envelop for wave method
% extended waves on the left
[lpksLoc, lpksVal, lbtmLoc, lbtmVal] = signalwavelet.internal.emd.emdWaveExtension(t(1), rsigL(1),...
    pksLoc(1), rsigL(pksIdx(1)),...
    btmsLoc(1), rsigL(btmsIdx(1)),...
    -1);

% extended waves on the right
[rpksLoc, rpksVal, rbtmLoc, rbtmVal] = signalwavelet.internal.emd.emdWaveExtension(t(end), rsigL(end),...
    pksLoc(end), rsigL(pksIdx(end)),...
    btmsLoc(end), rsigL(btmsIdx(end)),...
    1);

% append extended wave to extrema
uLoc = [lpksLoc;pksLoc;rpksLoc];
uVal = [lpksVal;rsigL(pksIdx);rpksVal];
bLoc = [lbtmLoc;btmsLoc;rbtmLoc];
bVal = [lbtmVal;rsigL(btmsIdx);rbtmVal];
end