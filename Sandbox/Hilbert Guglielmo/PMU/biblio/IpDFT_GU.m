function [f, A, ph, delta, eps] = IpDFT_GU(X,N,kmax,fs,window,interpolator)

df=fs/N;
if kmax==-1
    [~,kmax]=max(abs(X));
    kmax=kmax-1-1;
end

% if kmax==0
%     eps = -1;
%     kmax = kmax+1;
% elseif kmax==numel(X)-1
%     eps = 1;
%     kmax = kmax-1;
% else
if abs(X(kmax+1+1+1)) >= abs(X(kmax+1-1+1))
    eps = 1;
elseif abs(X(kmax+1+1+1)) <= abs(X(kmax+1-1+1))
    eps = -1;
end

indexkmax = kmax+1+1;
Xmax = X(indexkmax);

Xepsp = X(indexkmax+eps);
Xepsm = X(indexkmax-eps);

switch window
    case 'rect'
        switch interpolator
            case '2p'
                if abs(Xmax)+abs(Xepsp) == 0
                    delta = 0;
                else
                    delta = eps*(abs(Xepsp))/(abs(Xmax)+abs(Xepsp));
                end
            case '3p'
                if 2*abs(Xmax)-abs(Xepsp)+abs(Xepsm) == 0
                    delta = 0;
                else
                    delta = eps*(abs(Xepsp)+abs(Xepsm))/(2*abs(Xmax)+abs(Xepsp)-abs(Xepsm));
                end
        end
        if delta == 0
            A = abs(Xmax);
        else
            A = abs(Xmax)*abs((delta*pi)/sin(delta*pi));
        end
        
    case 'hann'
        switch interpolator
            case '2p'
                if abs(Xmax)+abs(Xepsp) == 0
                    delta = 0;
                else
                    delta = eps*(2*abs(Xepsp)-abs(Xmax))/(abs(Xmax)+abs(Xepsp));
                end
            case '3p'
                if abs(Xepsm)+2*abs(Xmax)+abs(Xepsp) == 0
                    delta = 0;
                else
                    delta = 2*eps*(abs(Xepsp)-abs(Xepsm))/(abs(Xepsm)+2*abs(Xmax)+abs(Xepsp));
                end
        end
        if delta == 0
            A = abs(Xmax)*abs(delta^2-1);
        else
            A = abs(Xmax)*abs((delta*pi)/sin(delta*pi))*abs(delta^2-1);
        end
        
    case 'cos'
        switch interpolator
            case '2p'
                if abs(Xmax)+abs(Xepsp) == 0
                    delta = 0;
                else
                    delta = eps*(1.5*abs(Xepsp)-0.5*abs(Xmax))/(abs(Xmax)+abs(Xepsp));
                end
            case '3p'
                if 2*abs(Xmax)+abs(Xepsp)+abs(Xepsm) == 0
                    delta = 0;
                else
                    delta = 1.5*eps*(abs(Xepsp)-abs(Xepsm))/(2*abs(Xmax)+abs(Xepsp)+abs(Xepsm));
                end
        end
        if abs(delta^2-0.25) == 0
            A = abs(Xmax)*4/pi; 
        else
            A = abs(Xmax)*4*abs(delta^2-0.25)/abs(cos(delta*pi));
        end
end
f = (kmax + delta)*df;
ph = wrapToPi(angle(Xmax)-pi*delta);
end