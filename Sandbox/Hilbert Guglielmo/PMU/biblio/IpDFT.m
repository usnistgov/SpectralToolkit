function [f, A, ph, delta, eps] = IpDFT(X,N,kmax,fs,interpolator,window)

df=fs/N;

if abs(X(kmax+1+1)) >= abs(X(kmax+1-1))
    eps = 1;
elseif abs(X(kmax+1+1)) <= abs(X(kmax+1-1))
    eps = -1;
end

indexkmax = kmax+1;
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
end

f = (kmax + delta)*df;
ph = wrapToPi(angle(Xmax)-pi*delta);
end