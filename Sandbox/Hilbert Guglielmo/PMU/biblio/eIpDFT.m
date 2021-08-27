function [f, A, ph, delta, eps] = eIpDFT(X,N,kmax,fs,n,f_old,A_old,ph_old,k,interpolator,window)

if f_old == 0
    [f, A, ph, delta, eps] = IpDFT(X,N,kmax,fs,interpolator,window);
else
    f=f_old;
    A=A_old;
    ph=ph_old;
end

for i=1:n
    Xn = spec_imag(-f,A,-ph,N,fs,k);
    Xp = X-Xn;
    [f, A, ph, delta, eps] = IpDFT(Xp,N,kmax,fs,interpolator,window);
end
end

