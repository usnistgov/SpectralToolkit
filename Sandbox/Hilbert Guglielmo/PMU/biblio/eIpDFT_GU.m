function [f, A, ph, delta, eps] = eIpDFT_GU(X,N,kmax,fs,window,P,interpolator,f_old,A_old,ph_old)

if f_old == 0 || P ==0
    [f, A, ph, delta, eps] = IpDFT_GU(X,N,kmax,fs,window,interpolator);
else
    f=f_old;
    A=A_old;
    ph=ph_old;
end

% for i=1:P
%     Xn = spec_imag(-f,A,-ph,N,fs,window);
%     Xp = X-Xn;
%     [f, A, ph, delta, eps] = IpDFT(Xp,N,kmax,fs,window,interpolator);
% end
end

