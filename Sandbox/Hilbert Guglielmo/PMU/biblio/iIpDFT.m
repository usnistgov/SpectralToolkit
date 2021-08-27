function [f, A, ph, delta, eps, fi, Ai, phi, deltai, epsi] = iIpDFT(s,N,kmax,fs,Q,P,interpolator)
k = 0:1:2*kmax+1;

w = hann(N,'periodic')';
K = 2/(0.5*N);

X = K*fftshift(fft(s.*w));
X = X(1+N/2+(k));

[f, A, ph, delta, eps] = eIpDFT(X,N,kmax,fs,P,0,0,0,k,interpolator);
fi = 0;
Ai = 0;
phi= 0;
deltai = 0;
epsi = 0;

for i=1:Q
    X0 = spec_imag(f,A,ph,N,fs,k)+spec_imag(-f,A,-ph,N,fs,k);
    Xi = X-X0;
    if sum(abs(Xi).^2)/sum(abs(X).^2) > 3.3e-3
        [fi, Ai, phi, deltai, epsi] = eIpDFT(Xi,N,-1,fs,P,fi,Ai,phi,k,interpolator);
        Xi = spec_imag(fi,Ai,phi,N,fs,k)+spec_imag(-fi,Ai,-phi,N,fs,k);
        X0 = X - Xi;
        [f, A, ph, delta, eps] = eIpDFT(X0,N,kmax,fs,P,f,A,ph,k,interpolator);
    else
        break
    end
end
end

