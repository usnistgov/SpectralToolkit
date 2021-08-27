function [DFT,k] = spec_imag_GU(f,A,ph,N,fs,window,k)

DFT  =  zeros(1,numel(k));

switch window
    case 'rect'
        for i=1:length(f)
            DFT=DFT+A(i)*DK_k(k-f(i)/fs*N,N,ph(i));
        end
    case 'hann'
        for i=1:length(f)
            DFT=DFT+2*A(i)*HnW_k(k-f(i)/fs*N,N,ph(i));
        end
    case 'cos'
        B=sum(sin(pi/N*(0:N-1)));
        for i=1:length(f)
            DFT=DFT + N/B*A(i)*cosW_k(k-f(i)/fs*N,N,ph(i));
        end
end
end


% DFTp=DFTp + A(i)*DK_k(k-f(i)/fs*N,N, ph(i));
% DFTn=DFTn + A(i)*DK_k(k+f(i)/fs*N,N,-ph(i));