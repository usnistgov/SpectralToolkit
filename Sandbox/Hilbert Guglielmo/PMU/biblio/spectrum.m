function [DFT,DFTp,DFTn,DTFT,DTFTp,DTFTn,k,kc] = spectrum(f,A,ph,N,Nc,fs,window)

k =     (0:1:N); kc = k;
% kc  =   (-N/2:0.01:N/2-1/100);

DFT =   zeros(1,numel(k));
DFTn =  zeros(1,numel(k));
DFTp =	zeros(1,numel(k));
DTFTn = zeros(1,numel(kc));
DTFTp = zeros(1,numel(kc));
DTFT =  zeros(1,numel(kc));

switch window
    case 'rect'
        for i = 1:length(f)
            DFTn =  DFTn + A(i)*DK_k(k+f(i)/fs*N,N,-ph(i));
            DFTp =  DFTp + A(i)*DK_k(k-f(i)/fs*N,N, ph(i));
            DFT =   DFTn + DFTp;
            DTFTn = DTFTn + A(i)*DK_k(kc+f(i)/fs*N,N,-ph(i));
            DTFTp = DTFTp + A(i)*DK_k(kc-f(i)/fs*N,N, ph(i));
            DTFT  = DTFTn + DTFTp;
        end
    case 'hann'
        for i = 1:length(f)
            DFTn =  DFTn + 1/0.5*A(i)*HnW_k(k+f(i)/fs*N,N,-ph(i));
            DFTp =  DFTp + 1/0.5*A(i)*HnW_k(k-f(i)/fs*N,N, ph(i));
            DFT =   DFTn + DFTp;
            DTFTn = DTFTn + 1/0.5*A(i)*HnW_k(kc+f(i)/fs*N,N,-ph(i));
            DTFTp = DTFTp + 1/0.5*A(i)*HnW_k(kc-f(i)/fs*N,N, ph(i));
            DTFT  = DTFTn + DTFTp;
        end
    case 'hamm'
        for i = 1:length(f)
            DFTn =  DFTn + 1/0.54*A(i)*HmW_k(k+f(i)/fs*N,N,-ph(i));
            DFTp =  DFTp + 1/0.54*A(i)*HmW_k(k-f(i)/fs*N,N, ph(i));
            DFT =   DFTn + DFTp;
            DTFTn = DTFTn + 1/0.54*A(i)*HmW_k(kc+f(i)/fs*N,N,-ph(i));
            DTFTp = DTFTp + 1/0.54*A(i)*HmW_k(kc-f(i)/fs*N,N, ph(i));
            DTFT  = DTFTn + DTFTp;
        end
    case 'cos'
        DFTw = cosW_k(k,N,0);
        C = max(abs(DFTw))/N;
        for i = 1:length(f)
            DFTn =  DFTn + 1/(C*N)*A(i)*cosW_k(k+f(i)/fs*N,N,-ph(i));
            DFTp =  DFTp + 1/(C*N)*A(i)*cosW_k(k-f(i)/fs*N,N, ph(i));
            DFT =   DFTn + DFTp;
            DTFTn = DTFTn + 1/(C*N)*A(i)*cosW_k(kc+f(i)/fs*N,N,-ph(i));
            DTFTp = DTFTp + 1/(C*N)*A(i)*cosW_k(kc-f(i)/fs*N,N, ph(i));
            DTFT  = DTFTn + DTFTp;
        end
end

end