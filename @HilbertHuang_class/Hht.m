function obj = Hht(obj)

IMFs = obj.IMFs;
N = size(IMFs.Data,2);
P = cell(4,N);

for i = 1:N
    sig = hilbert(IMFs.Data(:,i));
    energy = abs(sig).^2;
    phaseAngle = angle(sig);
    omega = gradient(unwrap(phaseAngle))/IMFs.TimeInfo.Increment;
    
    P(:,i) = {sig,energy,phaseAngle,omega};
end
obj.Hilbert = P;
end