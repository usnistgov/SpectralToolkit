function obj = Hht(obj)

IMFs = obj.IMFs;
N = size(IMFs.Data,2);
P = cell(4,N);

for i = 1:N
    
    % Window the data if requested
    Ns = length(IMFs.Data(:,i));
    switch (obj.Window)
        case {'none';'None'}
            win = rectwin(Ns);
        case {'Hann';'Hanning'}
            win = hann(Ns);
        case {'Hamming','Hamm','hamming','hamm'}
            win = hamming(Ns);
        case {'BH'}
            win = blackmanharris(Ns);
        case {'Triangular';'triangular';'Triang';'triang';'tri'}
            win = triang(Ns);
        case {'flattop';'Flattop';'FlatTop'}
            win = flattopwin(Ns);
        otherwise
            warning('Window Type %s not recognized',obj.Window)
    end
    
    Data = IMFs.Data(:,i).*win;  % apply the window
    
    sig = hilbert(Data);
    energy = abs(sig).^2;
    phaseAngle = angle(sig);
    omega = gradient(unwrap(phaseAngle))/IMFs.TimeInfo.Increment;
    
    P(:,i) = {sig,energy,phaseAngle,omega};
end
obj.Hilbert = P;
end