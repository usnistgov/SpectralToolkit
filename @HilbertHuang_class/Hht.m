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
    
    % frequency
    omega = gradient(unwrap(phaseAngle))/(IMFs.TimeInfo.Increment*2*pi);
    
    % limit the frequency range and resolution
    fRange = obj.hilOpts.FreqLimits;
    if isempty(fRange)
        fRange = [0,floor(0.5/(IMFs.TimeInfo.Increment))];
        obj.hilOpts.FreqLimits = fRange;
    end
    
    fRes = obj.hilOpts.FreqResolution;
    if isempty(fRes)
        fRes = (fRange(2)-fRange(1))/100;
        obj.hilOpts.FreqResolution = fRes;
    end        
    omega = fRes*round(omega/fRes);  % frequency resolution        
    omega = omega.*(omega>=fRange(1) & omega<=fRange(2)); % frequency range
    
    P(:,i) = {sig,energy,phaseAngle,omega};
end
obj.Hilbert = P;
end