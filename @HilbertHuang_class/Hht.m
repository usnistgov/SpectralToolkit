function P = Hht(obj)

IMFs = obj.IMFs;
N = size(IMFs,2);
P = cell(4,N);

for i = 1:N
    sig = hilbert(IMFs(:,i));
    energy = abs(sig).^2;
    phaseAngle = angle(sig);
    omega = gradient(unwrap(phaseAngle));
    
    % ----------- plots for debugging ------------------------------------
    figure()
    title(sprintf('IMF %d', i))
    subplot(4,1,1)
    plot(real(sig),'-b')
    pause
    hold on
    plot(imag(sig),'-r')
    hold off
    ylabel('analytic')
    pause
    subplot(4,1,2)
    plot(energy)
    ylabel('energy')
    pause
    subplot(4,1,3)
    plot(phaseAngle)
    ylabel('angle')
    pause
    subplot(4,1,4)
    plot(omega)
    ylabel('freq')
    pause
    %---------------------------------------------------------------------
    
    P(i,:) = {sig,energy,phaseAngle,omega};
end
obj.IMF = P;
end