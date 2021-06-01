function P = myHht(IMFs)

for i = 1:size(IMFs,2)
    sig = hilbert(IMFs(:,i));
    energy = abs(sig).^2;
    phaseAngle = angle(sig);
    omega = gradient(unwrap(phaseAngle));
    
    % ----------- plots
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
    %---------------

end