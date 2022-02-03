T = readtable('GridSearchMinValues');
Km = T{:,'Km'};
Fm = T{:,'Fm'};
zMin = T{:,'zMin'};

KmStart = Km(1);KmEnd=Km(end);KmIncr=mean(diff(unique(Km)));
FmStart = Fm(1);FmEnd=Fm(end);FmIncr=mean(diff(unique(Fm)));

% plot the surface
%nX = length(unique(x)); nY=length(unique(y));
[X,Y] = meshgrid(KmStart:KmIncr:KmEnd,FmStart:FmIncr:FmEnd);
zMat=reshape (zMin,length(FmStart:FmIncr:FmEnd),length(KmStart:KmIncr:KmEnd));
figure(4);
surf(X,Y,zMat)
colormap('parula')
set(gca,'ZScale','log')
set(gca,'FontSize',12)
xlabel('Km')
ylabel('Fm')
zlabel('Minimum objective function value')
title('Minimum objective value contour')

% Figure 4 has a z scale which is log base 10, so we will transpose zmin into log base 10
%  NOTE: Log of negative numbers are COMPLEX numbers undefined in the real numbers but very much
% existing in the imaginary plane!
zLog = log10(zMin);

% in this case, the imaginary part is a constant pi/ln(10).  Here is the "proof"
%plot(imag(zLog)-(pi/log(10)))

SF = fit([Km,Fm],real(zLog),'poly22');  % fit is our cfit object
figure(5);
plot(SF,[Km,Fm],real(zLog))
colormap(flipud(parula));
xlabel('Km');
ylabel('Fm');
zlabel('real (log (min objective function))')
title('Real part of the log of the minimum objective function values')
set(gca,'FontSize',12)

disp(SF);

p= coeffvalues(SF);
% check the values
act = zeros(length(KmStart:KmIncr:KmEnd)*length(FmStart:FmIncr:FmEnd),1);
exp = act;
i = 0;
for Km = KmStart:KmIncr:KmEnd
    for Fm = FmStart:FmIncr:FmEnd
  
        i = i+1;
        thrLog = p(1) + p(2)*Km + p(3)*Fm + p(4)*Km^2 + p(5)*Km*Fm + p(6)*Fm^2;
        %thrLog = p(1) + p(2)*Fm + p(3)*Km + p(4)*Fm^2 + p(5)*Km*Fm + p(6)*Km^2;

        ePoint = 10^(thrLog + 1i*(pi/log(10)));
        act(i) = real(ePoint);
        exp(i) = zMin(i);            
        
    end
end

figure(6)
plot(real(act)-exp);

