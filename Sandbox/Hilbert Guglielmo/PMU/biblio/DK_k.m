function Dk = DK_k(k,N,theta) 
% Hp DK all ones starting at n different from 0. 
Dk = 1/N*sin(pi*k)./(sin(pi*k/N));
for i=1:length(k)
    if k(i)==0
        Dk(i) = 1;
    end
end
Dk = Dk.*exp(-1i*pi*k*(N-1)/N).*exp(1i*theta);
end