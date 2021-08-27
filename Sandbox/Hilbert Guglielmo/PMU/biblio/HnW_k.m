function Wk = HnW_k(k,N,theta)

Wk = -0.25*DK_k(k-1,N,theta) + 0.5*DK_k(k,N,theta) - 0.25*DK_k(k+1,N,theta);

end