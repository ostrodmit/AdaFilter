function Lip = lip_const_upper_bound_1_p(y1,T,L,rho,pic)
% upper bound on |A|_{1,p}, p >= 2
    Y1 = dft(y1);
    Lip = rho * sqrt((T+L+1)/(T+1))^(pic+1) * norm(Y1(:),inf);
end