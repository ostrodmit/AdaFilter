function u = explicit_operator(phi, rho, y)
% computes the primal operator A(psi) in a brute-force manner in O(TL)
if iscolumn(phi) && iscolumn(y)
else
    error('Inputs must be column vectors')
end
phi = ifft(phi) * sqrt(length(phi));
uu = conv(y, phi, 'valid');
u = fft( rho * uu ) / length(phi);
end