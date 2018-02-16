function phi = adjoint_operator(u, rho, y)
% efficiently computes adjoint operator A^*(u)
if iscolumn(u) && iscolumn(y) % column vectors
    L = length(u)-1; % u is (L+1)x1
    T = length(y)-length(u); % y is (T+L+1)x1
    uu = [idft(u); zeros(T,1)]; % zero-padding
    yy = [y(T+1:T+L+1); y(1:T)]; % permutation
    pphi = idft( conj(dft(yy)) .* dft(uu) ); % fast circular convolution
    pphi = pphi(1:T+1); % truncation
    phi =  rho * sqrt(T+L+1) / sqrt(T+1) * dft(pphi); % DFT
elseif ismatrix(u) && ismatrix(y) ... % square matrices
        && prod(size(u) == size(u')) && prod(size(y) == size(y'))
    L = length(u)-1; % u is (L+1)x(L+1)
    T = length(y)-length(u); % y is (T+L+1)x(T+L+1)
    uu = [idft(u) zeros(L+1, T); zeros(T, T+L+1)]; % zero-padding
    yy = [y(T+1:T+L+1, T+1:T+L+1) y(T+1:T+L+1, 1:T); ...
          y(1:T, T+1:T+L+1) y(1:T, 1:T)]; % permutation
    pphi = idft( conj(dft(yy)) .* dft(uu) ); % fast circular convolution
    pphi = pphi(1:T+1, 1:T+1); % truncation
    phi =  rho * (T+L+1) / (T+1) * dft(pphi); % DFT
else
    error('Incorrect input');
end
end