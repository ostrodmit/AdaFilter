function u = direct_operator(phi, rho, y)
% efficiently computes primal operator A(phi)
if iscolumn(phi) && iscolumn(y) % column vectors
    T = length(phi)-1; % phi is (T+1)x1
    L = length(y)-length(phi); % y is (T+L+1)x1
    pphi = [idft(phi); zeros(L,1)]; % zero-padding
    yy = [y(T+1:T+L+1); y(1:T)]; % permutation
    uu = idft( dft(pphi) .* dft(yy) ); % fast circular convolution
    uu = uu(1:L+1); % truncation
    u = rho * sqrt(T+L+1) / sqrt(T+1) * dft(uu); % DFT
elseif ismatrix(phi) && ismatrix(y) ... % square matricess
        && prod(size(phi) == size(phi')) && prod(size(y) == size(y'))
    T = length(phi)-1; % phi is (T+1)x(T+1)
    L = length(y)-length(phi); % y is (T+L+1)x(T+L+1)
    pphi = [idft(phi) zeros(T+1, L); zeros(L, T+L+1)]; % zero-padding
    yy = [y(T+1:T+L+1, T+1:T+L+1) y(T+1:T+L+1, 1:T); ...
          y(1:T, T+1:T+L+1) y(1:T, 1:T)]; % permutation
    uu = idft( dft(pphi) .* dft(yy) ); % fast circular convolution
    uu = uu(1:L+1, 1:L+1); % truncation
    u = rho * (T+L+1) / (T+1) * dft(uu); % DFT
else
    error('Incorrect input');
end
end