function phi = init_warm_start(pic,y1,b,T,L,rho,p)
% warm start by the signal DFT for the primal solution
powSpect = abs(dft(y1(1:T+1, 1:(T+1)^pic)));
powSpect = powSpect(:);
[sortedPowSpect idx] = sort(powSpect, 'descend');
subDim = min(ceil(rho),(T+1)^(pic+1));
maxIdx = idx(1:subDim);
%phi(maxIdx) = 1/ceil(rho);
M = zeros((T+1)^(pic+1),subDim); % frequency mask
M(maxIdx,:) = eye(subDim);
% compute operator for least squares problem
A_sub = zeros((L+1)^(pic+1),subDim);
for j = 1:subDim
    M_j = M(:,j);
    if pic, M_j = vec2mat(M_j); end
    A_sub_j = direct_operator(M_j, rho, y1);
    A_sub(:,j) = A_sub_j(:);
end
% fit filter on fixed support
cvx_begin quiet
%cvx_solver mosek
cvx_precision best
    variable phi_sub(subDim) complex
    minimize ( norm (b(:) - A_sub * phi_sub, p) )
    subject to
        norm(phi_sub,1) <= 1;
cvx_end
phi = M * phi_sub;
if pic, phi = vec2mat(phi); end