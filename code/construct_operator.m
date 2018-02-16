function A = construct_operator(T, rho, y)
if iscolumn(y) % 1d
    L = length(y) - (T+1);
    A = zeros(L+1,T+1);
    for i = 1:T+1
        % take the i-th cannonical basis vector
        e = zeros(T+1, 1);
        e(i) = 1;
        % apply the operator
        A(:,i) = direct_operator(e, rho, y);
    end
elseif ismatrix(y) && prod(size(y) == size(y')) % 2d
    L = length(y) - (T+1);
    A = zeros((L+1)^2, (T+1)^2);
    for i = 1:T+1
        for j = 1:T+1
            % take the (i,j)-th cannonical basis matrix
            e = zeros(T+1, T+1);
            e(i,j) = 1;
            % apply the operator
            Ae = direct_operator(e, rho, y);
            A(:,i+(T+1)*(j-1)) = Ae(:);
        end
    end
else
    error('Wrong input format');
end
end