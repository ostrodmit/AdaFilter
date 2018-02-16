function x = over_idft(xx, n)
if ~ismatrix(xx)
    error('Input must be a matrix');
end
if iscolumn(xx) 
    N = length(xx);
else
    N = sqrt(numel(xx));
end
if N < n
    error('N must be >= n');
end
x = idft(xx);
if iscolumn(xx)
    x = x(1:n);
else
    x = x(1:n,1:n);
end
end