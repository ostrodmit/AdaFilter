function xx = over_dft(x, N)
if ~ismatrix(x)
    error('Input must be a matrix');
end
if iscolumn(x) 
    n = length(x);
else
    n = sqrt(numel(x));
end
if N < n
    error('N must be >= n');
end
if iscolumn(x) % do 1d-transform
    x0 = [x; zeros(N-n,1)];
    xx = fft(x0) / sqrt(N);
else % square matrix, do 2d-transform
    x0 = [x zeros(n, N-n); zeros(N-n,n) zeros(N-n,N-n)];
    xx = fft2(x0) / N;
end
end