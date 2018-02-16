function xx = dft(x)
if ~ismatrix(x)
    error('Input must be a matrix');
end
if iscolumn(x) % do 1d-transform
    xx = fft(x) / sqrt(length(x));
else % square matrix, do 2d-transform
    xx = fft2(x) / sqrt(numel(x));
end
end