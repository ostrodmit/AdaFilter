function x = idft(xx)
if ~ismatrix(xx)
    error('Input must be a matrix');
end
if iscolumn(xx) % do 1d-transform
    x = ifft(xx) * sqrt(length(xx));
else % 2x2 matrix, do 2d-transform
    x = ifft2(xx) * sqrt(numel(xx));
end
end