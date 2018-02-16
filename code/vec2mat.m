function X = vec2mat(x)
n = sqrt(length(x));
X = reshape(x,[n,n]);
end