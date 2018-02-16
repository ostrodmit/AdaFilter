function y = gen_soft_threshold(z, lambda, p)
pos = subplus(abs(z) - lambda);
y = sign(z) .* power(pos, p);
end