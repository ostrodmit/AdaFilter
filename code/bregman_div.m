function D = bregman_div (constrained, x, y, r, gamma)
% the Bregman divergence of two complex vectors
    if r==2 && gamma == 1;
        D = norm(x-y,2)^2/2;
    else
        constrained = 0;
        if constrained % based on |*|_r^r
            D = (norm(x,r)^r/r-norm(y,r)^r/r...
                 -real(dot(sign(y).*abs(y).^(r-1),x-y)))/gamma;
        else % based on |*|_r^2
%             nonzeros(x)
%             nonzeros(y)
%             nonzeros(x-y)
%             r
%             gamma
            D = (norm(x,r)^2/2 - norm(y,r)^2/2 ...
                 - norm(y,r)^(2-r)*real(dot(sign(y).*abs(y).^(r-1),x-y)))...
                 * length(y)^((r-1)*(2-r)/r)/gamma;
            %if D < 0, error(num2str(D)); end
        end
    end
end