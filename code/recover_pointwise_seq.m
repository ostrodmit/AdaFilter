function x = recover_pointwise_seq(y,T,control,verb,neighb)
% works both with 1d and 2d signals
% parallel
N = size(y);
if iscolumn(y), dim = 1; else dim = 2; end
if min(N(1:dim)) < 4*T
    error('Size must be at least 4*W');
end
result = zeros(N(1)*N(2),1);
if verb 
    disp(['Pointwise mode.  W = ' num2str(T) '.']);
    wb = waitbar(0,'Processing...','WindowStyle','modal');
end
for idx = 1:N(1)*N(2)
    % get a 2-d pixel index
    i1 = idivide(idx-1,uint32(N(2)))+1;
    i2 = mod(idx-1,N(2))+1;
    i = [i1;i2];
    y1_t=zeros(dim,4*T); y2_t=zeros(dim,2*T); y3_t=zeros(dim,2*T+1);
    zone_code = zeros(dim,1);
    for ax = 1:dim
        zone_code(ax) = -(i(ax)<=T)+(i(ax)>=N(ax)-T+1)...
                        -(i(ax)<=2*T)+(i(ax)>=N(ax)-2*T+1);
        if zone_code(ax)==-2 % leftmost
            y1_t(ax,:) = 1:4*T;
            y2_t(ax,:) = i(ax):i(ax)+2*T-1;
            y3_t(ax,:) = 1:2*T+1;
        elseif zone_code(ax)==-1 % left
            y1_t(ax,:) = 1:4*T;
            y2_t(ax,:) = T+1:3*T;
            y3_t(ax,:) = i(ax)-T:i(ax)+T;
        elseif zone_code(ax)==0 % middle
            y1_t(ax,:) = i(ax)-2*T+1:i(ax)+2*T;
            y2_t(ax,:) = i(ax)-T+1:i(ax)+T;
            y3_t(ax,:) = i(ax)-T:i(ax)+T;
        elseif zone_code(ax)==1 % right
            y1_t(ax,:) = N(ax)-4*T+1:N(ax);
            y2_t(ax,:) = N(ax)-3*T+1:N(ax)-T;
            y3_t(ax,:) = i(ax)-T:i(ax)+T;
        elseif zone_code(ax)==2 % rightmost
            y1_t(ax,:) = N(ax)-4*T+1:N(ax);
            y2_t(ax,:) = i(ax)-2*T+1:i(ax);
            y3_t(ax,:) = N(ax)-2*T:N(ax);
        end
    end
    if dim == 1
        y1 = y(y1_t(1,:),1);
        y2 = y(y2_t(1,:),1);
        y3 = y(y3_t(1,:),1);
    else
        y1 = y(y1_t(1,:),y1_t(2,:));
        y2 = y(y2_t(1,:),y2_t(2,:));
        y3 = y(y3_t(1,:),y3_t(2,:));
    end
    sol = fit_filter(y1,y2,control);
    if sol.status, 
        fprintf('Problem unsolved. Status: %d\n', sol.status);
    end
    if dim == 1
        result(idx) = conv(y3, sol.filter,'valid');
    else
        result(idx) = conv2(y3, sol.filter,'valid');
    end
    if neighb && idx > 1 % initialize by the previous result
        % phase adjustment
        phase = ones(2*T+1,dim);
        for ax = 1:dim
            if zone_code(ax)==-2
                phase(:,ax) = exp(2*pi*1i*(0:2*T)'./(2*T+1));
            elseif zone_code(ax)==2
                phase(:,ax) = exp(2*pi*1i*(0:2*T)'./(2*T+1));
            end
        end
        Phase2d = phase(:,1);
        if dim == 2, Phase2d = Phase2d * phase(:,2).'; end
        %Phase2d = ones(2*T+1,dim);
        control.phi0 = sol.phi .* Phase2d;
        control.u0 = sol.u;
    end
    %sol.iter
    if verb, waitbar(idx/(N(1)*N(2))); end
end
if verb, close(wb); end
x = conj(reshape(result,N(2),N(1))');
end