function x = recover_lepski(y,sigm,mode,control,parallel,verb,neighb)
% mode == 0 - pixelwise, mode == 1 - blocks, mode == 2 - overlapping blocks
% alpha = 0.2; % confidence = 1 - alpha
ratio = 1.5; % NIPS arxiv paper
ratio = 2; % algo-optfilter
if iscolumn(y), dim = 1; else dim = 2; end
N = size(y);
% producing bandwidths
bws = [];
w = 2;
if strcmp(control.accuracy,'rel')
    w = 6; % in small samples the disrepancy may be 0 sometimes
end
if mode == 0
    wmax = min(N(1:dim))/4;
elseif mode == 1
    wmax = min(N(1:dim))/2;
elseif mode == 2
    wmax = min(N(1:dim))/3;
else 
    error('Incorrect mode.');
end
while w <= wmax
    bws = [bws;w];
    w = ceil(w*ratio);
end
wmax = w;
if mode < 2
    rec = zeros(length(bws),N(1),N(2));
else
    rec = zeros(2,dim,length(bws),N(1),N(2));
end
% eps = sigm*(1+control.rho)*(sqrt(log(wmax+1)))./(sqrt(bws+1)).^(dim); % concervative 
eps = sigm*(control.rho)*(sqrt(log(bws+1))) ./(sqrt(bws+1)).^(dim)/2; % heuristic, used in arxiv (Brodatz)
eps = sigm*(control.rho)*(sqrt(log(wmax+1))) ./(sqrt(bws+1)).^(dim)/2; % heuristic, used in algo-optfilter

if verb
    disp(['Lepski recovery. Bandwidths: ' num2str(bws')])
end
for k = 1:length(bws)
    w = bws(k);
    if mode == 0
        rec(k,:,:) = recover_pointwise(y,w,control,parallel,verb,neighb);
    elseif mode == 1
        rec(k,:,:) = recover_blockwise(y,w,control,parallel,verb,neighb);
    else
        offset = w; % half a block size
        for ifDown = 0:1
            for ifRight = 0:dim-1
                if ~ifDown && ~ifRight
                    rec(1,1,k,:,:) = zeros(N(1),N(2));
                else 
                    rec(ifDown+1,ifRight+1,k,:,:) = rec(1,1,k,:,:);
                end
                rec(ifDown+1,ifRight+1,k,...
                    offset*ifDown+1:end,offset*ifRight+1:end)...
                    = recover_blockwise(...
                        y(offset*ifDown+1:end,...
                          offset*ifRight+1:end),...
                        w,control,parallel,verb,neighb);
            end
        end
    end
end
x = zeros(size(y));
% counter = zeros(length(bws),1);
% it = 0;
for t1 = 1:N(1)
    for t2 = 1:N(2)
        % find the largest admissible bandwidth
        k = 1;
        if mode < 2
            admissible = 1;
        else
            admissible = ones(2,dim);
            precision = ones(2,dim);
            K = ones(2,dim);% largest admissible
        end
        while sum(sum(admissible)) && k <= length(bws)
            k = k + 1;
            if mode < 2
                x(t1,t2) = rec(k-1,t1,t2);
                if k <= length(bws)
                    admissible = prod(abs(rec(k,t1,t2)-rec(1:k,t1,t2))...
                        <= eps(k) + eps(1:k));
                else
                    admissible = 0;
                end
            else
                for ifDown = 0:1
                    for ifRight = 0:dim-1
                        if admissible(ifDown+1,ifRight+1)
                            K(ifDown+1,ifRight+1) = k-1;
                            precision(ifDown+1,ifRight+1) = 1/eps(k-1)^2;
                        end
                        if k <= length(bws)
                            admissible(ifDown+1,ifRight+1) = prod(...
                                squeeze(abs(rec(ifDown+1,ifRight+1,1:k,t1,t2)...
                                    -rec(ifDown+1,ifRight+1,k,t1,t2)))...
                                <= eps(k) + eps(1:k));
                        else
                            admissible(ifDown+1,ifRight+1) = 0;
                        end
                    end
                end
            end
        end
        if mode == 2
%             if sum(range(precision))
%                 it = it+1;
%             end
            x(t1,t2) = 0;
            for ifDown = 0:1
                    for ifRight = 0:dim-1
                        x(t1,t2) = x(t1,t2) + precision(ifDown+1,ifRight+1)...
                            * rec(ifDown+1,ifRight+1,K(ifDown+1,ifRight+1),t1,t2);
                    end
            end
            x(t1,t2) = x(t1,t2) / sum(sum(precision));
        end
        % counter(k-1)=counter(k-1)+1;
    end
end
% counter
% it
end