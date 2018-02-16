function x = recover_overblock(y,w,control,parallel,verb,neighb)
N = size(y);
if iscolumn(y), dim = 1; else dim = 2; end
nShifts = 16;
if min(N(1:dim)) < 3*w;
    error('Size must be at least 3*W');
end
if verb,
    disp(['Overlapping blocks. # of shifts: ' num2str(nShifts^(dim))]);
end
udShifts = floor(linspace(0,w,nShifts));
lrShifts = floor(linspace(0,w,nShifts));
if dim==1, lrShifts = 0; end
J = length(udShifts)*length(lrShifts);
rec = zeros(J,N(1),N(2));
fact = zeros(N(1),N(2)); % normalization factor when averaging
%base_rec = recover_blockwise(y,w,control,parallel,verb,neighb);
j=1;
for udShift = udShifts
    for lrShift = lrShifts
        %rec(j,:,:) = base_rec;
        rec(j,:,:)=0;
        rec(j,udShift+1:end-w+udShift,lrShift+1:end-w+lrShift)...
            = recover_blockwise(...
                y(udShift+1:end-w+udShift,lrShift+1:end-w+lrShift),...
                    w,control,parallel,verb,neighb);
        fact(udShift+1:end-w+udShift,lrShift+1:end-w+lrShift)...
            = fact(udShift+1:end-w+udShift,lrShift+1:end-w+lrShift)+1;
        j=j+1;
    end
end
x = squeeze(sum(rec(:,:,:),1)) ./ fact;