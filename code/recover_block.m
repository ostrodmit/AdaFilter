function [rec,new_initial] = recover_block(y,pointwise,T,control,initial)
rec = zeros(size(y));
if iscolumn(y), dim = 1; else dim = 2; end
for ifFliplr = 0:(dim-1)
    for ifFlipud = 0:1
        yflip = y;
        if ifFliplr, yflip = fliplr(yflip); end
        if ifFlipud, yflip = flipud(yflip); end
        part_rec = zeros(size(y));
        if ~pointwise
            if dim == 1
                y1 = yflip(1:2*T,1);
                y2 = yflip(T+1:2*T,1);
            else
                y1 = yflip(1:2*T,1:2*T);
                y2 = yflip(T+1:2*T,T+1:2*T);
            end
            if ~isempty(initial)
                control.phi0 = initial.phi{ifFliplr+1,ifFlipud+1};
                control.u0 = initial.u{ifFliplr+1,ifFlipud+1};
            end
            sol = fit_filter(y1,y2,control);
%             if sol.status, 
%                 fprintf('Problem unsolved. Status: %d\n', sol.status);
%             end
            new_initial.phi{ifFliplr+1,ifFlipud+1} = sol.phi;
            new_initial.u{ifFliplr+1,ifFlipud+1} = sol.u;
            if dim == 1
                part_rec(T+1:2*T,1) = conv(y1, sol.filter, 'valid');
            else
                part_rec(T+1:2*T,T+1:2*T) = conv2(y1, sol.filter, 'valid');
            end
%         else % currently not used: process pointwise inside 4Tx4T blocs
%             for i = (1:2*T)
%                 for j = (1:2*T)
%                     y1 = yflip((0:2*T-1)+i,(0:2*T-1)+j);
%                     y2 = yflip((T:2*T-1)+i,(T:2*T-1)+j);
%                     sol = fit_filter(y1,y2,control);
%                     if sol.status, disp(sol.status); end
%                     quart_rec(2*T+i,2*T+j) = conv2(...
%                         yflip((T:2*T)+i,(T:2*T)+j),...
%                         sol.filter, 'valid');
%                     if verb
%                         progress = 100*(0.25*((i-1)*2*T+j)/(2*T)^2 ... 
%                             + 0.5 * ifFliplr + 0.25 * ifFlipud);
%                         curCent = idivide(progress,uint32(1));
%                         if  curCent>0 && ~cents(curCent)
%                             fprintf('*');
%                             cents(curCent) = 1;
%                         end
%                     end
%                 end
%             end
        end
        if ifFlipud, part_rec = flipud(part_rec); end
        if ifFliplr, part_rec = fliplr(part_rec); end
        rec = rec + part_rec;
    end
end
end