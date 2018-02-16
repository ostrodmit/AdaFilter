function Ax=raha_oracle(x,rdata,h_ind)
% multiplication oracle for the  Recht et Co stuff
ni=nargin;
if ni<3, h_ind=0;end
if (ni==3)&&ischar(h_ind)
    h_ind=strcmpi(h_ind,'y');
end

if ~h_ind
    if length(rdata.lsize)==1 %scalar signal
        temp=ifft(x)*sqrt(rdata.lsize);
        Ax=temp(1:rdata.ssize);
    else
        temp=reshape(x,rdata.lsize(1),rdata.lsize(2));
        temp=ifft2(temp)*sqrt(rdata.lsize(1)*rdata.lsize(2));
        Ax=temp(1:rdata.ssize(1), 1:rdata.ssize(2));
    end
else
    if length(rdata.lsize)==1 %scalar signal
        temp=[x;zeros(rdata.lsize-rdata.ssize,1)];
        Ax=fft(temp)/sqrt(rdata.lsize);
    else
        temp=zeros(rdata.lsize); 
        temp(1:rdata.ssize(1),1:rdata.ssize(2))=x;
        temp=fft2(temp)/sqrt(rdata.lsize(1)*rdata.lsize(2));
        Ax=reshape(temp,rdata.lsize(1)*rdata.lsize(2),1);
    end
end
end % end of raha_oracle.m