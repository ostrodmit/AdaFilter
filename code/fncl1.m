function sol=fncl1(data,lambda,control)
% For description of the algorithm, see .pdf file in the package.
% Solves the problem
% \min_x [\lambda*\|x\| + \|A*x-y\|_2^2,
% where x is complex-valued, f is a smooth convex function, \lambda >= 0, 
% and \|x\| is \ell_1 norm
% Note that b may be complex too (!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% data (structure): problem's data
%       data.n (integer): dimension of x;
%       data.y (complex array) array y 
%       data.FOracle handle to procedure with calling sequence 
%                    [f,df]=<...>(arg,data) which should return the value 
%                    (f) and the gradient (df) of f at the point arg
%       data.Lf (optional, positive real): \|\cdot\|_2-Lipschitz constant
%                                         of \nabla f(\cdot)
%       For sample data generator, see GetDataBlockEll1.m
%       
% lambda (double): parameter \lambda
% control (structure): control parameters
%        control.Eucl (char): 'y' for Euclidean proximal setup, any other 
%                              character for non-Euclidean setup
%        control.phistar (optional, double): true optimal value
%        control.xstar (optinal, 1D double array): true optimal solution        
%        control.ItrMax (integer): # of steps 
%        control.print (integer): zero value suppress printout, otherwise
%                                 some information is displayed every 
%                                 control.print steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% sol (structure):
%     sol.x (double 1D array): approximate solution
%     sol.obj (double): objective's value at the approximate solution
%     sol.Calls (integer): total number of computations of value and 
%                          derivative of f(A*x-b)     
%     sol.cpu (double): CPU time, sec
%     sol.Steps (integer): # of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrep_Lf=100;
if isfield(control,'Eucl'),
    if ischar(control.Eucl), 
        control.Eucl=strcmpi(control.Eucl(1),'y');
    end
else control.Eucl=1; 
end

scale=scale_block_l1(data,control.Eucl);
%
n=data.n;
%%%
sol.n=n;
if isfield(data,'Lf'), Lf=data.Lf; 
else % compute an estimation of the Lipschitz constant of the gradient
    Lf=0;
    if control.Eucl
        temp=randn(n,1)+sqrt(-1)*randn(n,1); temp=temp/norm(temp);
        for irep=1:nrep_Lf,
            %        f=full(data.A'*(data.A*e)); % compute A'*A*e on random e
            temp2=raha_oracle(temp,data);
            temp=raha_oracle(temp2,data,1);
            nt=norm(temp); temp=temp/nt;
            Lf=max(Lf,nt);
        end
    else
        for ibl=1:scale.n,
            tt=rand(1)+sqrt(-1)*rand(1); tt=tt/abs(tt);
            e=zeros(n,1);
            e(ibl)=tt;
            temp=raha_oracle(e,data);
            Lf_loc=norm(temp(:));
            Lf=max(Lf,Lf_loc);
        end
    end
end
sol.Lf=Lf;
%
y=zeros(n,1);
yp=y;
sol.x=y;
f=r_oracle(y,data);
sol.Calls=1;
sol.obj=f;
sol.Lmax=0;
sol.Lav=0;
Aold=0;
LinTerm=zeros(n,1);
ConstTerm=0;
L=Lf;
tstart=clock;
for Step=0:control.ItrMax,
    [zt,omegazt,domegazt]=prox_block_l1(LinTerm,Aold,lambda,scale);
    
        
    while(1==1)
        anew=(0.5/L)+sqrt(0.25/L^2+Aold/L);
        Anew=Aold+anew;
%
        taut=anew/Anew;
        xtp1=taut*zt+(1-taut)*yp;
        [fxtp1,dfxtp1]=r_oracle(xtp1,data);
        sol.Calls=sol.Calls+1;
        obj=fxtp1+lambda*sum(abs(xtp1));
        if obj<sol.obj,
            sol.obj=obj;
            sol.x=xtp1;
        end;
        %
        eta=anew*dfxtp1-domegazt;
        alpha=anew;
        [hatxtp1,omegahatxtp1]=prox_block_l1(eta,alpha,lambda,scale);
        ytp1=taut*hatxtp1+(1-taut)*yp;
        fytp1=r_oracle(ytp1,data);
        sol.Calls=sol.Calls+1;

        obj=fytp1+lambda*sum(abs(ytp1));
        if obj<sol.obj,
            sol.obj=obj;
            sol.x=ytp1;
        end;
        deltat=(omegahatxtp1-omegazt-(hatxtp1-zt)'*domegazt)/Anew...
            +(ytp1-xtp1)'*dfxtp1+fxtp1-fytp1;
        if deltat>=-1.e-6/Anew,
            Aold=Anew;
            LinTerm=LinTerm+anew*dfxtp1;
            ConstTerm=ConstTerm+anew*fxtp1-anew*dfxtp1'*xtp1;
            sol.Lmax=max(sol.Lmax,L);
            sol.Lav=(Step*sol.Lav+L)/(Step+1);
            break;
        else
            L=2*L;
        end;
    end;
    yp=ytp1;
    if control.print,
        if mod(Step,control.print)==0,
            tused=etime(clock,tstart);
            fprintf('%5.2f Steps=%3d Calls=%3d obj=%e L=%.3e\n',...
                tused,Step+1,sol.Calls,sol.obj,L);
        end;
    end;
    L=L/2;
end;
sol.Steps=Step;
tused=etime(clock,tstart);
sol.cpu=tused;
sol.Eucl=control.Eucl;
% 
sol.lambda=lambda;
%
sol.nnz=0;
nrm=norm(sol.x);
sol.nnz=sum(abs(sol.x)>1.e-6*nrm);
end % endof fncl1

