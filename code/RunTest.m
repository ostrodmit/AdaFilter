
%% Generate data
addpath('../minConf');
addpath('../minConf/minConf');
addpath('../minFunc_2012');
addpath('../minFunc_2012/minFunc');
fprintf('Generate Data...\n')
m = input('row dimension >');                      
n = input('column dimension >');                     
r = input('rank >');
p = input('observation density>');
%R = input('bound for nuclear norm>');
data.sol=zeros(m,n);
for ir=1:r,
    e=randn(m,1);
    f=randn(1,n);
    data.sol=data.sol+e*f;
end;
data.p = p;
data.n = n;
data.m = m;
data.r = r;
%data.R=R;
data.mask=double(rand(m,n)<=p);
data.mult = @(x,Tflag,data)PlainMult(x,Tflag,data,'partial');
noise = input('noise level>');
data.noise = noise*mean(mean(abs(data.sol)));
data.obs = data.mult(data.sol+ data.noise*randn(m,n),1,data);
data.reg = data.noise/10;
data.opt = -Inf;
data.R_nuc=10000;
clear data.sol;

%% Warm Start
s=input('Warm start? Yes[1]/No[0]>');
if s==1
    fprintf('Computing warm-up solution via smoothing CG...\n')
    control.SMOOTH = 1e-3;
    control.PRINT = 1;
    control.MAXITER = 50;
    Warm_Start = Smoothing_CG(data,control);
    initial.x= Warm_Start.x;
    initial.zeta= Warm_Start.zeta;
    fprintf('Initial solution generated with obj=%5.2f\n',Warm_Start.opt)
else
    initial.x = zeros(data.m, data.n);
    initial.zeta = 0;
end

%% Set control parameters
control.PRINT = 1;
control.PRINT_STEP = 1;
control.RELTOL = 1e-3;
control.MAXCPU=200;
warning('off','all')

%%
while(1==1)
    s = input('CVX[0]/FullProx[1]/SemiProx[2]/SmoothingCG[3]/SmoothingAPG[4]/Exit[-1]>');
    
    %% CVX
    if s==0  
        cvx_begin
        variable Z(data.m,data.n)
        minimize (norm(data.mask.*Z-data.obs,'fro') + data.reg*norm_nuc(Z))
        cvx_end
        data.R_fro=norm(Z,'fro');
        data.R_nuc = norm_nuc(Z);
        data.opt=cvx_optval;
        data.sol=Z;  
    end
       
    %% FullProx
    if s==1
        control.MAXITER = input('MaxIterations>');
        control.PROX_TYPE ='exact';
        Sol_CMP = CoMP_CG(data,control)
    end
    
    %% SemiProx CG
    if s==2
        control.PROX_TYPE = 'CondG';
        control.MAXITER = input('MaxIterations>');
        control.CGMaxSteps = input('MaxInnerConGSteps>');
        t = 1:(control.MAXITER+1);
        control.CGAccuracy = 1e-4*t./t;
        %control.CGAccuracy = 5./(sqrt(t)+5);
        Sol_CMP_CG = CoMP_CG(data,control)
        %Sol_CMP_CG_refine = CoMP_CG_Refine(data,control)
    end
    
    %% Smoothing CG
    if s==3
        control.MAXITER = input('MaxIterations>');
        control.SMOOTH = input('Smooth parameter>'); %{100,10,1,0.1,0.01,0.001}
        Sol_Smooth_CG = Smoothing_CG(data,control)
    end
    
     %% Smoothing APG-CG
    if s==4
        control.PROX_TYPE='exact';
        control.STEPSIZE='backtracking';
        control.MAXITER = input('MaxIterations>');
        control.SMOOTH = input('Smooth parameter>'); %{100,10,1,0.1,0.01,0.001}
        t = 1:(control.MAXITER+1);
        control.CGAccuracy = 0.2*t./t;   
        Sol_APG_CG = Smoothing_AccProxGrad_CG(data,control)
    end   
    %% Exit
    if s ==-1
        return;
    end
end

%% Compare
if isfield(data,'opt')
plot(ones(control.MAXITER,1)*data.opt,'k','Linewidth',2);
end
hold on;
plot([Sol_CMP.obj_init,Sol_CMP.objvals],'b','Linewidth',2)
hold on;
plot([Sol_CMP_CG.obj_init,Sol_CMP_CG.objvals],'r','Linewidth',2)
hold on;
plot([Sol_Smooth_CG.obj_init,Sol_Smooth_CG.objvals],'m','Linewidth',2)
hold on;
plot([Sol_APG_CG.obj_init,Sol_APG_CG.objvals],'g','Linewidth',2)


%%
loglog(([Sol_CMP.obj_init,Sol_CMP.objvals]-data.opt),'b','Linewidth',2)
hold on;
loglog(([Sol_CMP_CG.obj_init,Sol_CMP_CG.objvals]-data.opt),'r','Linewidth',2)
hold on;
loglog(([Sol_Smooth_CG.obj_init,Sol_Smooth_CG.objvals]-data.opt),'m','Lineidth',2)
hold on;
loglog([Sol_APG_CG.obj_init,Sol_APG_CG.objvals]-data.opt,'g','Linewidth',2)
xlabel('Number of Iterations')
ylabel('Objective Accuracy')
legend('CMP','CMP-CG','Smooth-CG','APG-CG')

