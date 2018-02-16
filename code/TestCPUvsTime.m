%% Generate data
diary on;
diary 'History.txt'
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
data.obs=zeros(m,n);
for ir=1:r,
    e=randn(m,1);
    f=randn(1,n);
    data.obs=data.obs+e*f;
end;
data.p = p;
data.n = n;
data.m = m;
data.r = r;
%data.R=R;
data.mask=double(rand(m,n)<=p);
data.mult = @(x,Tflag,data)PlainMult(x,Tflag,data,'partial');
noise = input('noise level>');
data.noise = noise*mean(mean(abs(data.obs)));
data.obs = data.mult(data.obs+ data.noise*randn(m,n),1,data);
data.reg = data.noise/10;
data.R_nuc=10000;
data

%% Compare function value vs time
control.PRINT = 1;
control.PRINT_STEP = 20;
control.RELTOL = 1e-3;
control.MAXCPU=200;

%% CMP
control.PROX_TYPE ='exact';
control.MAXITER=5000;
Sol_CMP = CoMP_CG(data,control)

%% CMP_CG
control.MAXITER=1000;
control.PROX_TYPE = 'CondG';
control.CGMaxSteps = 12;
t = 1:(control.MAXITER+1);
control.PRINT_STEP = 1;
control.CGAccuracy = 1e-4*t./t;
Sol_CMP_CG = CoMP_CG(data,control)

%% Tune Smoothness
% Tune the smooth parameter through trial-and-error method
smoothness=[100,10,1,0.1,0.01];
control.MAXITER=100;
control.PRINT = 0;
fprintf('Searching for best smooth parameter...\n')
for trial = 1:5;
    control.SMOOTH = smoothness(trial);
    fprintf('Choose smoothness = %5.3f\n',control.SMOOTH);
    Sol_Smooth_CG = Smoothing_CG(data,control);
    fprintf('Objective after 100 iterations is %5.3f\n', Sol_Smooth_CG.opt);
    obj(trial)=Sol_Smooth_CG.opt;
end
[void,best]=max(obj);
best_smooth = smoothness(best);
fprintf('Found best smooth parameter %5.3f\n',best_smooth)

%% Smooth-CG
control.SMOOTH = best_smooth;
control.MAXITER=1000;
control.PRINT = 1;
Sol_Smooth_CG = Smoothing_CG(data,control)


%% APG
control.PROX_TYPE = 'exact';
control.MAXITER=5000;
control.STEPSIZE='backtracking';
 Sol_APG = Smoothing_AccProxGrad_CG(data,control)
 
%% APG-CG
control.PROX_TYPE = 'CondG';
control.MAXITER=1000;
control.CGMaxSteps = 12;
t = 1:(control.MAXITER+1);
control.PRINT_STEP = 1;
control.CGAccuracy = 1e-5*t./t;
Sol_APG_CG = Smoothing_AccProxGrad_CG(data,control)


%% Time
TimeMax = control.MAXCPU;
TimeMax = ceil(TimeMax/100)*100;
items = 21;
intervals = TimeMax/(items-1);
t = 1:1:items;
times = intervals*(t-1);

%% CMP
CMP_objs= zeros(items,1);
CMP_objs(1) = Sol_CMP.obj_init;
i0 = 1;
for j = 1:items
    for i = i0:length(Sol_CMP.objvals)
        if Sol_CMP.tocvals(i)<times(j)
            CMP_objs(j) = Sol_CMP.objvals(i);
        else
            i0=i;
            break;
        end
    end
end

%% CMP-CG
CMP_CG_objs= zeros(items,1);
CMP_CG_objs(1) = Sol_CMP_CG.obj_init;
i0 = 1;
for j = 1:items
    for i = i0:length(Sol_CMP_CG.objvals)
        if Sol_CMP_CG.tocvals(i)<=times(j)
            CMP_CG_objs(j) = Sol_CMP_CG.objvals(i);
        else
            i0=max(i-1,1);
            break;
        end
    end
end

%% Smooth-CG
Smooth_CG_objs= zeros(items,1);
Smooth_CG_objs(1) = Sol_Smooth_CG.obj_init;
i0 = 1;
for j = 1:items
    for i = i0:length(Sol_Smooth_CG.objvals)
        if Sol_Smooth_CG.tocvals(i)<times(j)
            Smooth_CG_objs(j) = Sol_Smooth_CG.objvals(i);
        else
            i0=i;
            break;
        end
    end
end

%% APG
APG_objs= zeros(items,1);
APG_objs(1) = Sol_APG.obj_init;
i0 = 1;
for j = 1:items
    for i = i0:length(Sol_APG.objvals)
        if Sol_APG.tocvals(i)<times(j)
            APG_objs(j) = Sol_APG.objvals(i);
        else
            i0=i;
            break;
        end
    end
end

%% APG-CG
APG_CG_objs= zeros(items,1);
APG_CG_objs(1) = Sol_APG_CG.obj_init;
i0 = 1;
for j = 1:items
    for i = i0:length(Sol_APG_CG.objvals)
        if Sol_APG_CG.tocvals(i)<times(j)
            APG_CG_objs(j) = Sol_APG_CG.objvals(i);
        else
            i0=i;
            break;
        end
    end
end

%% Function value vs Time
figure;
if isfield(data,'opt')
plot(ones(length(times),1)*data.opt,'k','Linewidth',2);
end
hold on
plot(times, CMP_objs,'r--','MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(times, CMP_CG_objs,'r','MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(times, Smooth_CG_objs,'b','MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(times, APG_objs,'g--','MarkerSize', 10, 'LineWidth', 1.5)
hold on
plot(times, APG_CG_objs,'g','MarkerSize', 10, 'LineWidth', 1.5)
xlabel('Elapsed Time')
ylabel('Function value')
if isfield(data,'opt')
    legend('true','CMP','CMP-CG','Smooth-CG','APG','APG-CG')
else
    legend('CMP','CMP-CG','Smooth-CG','APG','APG-CG')
end 
box on;


%% Function Accuracy vs Time
if isfield(data,'opt')
    figure;
    semilogy(times, CMP_objs-data.opt,'m','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    semilogy(times, CMP_CG_objs-data.opt,'r','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    semilogy(times, Smooth_CG_objs-data.opt,'b','MarkerSize', 10, 'LineWidth', 1.5)
    xlabel('Elapsed Time')
    ylabel('Objective Accuracy')
    legend('CMP','CMP-CG','Smooth-CG')
    box on;
    
    %%
    times(1)=1;
    figure;
    data.opt=Sol_CMP.opt-5;
    loglog(times, CMP_objs-data.opt,'r--','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    loglog(times, CMP_CG_objs-data.opt,'r','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    loglog(times, Smooth_CG_objs-data.opt,'b','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    loglog(times, APG_objs-data.opt,'g--','MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    loglog(times, APG_CG_objs-data.opt,'g','MarkerSize', 10, 'LineWidth', 1.5)
    xlabel('Elapsed Time')
    ylabel('Objective Accuracy')
    legend('CMP','CMP-CG','Smooth-CG','APG-CG')
    box on;
    
end

diary off;