% *************************************************************************
% Optimum parameters from fitting to data:

params=[1.4695
0.3084
0.0599
3.7486
1.2938
2.1738
15.2210
1.1009
0.3576
9.9831
0.4096
0.3405
1.0000
3.8911];

% *************************************************************************
% Define the Penalty function with reference to the ODE-SIR model:
ObjectiveFunction=@ode_solve_seirX;

% Maximum realisable lockdown:
% umax=1-(params(9)/params(1));
% Round up to 0.8
umax=0.64;

% Full on bang-bang
n=12;

lb= zeros(1,2*n);
ub=[zeros(1,n)+umax,zeros(1,n)+1];

% u0=(lb+ub)/2;

u0=[    0.000346878727985
   0.412212365315978
   0.000536845525725
   0.185089252168063
   0.639997168503860
   0.639368516931177
   0.565739832938898
   0.629707648915150
   0.460049986869823
   0.293755948320455
   0.000000183696818
   0.321423635973558
   0.220998899157731
   0.013590459721689
   0.010547755001967
   0.000679022625668
   0.234134488130270
   0.678661121717337
   0.893507918321827
   0.000319664401889
   0.861169475644098
   0.818511888212091
   0.865523966517082
   0.000872048392506];


% *************************************************************************
% Simulated Annealing
% options = optimoptions('simulannealbnd','MaxFunctionEvaluations',120000,...
%         'FunctionTolerance',1e-6,'ReannealInterval',10,...
%         'PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf},...
%         'HybridFcn',{@fmincon,hybridopts});

hybridopts=optimoptions('fmincon','Display','iter','Algorithm','interior-point');
  
options = optimoptions('simulannealbnd','FunctionTolerance',1e-6,...
          'MaxFunctionEvaluations',1200000,'PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},...
          'HybridFcn',{@fmincon,hybridopts});
            
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,u0,lb,ub,options);

% *************************************************************************
% Global Search

% gs = GlobalSearch('Display','iter', ...
%     'StartPointsToRun','bounds-ineqs');
% 
% 
% problem = createOptimProblem('fmincon','x0',u0,...
%     'objective',ObjectiveFunction,'lb',lb,'ub',ub);
% 
% [x,fval] = run(gs,problem);

% *************************************************************************
% Surrogate Optimization

% mf=10000;
% options = optimoptions('surrogateopt','MaxFunctionEvaluations',mf,'PlotFcn','surrogateoptplot');
% [x,fval,~,~,pop] = surrogateopt(ObjectiveFunction,lb,ub,options);

% *************************************************************************
% Particle Swarm

% options = optimoptions('particleswarm','Display','iter');%,'SwarmSize',150);
% [x,fval] = particleswarm(ObjectiveFunction,length(lb),lb,ub,options);
% 
% display(fval)

% *************************************************************************
% Genetic Algorithm

% options = optimoptions('ga','PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});%,...
% %         % 'HybridFcn',{@fmincon,hybridopts});
% 
% [x,fval] = ga(ObjectiveFunction,n,[],[],[],[],lb,ub);
% 
% display(output)

% *************************************************************************
% Fmincon: 
% 
% options = optimoptions('fmincon','MaxFunctionEvaluations',50000,'Display','iter');
%  
% [x,fval] = fmincon(ObjectiveFunction,u0,[],[],[],[],lb,ub,[],options);

% *************************************************************************


display('optimal control policy:')
display(x)
display(strcat('minimum of cost function is...',num2str(fval)));

% *************************************************************************
% Plotting

[t,y,dydt,penalty]=ode_solve_seir(x);
ut_vec=plot_ut(x,t);

Is=y(:,5);

subplot(2,1,1), plot(t,0.016*Is,'linewidth',3,'color','blue')
hold on
subplot(2,1,1), plot(t,0*t+300,'linewidth',1,'color','black')
hold off
subplot(2,1,2), plot(t,ut_vec,'linewidth',3,'color','red')
drawnow