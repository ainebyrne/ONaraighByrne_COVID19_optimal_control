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
umax=0.8;

lb = [0,0,0];
ub = [umax,umax,1];

u0=(lb+ub)/2;

% *************************************************************************

%options = optimoptions('simulannealbnd','MaxFunctionEvaluations',10000,'PlotFcns',...
%          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});

hybridopts=optimoptions('fmincon','Display','iter','Algorithm','interior-point');
  
options = optimoptions('simulannealbnd','FunctionTolerance',1e-6,'PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},...
          'HybridFcn',{@fmincon,hybridopts});
            
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,u0,lb,ub,options);

% *************************************************************************

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