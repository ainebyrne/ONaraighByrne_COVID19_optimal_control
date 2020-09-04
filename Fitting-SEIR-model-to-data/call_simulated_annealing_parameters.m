T = readtable('CovidStatisticsProfileHPSCIrelandOpenData.csv');
confirmed_cases=table2array(T(1:end,5));
recovered=table2array(T(1:end,9));
deaths=table2array(T(1:end,7));


case_data=zeros(3,length(confirmed_cases));
case_data(1,:)=confirmed_cases;
case_data(2,:)=recovered;
case_data(3,:)=deaths;



% Define the Penalty function with reference to the ODE-SIR model:
ObjectiveFunction=@(x)ode_solve_seir_parameters_betasX(x,case_data);

% First parameter  -  x(1)=beta0
% Second parameter -  x(2) = f
% Third parameter  -  x(3) = g
% Fourth parameter -  x(4) = tau_L
% Fifth parameter  -  x(5) = tau_IP
% Sixth parameter  -  x(6) = tau_I
% Seventh parameter-  x(7) = tau_D
% Eith parameter    - x(8) = beta1
% Ninth parameter   - x(9)= beta2
% Tenth parameter   - x(10)= t_offset

% Eleventh parameter -x(11)=ia;
% Twelth parameter   -x(12)=is;
% Thirteenth paramter-x(13)=symptomatic tested.

% Incubation period C - between 4.5 and 5.8 days, see
% https://www.acpjournals.org/doi/10.7326/M20-0504


% Lower and upper bounds on parameters
lb = [ 0, 0.1, 0.05 ,  3.7-0.65, 1.5-0.65,2.3-1,10, 0,0,1,   0,  0,  1, 1   ];
ub = [ 5, 0.5, 0.07 ,  3.7+0.65, 1.5+0.65,2.3+1,30, 5,5,15,  1,  1,  1, 10  ];

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

% options = optimoptions('fmincon','MaxFunctionEvaluations',50000,'Display','iter');
% %  
% [x,fval] = fmincon(ObjectiveFunction,u0,[],[],[],[],lb,ub,[],options);

% *************************************************************************


display(output)
display('optimal control policy:')
display(x)
display(strcat('minimum of cost function is...',num2str(fval)));

R0 = calculate_R0(x);

display('With case isolation, R0=')
display(R0)

R0x = calculate_R0([x(1:11),1,x(13)]);
display('Without case isolation, R0=')
display(R0x)

R0_schools=R0*(x(8)/x(1));
display('With case isolation and schools closed, R0=')
display(R0_schools)

R0_lockdown=R0*(x(9)/x(1));
display('With case isolation and lockdown, R0=')
display(R0_lockdown)

% *************************************************************************

% Plotting
[t,y,dydt,case_data_x,penalty_temp,penalty]=ode_solve_seir_parameters_betas(x,case_data);

confirmed_cases_x=case_data_x(1,:);
recovered_x=case_data_x(2,:);
deaths_x=case_data_x(3,:);

d_dt_confirmed_x=confirmed_cases_x(2:end)-confirmed_cases_x(1:end-1);
d_dt_confirmed_x(end+1)=2*d_dt_confirmed_x(end)-d_dt_confirmed_x(end-1);

d_dt_deaths_x=deaths_x(2:end)-deaths_x(1:end-1);
d_dt_deaths_x(end+1)=2*d_dt_deaths_x(end)-d_dt_deaths_x(end-1);

confirmed_index=10;
deaths_index=8;


subplot(2,2,1), plot(t,y(:,confirmed_index),'linewidth',3,'color','blue')
hold on
subplot(2,2,1), plot(t,confirmed_cases_x,'o');
xlim([0 length(t)])
hold off
%
%
subplot(2,2,2), plot(t,y(:,deaths_index),'linewidth',3,'color','red')
hold on
subplot(2,2,2), plot(t,deaths_x,'o','color','red')
%
%
subplot(2,2,3), plot(t,dydt(:,confirmed_index),'linewidth',3,'color','blue')
hold on
subplot(2,2,3), plot(t,d_dt_confirmed_x,'o','color','blue')
xlim([0 length(t)])
hold off
%
%
subplot(2,2,4), plot(t,dydt(:,deaths_index),'linewidth',3,'color','red')
hold on
subplot(2,2,4), plot(t,d_dt_deaths_x,'o','color','red')

drawnow