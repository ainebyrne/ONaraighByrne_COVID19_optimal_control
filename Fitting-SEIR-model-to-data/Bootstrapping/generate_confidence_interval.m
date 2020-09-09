

rng('default') % for reproducibility

T = readtable('CovidStatisticsProfileHPSCIrelandOpenData.csv');
confirmed_cases=table2array(T(1:end,5));
recovered=table2array(T(1:end,9));
deaths=table2array(T(1:end,7));


case_data=zeros(3,length(confirmed_cases));
case_data(1,:)=confirmed_cases;
case_data(2,:)=recovered;
case_data(3,:)=deaths;

% *************************************************************************

u0 = [1.4695
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
% Now I compute the parameters of the best-fit model:

u=my_nl_fit(case_data,u0);

% *************************************************************************
% Now I compute the model predictions arising from the best-fit model:
[t,y,dydt,case_data_x,penalty_temp,penalty]=ode_solve_seir_parameters_betas(u,case_data);

% I also extract the data from the relevant matlab codes:
confirmed_index=10;
deaths_index=8;

confirmed_cases_x=case_data_x(1,:);
deaths_x=case_data_x(3,:);

d_dt_confirmed_x=confirmed_cases_x(2:end)-confirmed_cases_x(1:end-1);
d_dt_confirmed_x(end+1)=2*d_dt_confirmed_x(end)-d_dt_confirmed_x(end-1);

d_dt_deaths_x=deaths_x(2:end)-deaths_x(1:end-1);
d_dt_deaths_x(end+1)=2*d_dt_deaths_x(end)-d_dt_deaths_x(end-1);


% *************************************************************************
% Now I compute the residuals.

residual_confirmed=log(confirmed_cases_x)-log(y(:,confirmed_index))';
residual_d_dt_confirmed=log(d_dt_confirmed_x)-log(dydt(:,confirmed_index))';

% Throw out the outlier.
% [~,ix]=max(abs(residual_d_dt_confirmed));
% residual_d_dt_confirmed(ix)=0;

% Deaths start on day 20.
temp0=[d_dt_deaths_x(20:end-1),1];
temp=log(temp0)-log(dydt(20:end,deaths_index))';


% *************************************************************************
% Now I do the bootstrapping.  I generate nbootstrap bootstrap samples.

nbootstrap=1000;

U=zeros(nbootstrap,length(u));

for i=1:nbootstrap
    
    residual_d_dt_confirmed_synthetic=randsample(residual_d_dt_confirmed,length(t),'true');
    d_dt_confirmed_synthetic=exp(residual_d_dt_confirmed_synthetic).*dydt(:,confirmed_index)';
    
    temp_synthetic=randsample(temp,length(temp),'true');
    residual_d_dt_deaths_synthetic=[0*(1:19),temp_synthetic];
    
    d_dt_deaths_synthetic=exp(residual_d_dt_deaths_synthetic).*dydt(:,deaths_index)';
    
    confirmed_synthetic=0*(1:length(t));
    deaths_synthetic=0*(1:length(t));
    
    % *********************************************************************
    
    my_seed0=randsample(residual_confirmed,1,'true');
    my_seed=y(1,confirmed_index)*exp(my_seed0);

    % *********************************************************************

    confirmed_synthetic(1)=my_seed;
    
    for ii=2:length(t)
        confirmed_synthetic(ii)=confirmed_synthetic(ii-1)+d_dt_confirmed_synthetic(ii);
        deaths_synthetic(ii)=deaths_synthetic(ii-1)+d_dt_deaths_synthetic(ii);
    end
    
    synthetic_case_data=zeros(4,length(t));
    synthetic_case_data(1,:)=confirmed_synthetic;
    synthetic_case_data(2,:)=deaths_synthetic;
    synthetic_case_data(3,:)=d_dt_confirmed_synthetic;
    synthetic_case_data(4,:)=d_dt_deaths_synthetic;
    
    u_temp=my_nl_fit1(synthetic_case_data,u0);
    U(i,:)=u_temp;
    
    display(i)
end

bootCI=prctile(U,[2.5,97.5]);


% *************************************************************************

function u=my_nl_fit(case_data,u0)
    lb = [ 0, 0.1, 0.05 ,  3.7-0.65, 1.5-0.65,2.3-1,10, 0,0,1,   0,  0,  1, 1   ];
    ub = [ 5, 0.5, 0.07 ,  3.7+0.65, 1.5+0.65,2.3+1,30, 5,5,15,  1,  1,  1, 10  ];
    ObjectiveFunction=@(z)ode_solve_seir_parameters_betasX(z,case_data);
    options = optimoptions('fmincon','MaxFunctionEvaluations',50000,'Display','iter'); 
    [u,~] = fmincon(ObjectiveFunction,u0,[],[],[],[],lb,ub,[],options);
end

function u=my_nl_fit1(synthetic_case_data,u0)
    lb = [ 0, 0.1, 0.05 ,  3.7-0.65, 1.5-0.65,2.3-1,10, 0,0,1,   0,  0,  1, 1   ];
    ub = [ 5, 0.5, 0.07 ,  3.7+0.65, 1.5+0.65,2.3+1,30, 5,5,15,  1,  1,  1, 10  ];
    ObjectiveFunction=@(z)ode_solve_synthetic_wrapper(z,synthetic_case_data);
    options = optimoptions('fmincon','MaxFunctionEvaluations',50000); 
    [u,~] = fmincon(ObjectiveFunction,u0,[],[],[],[],lb,ub,[],options);
end

