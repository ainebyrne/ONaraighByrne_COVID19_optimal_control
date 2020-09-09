function [t,y,dydt,case_data_x,penalty_temp,penalty]=ode_solve_seir_parameters_betas(u,case_data)

% Parameters:
% beta0,
% f,g,h,tau_L,tau_IP,tau_I,tau_D
% beta1,beta2,
% t_offset
% ip,ia,is,

beta0=u(1);
f=u(2);
g=u(3);
% h=1;
tau_L=u(4);
tau_IP=u(5);
tau_I=u(6);
tau_D=u(7);

beta1=u(8);
beta2=u(9);

t_offset=u(10);

ip=1;
ia=u(11);
is=u(12);

h=u(13);

tau_T=u(14);


confirmed_index=10;
deaths_index=8;

% *************************************************************************
% Population of country

N = 4.9*10^6;

% *************************************************************************
% Import case data

confirmed_cases=case_data(1,:);
recovered=case_data(2,:);
deaths=case_data(3,:);

% *************************************************************************
% First confirmed case: 29th Feb 2020 = day zero
% Schools closed: effective 13th March 2020 = day 13
% Pubs closed: 15th March 2020
% Lockdown, 2km limit:  effective 28th March 2020  = day 28
% Lockdown eased, 5km limit: 5th May 2020

% Integrate until April 30th 2020:
%tfinal =61;%12;
tfinal=77;

% Patient zero is introduced at time t=tstart.
% Day t=0 is 29th Februrary (first recorded case).
tstart = 0-t_offset;

x0=0*(1:confirmed_index);
x0(1)=N-1;
x0(3)=1;

% *************************************************************************
% ODE solver options

options = odeset('OutputSel',1,'Refine',4);

% *************************************************************************
% Solve ODE equations
sol = ode45(@SEIR_equations,[tstart,tfinal],x0,options,N,beta0,beta1,beta2,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D);


% Time in days:
t = 0:1:tfinal;

y = deval(sol,t)';
dydt=y;

for ii=1:length(t)
    dydt(ii,:)=SEIR_equations(t(ii),y(ii,:),          N,beta0,beta1,beta2,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D);
end

% *************************************************************************

confirmed_cases_x=confirmed_cases(1:length(t));
recovered_x     =recovered(1:length(t));
deaths_x         =deaths(1:length(t));

% requires curve-fitting toolbox:
confirmed_cases_x=smooth(confirmed_cases_x,1);
recovered_x =smooth(recovered_x ,1);
deaths_x  =smooth(deaths_x  ,1);

%**************************************************************************

my_length=length(confirmed_cases_x);
data=reshape(confirmed_cases_x,1,my_length);
theory=reshape(y(:,confirmed_index),1,my_length);

norm1=(1/my_length)*sum( ( log(data)-log(theory) ).^2);
% norm1= ( log(data(end))-log(theory(end)) )^2;

%**************************************************************************

d_dt_confirmed_x=confirmed_cases_x(2:end)-confirmed_cases_x(1:end-1);
d_dt_confirmed_x(end+1)=2*d_dt_confirmed_x(end)-d_dt_confirmed_x(end-1);

[~,ix]=max(d_dt_confirmed_x);
d_dt_confirmed_x(ix)=(d_dt_confirmed_x(ix+1)+d_dt_confirmed_x(ix-1))/2;

my_length=length(d_dt_confirmed_x);
data=reshape(d_dt_confirmed_x,1,my_length);
theory=reshape(dydt(:,confirmed_index),1,my_length);

norm2=(1/my_length)*sum( (log(data)-log(theory) ).^2);

%**************************************************************************

ctr=1;
while(deaths_x(ctr)==0)
    ctr=ctr+1;
end

my_length=length(t)-ctr+1;
data=10*reshape(deaths_x(ctr:end),1,my_length);
theory=10*reshape(y(ctr:end,deaths_index),1,my_length);

norm3=(1/my_length)*sum( ( log(data)-log(theory) ).^2);
% norm3= ( log(data(end))-log(theory(end)) )^2;


%**************************************************************************

d_dt_deaths_x=deaths_x(2:end)-deaths_x(1:end-1);
d_dt_deaths_x(end+1)=2*d_dt_deaths_x(end)-d_dt_deaths_x(end-1);

[~,ix]=max(d_dt_deaths_x);
d_dt_deaths_x(ix)=(d_dt_deaths_x(ix+1)+d_dt_deaths_x(ix-1))/2;

% Deaths start on day 20.
t_start=20;
t_end=min(length(t),75);

my_length=t_end-t_start+1;
data=10*reshape(d_dt_deaths_x(t_start:t_end),1,my_length);
theory=10*reshape(dydt(t_start:t_end,deaths_index),1,my_length);

norm4=(1/my_length)*sum( (log(data)-log(theory) ).^2);

%**************************************************************************


penalty_temp=norm1+norm2+norm3+norm4;

% [mx,ix]=max(dydt(:,10));
% t_peak=t(ix);

% display(t_peak)

% temp_vec=dydt(:,5);
% max_temp=max(temp_vec(40:end));

% if( (beta2>beta1)||(beta1>beta0))
%     penalty=1e5;
% elseif(max_temp>0)
%     penalty=1e5;
% else
%     penalty=penalty_temp;
% end
% 

penalty=penalty_temp; 


% *************************************************************************

case_data_x=zeros(3,length(confirmed_cases_x));
case_data_x(1,:)=confirmed_cases_x;
case_data_x(2,:)=recovered_x;
case_data_x(3,:)=deaths_x;

% *************************************************************************

    function y = SEIR_equations(t,x,N,beta0,beta1,beta2,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D)
        % Equations for SEIR model, with multiple infected classes
        
        if t<13
            beta_t = beta0;
            h_t=h;
        elseif t<28
            beta_t = beta1;
            h_t=1;
        else
            beta_t = beta2;
            h_t=1;
        end
       
       
        % *****************************************************************
        % S, susceptible
        y(1) = -(beta_t/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5));
        
        % *****************************************************************
        % E, exposed
        y(2) =  (beta_t/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5))-(1/tau_L)*x(2);
        
        % *****************************************************************
        % I, infectious
        
        % I_p, pre-symptomatic
        y(3) =  (1/tau_L)*x(2)-(1/tau_IP)*x(3);
        
        % I_a, asymptomatic, infectious
        y(4) =      f*(1/tau_IP)*x(3)-(1/tau_I)*x(4);
        % I_s, sympomatic, infectious
        y(5) =  (1-f)*(1/tau_IP)*x(3)-(1/tau_I)*x(5);
       
        % *****************************************************************
        % Removed / immune:
        
        y(6)=(1/tau_I)*x(4)+((1-g)/tau_I)*x(5);
        
        % *****************************************************************
        % Dying
        
        y(7)=(g/tau_I)*x(5)-(1/tau_D)*x(7);
        
        % *****************************************************************
        % Dead
        
        y(8)=(1/tau_D)*x(7);
        
        % *****************************************************************
        % Awaiting test
        y(9)=h_t*(1-f)*(1/tau_IP)*x(3)-(1/tau_T)*x(9);
        
        y(10)=(1/tau_T)*x(9);
        
        y = y';
        
    end

% *************************************************************************

end