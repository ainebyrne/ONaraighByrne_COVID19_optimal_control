function [t,y,dydt,penalty_temp,penalty]=ode_solve_synthetic(u,synthetic_case_data)

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

my_length=length(t);

confirmed_cases_x=synthetic_case_data(1,:);
deaths_x=synthetic_case_data(2,:);
d_dt_confirmed_x=synthetic_case_data(3,:);
d_dt_deaths_x=synthetic_case_data(4,:);

%**************************************************************************

data=reshape(confirmed_cases_x,1,my_length);
theory=reshape(y(:,confirmed_index),1,my_length);

norm1=(1/my_length)*sum( ( log(data)-log(theory) ).^2);
% norm1= ( log(data(end))-log(theory(end)) )^2;

%**************************************************************************

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
nn=length(d_dt_deaths_x);

if(ix==1)
    d_dt_deaths_x(ix)=2*d_dt_deaths_x(ix+1)-d_dt_deaths_x(deaths(ix+2));
elseif(ix==nn)
    d_dt_deaths_x(ix)=2*d_dt_deaths_x(nn-1)-d_dt_deaths_x(nn-2);
else
    d_dt_deaths_x(ix)=(d_dt_deaths_x(ix+1)+d_dt_deaths_x(ix-1))/2;
end

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