function [t,y,dydt,penalty]=ode_solve_seir(u)

% Parameters:
% beta0,
% f,g,h,tau_L,tau_IP,tau_I,tau_D
% beta1,beta2,
% t_offset
% ip,ia,is,

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

beta0=params(1);
f=params(2);
g=params(3);
tau_L=params(4);
tau_IP=params(5);
tau_I=params(6);
tau_D=params(7);

% beta1=params(8);
% beta2=params(9);

t_offset=params(10);

ip=1;
ia=params(11);
is=params(12);

h=params(13);

tau_T=params(14);

% *************************************************************************
% Choose optimum control horizon to be one year.
tfinal=365;

% Controls:

nn=length(u);
n=nn/2;

u_control=u(1:n);

tau_vec=tfinal*u(n+1:end)/sum(u(n+1:end));

ts_vec=0*tau_vec;
for ii=1:n
    ts_vec(ii)=sum(tau_vec(1:ii));
end

% *************************************************************************
% Population of country
N = 4.9*10^6;

% Patient zero is introduced at time t=tstart.
% Day t=0 is 29th Februrary (first recorded case).
tstart = 0-t_offset;

x0=0*(1:10);
x0(1)=N-1;
x0(3)=1;

% *************************************************************************
% ODE solver options

options = odeset('OutputSel',1,'Refine',4);

% *************************************************************************
% Solve ODE equations

% Solve up until the implementation of some NPIss
sol =ode45(@SEIR_equations0,[tstart,13],x0,     options,N,beta0,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D);

t_temp=tstart:1:13;
y_temp=deval(sol,t_temp);
x0_temp=y_temp(:,end);

% Thereafter, explroe what is the optimal sequence of NPIs:
sol = ode45(@SEIR_equations1,[1 tfinal], x0_temp,options,N,beta0,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D);

% Time in days:
t = 1:1:tfinal;

y = deval(sol,t)';
dydt=y;

for ii=1:length(t)
    dydt(ii,:)=SEIR_equations1(t(ii),y(ii,:),          N,beta0,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D);
end

% *************************************************************************
% Working out the penalty function.

Is=y(:,5);
max_Is=max(Is);

penalty_temp= sum(u_control.*tau_vec);

z0=smooth(y(:,5),5);
for i=1:length(z0)
    if(z0(i)<0.1)
        z0(i)=0.1;
    end
end
z=smooth(z0);
TF = islocalmax(z);
zz=t(TF);
ctr=length(zz);

dy_val=dydt(end,5);

penalty_temp=penalty_temp+(tanh((0.016*max_Is-300)/0.1)+1)*( (0.016*max_Is-300)^2)...
   +(tanh(dy_val/0.1)+1)*(dy_val^2)...
    +0*100/ctr;

% dy_val=dydt(end,5);
% penalty_temp=sum( (0.016*Is-300).^2)+(tanh(dy_val/0.1)+1)*(dy_val^2);

%if((0.016*max_Is>300.1)||(tswitch3>tfinal)||(dydt(end,5)>0)||(ctr<3)||(ctr>5))
% if((dydt(end,5)>0))
%     penalty=penalty_max;
% else
%     penalty=penalty_temp;
% end

penalty=penalty_temp;

% *************************************************************************

    function y = SEIR_equations0(~,x,N,beta0,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D)
        % Equations for SEIR model, with multiple infected classes

        % *****************************************************************
        % S, susceptible
        y(1) = -(beta0/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5));
        
        % *****************************************************************
        % E, exposed
        y(2) =  (beta0/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5))-(1/tau_L)*x(2);
        
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
        y(9)=h*(1-f)*(1/tau_IP)*x(3)-(1/tau_T)*x(9);
        
        y(10)=(1/tau_T)*x(9);
        
        y = y';
        
    end

% *************************************************************************

    function y = SEIR_equations1(t,x,N,beta0,ip,ia,is,f,g,h,tau_L,tau_IP,tau_I,tau_D)
        % Equations for SEIR model, with multiple infected classes

        ut=0;
        for jj=2:n

            Hp=0.5*(sign(t-ts_vec(jj-1))+1);
            Hm=0.5*(sign(ts_vec(jj)-t)  +1);

            ut=ut+u_control(jj)*Hp.*Hm;
        end

        Hm=0.5*(sign(ts_vec(1)-t)  +1);
        ut=u_control(1)*Hm+ut;

        beta_t=beta0*(1-ut);
        
        indicator=1;
        
        if((ip*x(3)+ia*x(4)+is*x(5))<1)
            indicator=0;
        end

        % *****************************************************************
        % S, susceptible
        y(1) = -indicator*(beta_t/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5));
        
        % *****************************************************************
        % E, exposed
        y(2) =  indicator*(beta_t/N)*x(1)*(ip*x(3)+ia*x(4)+is*x(5))-(1/tau_L)*x(2);
        
        % *****************************************************************
        % I, infectious
        
        % I_p, pre-symptomatic
        y(3) =  (1/tau_L)*x(2)-(1/tau_IP)*x(3);
        
        % I_a, asymptomatic, infectious
        y(4) =        f*(1/tau_IP)*x(3)-(1/tau_I)*x(4)   ;
        % I_s, sympomatic, infectious
        y(5) =    (1-f)*(1/tau_IP)*x(3)-(1/tau_I)*x(5)   ;
       
        % *****************************************************************
        % Removed / immune:
        
        y(6)=      (1/tau_I)*x(4)+((1-g)/tau_I)*x(5)    ;
        
        % *****************************************************************
        % Dying
        
        y(7)=(g/tau_I)*x(5)-(1/tau_D)*x(7);
        
        % *****************************************************************
        % Dead
        
        y(8)=(1/tau_D)*x(7);
        
        % *****************************************************************
        % Awaiting test
        y(9)=h*(1-f)*(1/tau_IP)*x(3)-(1/tau_T)*x(9);
        
        y(10)=(1/tau_T)*x(9);
        
        y = y';
        
    end

% *************************************************************************

end