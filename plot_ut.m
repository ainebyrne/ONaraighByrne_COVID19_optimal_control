function ut_vec=plot_ut(u,t_vec)

% *************************************************************************
% Controls:
tfinal=t_vec(end);

nn=length(u);
n=nn/2;

u_control=u(1:n);

n_hat=u(n+1:end)/norm(u(n+1:end));
xi_s=n_hat.^2;
tau_vec=tfinal*xi_s;

ts_vec=0*tau_vec;
for ii=1:n
    ts_vec(ii)=sum(tau_vec(1:ii));
end

u_vec=0*t_vec;

for jj=2:n

    Hp=0.5*(sign(t_vec-ts_vec(jj-1))+1);
    Hm=0.5*(sign(ts_vec(jj)-t_vec)  +1);

    u_vec=u_vec+u_control(jj)*Hp.*Hm;
end

Hm=0.5*(sign(ts_vec(1)-t_vec)  +1);
ut_vec=u_control(1)*Hm+u_vec;

% *************************************************************************


end