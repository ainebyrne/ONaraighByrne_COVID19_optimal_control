function y = calculate_R0(u)

% Maximum allowed value of R0:
% https://academic.oup.com/jtm/article/27/2/taaa021/5735319

V = zeros(4,4);
F = zeros(4,4);

beta0=u(1);
f=u(2);
g=u(3);
h=1;
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

% *************************************************************************
% Calculate next-generation matrix G.

F(1,2)=beta0*ip;
F(1,3)=beta0*ia;
F(1,4)=beta0*is;

V(1,1)=1/tau_L;

V(2,1)=-1/tau_L;
V(2,2)=1/tau_IP;

V(3,2)=-f/tau_IP;
V(3,3)=1/tau_I;

V(4,3)=-(1-f)/tau_IP;
V(4,4)=1/tau_I;

G=F*inv(V);

% *************************************************************************


y = max(eigs(G));

end

