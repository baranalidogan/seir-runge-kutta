%First, run modelcovidSEIR.m, then the rungekutta.m
function dydt = modelcovidSEIR(t,y)
w=391347.066;
beta1=2e-11;
beta2=2.2e-09;
beta=1/14;
mu=2.08547e-05;
lambda=2.35e-05;
too=0.03;
delta=3.4e-02;
S0=7610026000;
E0=80000;
I0=24545;
R0=907;
N=S0+E0+I0+R0;
S = y(1);
E = y(2);
I = y(3);
R = y(4);
dS_dt = w-(beta1*E+beta2*I)*S-mu*S;
dE_dt = (beta1*E+beta2*I)*S-(lambda+mu)*E;
dI_dt = lambda*E-(too+mu+delta)*I;
dR_dt = too*I-mu*R;
dydt = [dS_dt,dE_dt,dI_dt,dR_dt]; %[S, E, I, R];
end