% script: continuous_pH_model

S = 60;
b1 = 20*(rand(S,1)-1/2); % 3; 
b2 = 1; % b1/3; 
pc = 0.6; sig2 = 1;
K = 1;
al = 5;
fl = 1; % stabilize p
p_opt_vec = 1.4*ones(S,1) + 1*randn(S,1);
lam = 1e-5;

Gamma_p = @(p,p_opt,pc,sig2) ...
   1 - (1-exp(-(p-p_opt).^2/2./sig2))./( 1-exp(-pc^2/2./sig2) );
dn_dt_p = @(n,p,al,K,p_opt_vec,pc,sig2,lam) ...
  al*n.*(1 - n/K + Gamma_p(p,p_opt_vec,pc,sig2)) + lam;
dp_dt = @(n,p,b1,b2,al,K,p_opt_vec,pc,sig2,fl) ...
  al/S*sum( b1.*dn_dt_p(n,p,al,K,p_opt_vec,pc,sig2,lam).*...
   (dn_dt_p(n,p,al,K,p_opt_vec,pc,sig2,lam)>0) ) + al/S*b2*sum(n) - fl*(p-1);

dnp_dt = @(t,np) ...
   [(dn_dt_p(np(1:S),np(S+1),al,K,p_opt_vec,pc,sig2,lam));...
   dp_dt(np(1:S),np(S+1),b1,b2,al,K,p_opt_vec,pc,sig2,fl)];
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

n0 = rand(S,1); n0 = n0/sum(n0) * K/2;
np0 = [n0 ; 1];
[T_now,NPt_now] = ode45(dnp_dt,0:0.1:240,np0,options);
plot(T_now,NPt_now(:,1:S)), hold on
plot(T_now,NPt_now(:,end),'k','LineWidth',1), hold off