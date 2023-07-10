% script: daily_dilutions_pH_model

% couple species through pH:
Gamma_p = @(p,p_opt,pc,sig2) ...
   1 - (1-exp(-(p-p_opt).^2/2./sig2))./( 1-exp(-pc^2/2./sig2) );
dn_dt_p = @(n,p,al,K,p_opt_vec,pc,sig2) al*n.*(1-sum(n)/K).*Gamma_p(p,p_opt_vec,pc,sig2);
dp_dt = @(n,p,b1,b2,al,K,p_opt_vec,pc,sig2) al*b1*sum( dn_dt_p(n,p,al,K,p_opt_vec,pc,sig2).*...
   (dn_dt_p(n,p,al,K,p_opt_vec,pc,sig2)>0) ) + al*b2*sum(n);

S = 60; 
b1 = 3; 
b2 = b1/3; pc = 0.6; sig2 = 1;
K = 1; % total capacity
al = 1;
p_opt_vec = 1.4*ones(S,1)+0.3*randn(S,1);

dnp_dt = @(t,np) ...
   [(dn_dt_p(np(1:S),np(S+1),al,K,p_opt_vec,pc,sig2).*(np(1:S)>1e-16));...
   dp_dt(np(1:S),np(S+1),b1,b2,al,K,p_opt_vec,pc,sig2)];
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% daily dilution:
dilution = 100;
days = 200;
T0 = 0; np0 = [1e-5*rand(S,1);1];
T = []; NPt = [];
for kk=1:days
   [T_now,NPt_now] = ode45(dnp_dt,T0+(0:0.1:24),np0,options);
   if T_now < T0+24, break, end % failed to complete integration
   T = [T ; T_now]; NPt = [NPt ; NPt_now]; %#ok<AGROW>
   np0(1:S) = NPt_now(end,1:S)/dilution; np0(S+1) = 1; % dilution sets p=1(?)
   np0(1:S) = np0(1:S).*(np0(1:S)>1e-16);
   T0 = T_now(end);
end
inds = 240:241:length(T);

subplot(2,1,1)
plot(T(inds)/24,NPt(inds,1:S),'-')
set(gca,'YScale','log')
xlabel('t'), ylabel('N_i')
subplot(2,1,2)
plot(T(inds)/24,1-log(NPt(inds,S+1)))
xlabel('t'), ylabel('pH')