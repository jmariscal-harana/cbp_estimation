function [Pout, P0, tau, Est_hist_P] = Estimate_Pout_P0_and_tau(t, t_ini, t_end, P, Est_old_P)
%% Calculation of Pout and tau (for RT and CT) based on the exponential decay of P
% P = Pout + (P_0 - Pout)*exp(-t/tau)
t_exp = t(t_ini:t_end)-t(t_ini);
P_exp = P(t_ini:t_end);


Est_hist_P(1,1:3) = Est_old_P;
dt = t_exp(2,1)-t_exp(1,1);

% Parameter estimation
for it=1:15
    if it ~= 1
    Est_old_P = Est_new_P;
    end
    grad = zeros(length(Est_old_P),1);
    H = zeros(length(Est_old_P),length(Est_old_P));
    coste = 0.0;
    for k = 1:length(t_exp)
        p = P_exp(k);
        Po = Est_old_P(1);
        P0 = Est_old_P(2);
        tau = Est_old_P(3);
        tim = t_exp(k);
        y = p;
        f = Po + (P0 - Po)*exp(-tim/tau);
        dfdPo  = 1-exp(-tim/tau);
        dfdP0  = exp(-tim/tau); 
        dfdTau = tim*(P0 - Po)*tau^(-2)*exp(-tim/tau);

        dfdPoPo  = 0;                        dfdPoP0  = 0;                       dfdPoTau  = -tim/tau^2*exp(-tim/tau);
        dfdP0Po  = 0;                        dfdP0P0  = 0;                       dfdP0Tau  = tim/tau^2*exp(-tim/tau);
        dfdTauPo = -tim/tau^2*exp(-tim/tau); dfdTauP0 = tim/tau^2*exp(-tim/tau); dfdTauTau = (tim-2*tau)*tim*(P0 - Po)*tau^(-4)*exp(-tim/tau);

        grad(1,1) = grad(1,1) + (y-f)*dfdPo;
        grad(2,1) = grad(2,1) + (y-f)*dfdP0;
        grad(3,1) = grad(3,1) + (y-f)*dfdTau;

        H(1,1) = H(1,1) - dfdPo *dfdPo  + (y-f)*dfdPoPo  ;
        H(1,2) = H(1,2) - dfdP0 *dfdPo  + (y-f)*dfdPoP0  ;
        H(1,3) = H(1,3) - dfdTau*dfdPo  + (y-f)*dfdPoTau ;
        H(2,1) = H(2,1) - dfdPo *dfdP0  + (y-f)*dfdP0Po  ;
        H(2,2) = H(2,2) - dfdP0 *dfdP0  + (y-f)*dfdP0P0  ;
        H(2,3) = H(2,3) - dfdTau*dfdP0  + (y-f)*dfdP0Tau ;
        H(3,1) = H(3,1) - dfdPo *dfdTau + (y-f)*dfdTauPo ;
        H(3,2) = H(3,2) - dfdP0 *dfdTau + (y-f)*dfdTauP0 ;
        H(3,3) = H(3,3) - dfdTau*dfdTau + (y-f)*dfdTauTau;

        coste = coste + (y-f)^2;
    end
    grad = -2*grad;
    H = -2*H;
    b = inv(H)*grad;
    Est_new_P = Est_old_P - b;
    Est_hist_P(it+1,1:3) = Est_new_P;
    Coste_hist_P(it,1) = coste;
end
Pout = Est_hist_P(end,1); %(mmHg)
P0 = Est_hist_P(end,2); %(mmHg)
tau = Est_hist_P(end,3); %(s)
end