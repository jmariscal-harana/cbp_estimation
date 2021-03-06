function [R1, R2, C, Est_hist, FinalCost, fig_num] = Estimate_R2_and_C(t, P_raw, Q_raw, Pout, Est_old, Lim_t, Lim, plots, fig_num)
% This function estimates parameters R1 and C, given that Pout is known and
% (Pavg - Pout)/Qavg = Rt = R1 + R2 -> R2 = Rt - R1
% P = (Po - R1想o - Pout)搪xp((t-to)/(R2C)) + R1想 + Pout +
% exp(-t/(R2嵩))/C*int(exp(t'/(R2嵩))想搞t',to,t)

%% First step - Filtering P and Q functions
num = 15; %number of harmonics considered
P_F = FourierSeries(t, P_raw, num);
Q_F = FourierSeries(t, Q_raw, num);
P = zeros(length(t),1);
Q = zeros(length(t),1);
T = t(end);
w = 2*pi/T;
for i = 1:length(t)
        P(i,1) = P_F(1,1)/2; Q(i,1) = Q_F(1,1)/2;
        for v = 1:length(P_F)-1
            P(i,1) = P(i,1) + P_F(v+1,1)*cos(v*w*t(i,1)) + P_F(v+1,2)*sin(v*w*t(i,1));
            Q(i,1) = Q(i,1) + Q_F(v+1,1)*cos(v*w*t(i,1)) + Q_F(v+1,2)*sin(v*w*t(i,1));
        end
end

%% Second step - Estimation
Rt = (mean(P)-Pout)/mean(Q);
Est_hist(1,1:2) = Est_old;
t0 = t(1,1);
p0 = P(1,1); q0 = Q(1,1);
% Cost function
% Calculation
eval(['X_1 = linspace(' num2str(Lim_t(1,1)) ',' num2str(Lim_t(1,2)) ',50);']);
eval(['Y_1 = linspace(' num2str(Lim_t(2,1)) ',' num2str(Lim_t(2,2)) ',50);']);
Z_1 = Calculate_Cost_Matrix(t,X_1,Y_1,P,Q,Pout);
eval(['X_2 = linspace(' num2str(Lim(1,1)) ',' num2str(Lim(1,2)) ',50);']);
eval(['Y_2 = linspace(' num2str(Lim(2,1)) ',' num2str(Lim(2,2)) ',50);']);
Z_2 = Calculate_Cost_Matrix(t,X_2,Y_2,P,Q,Pout);

if plots == 1
	linewidth = 1.5;
	markersize = 10;
	linewidth = 1.5;
	markersize = 10;
	fig_num  = fig_num + 1;
	eval(['figure(' num2str(fig_num) ')'])
	eval('surf(X_1,Y_1,Z_1)')
	xlabel('R_2 (mmHg新/mL)')
	ylabel('C (mL/mmHg)')
	title(['Cost function'])
	fig_num  = fig_num + 1;
	levels = 10;
	eval(['figure(' num2str(fig_num) ')'])
	eval(['contourf(X_2,Y_2,Z_2,' num2str(levels) ')'])
	xlabel('R_2 (mmHg新/mL)')
	ylabel('C (mL/mmHg)')
	title(['Cost function'])
end

% Parameter estimation
for it=1:5
    if it ~= 1
       Est_old = Est_new;
    end
    grad = zeros(length(Est_old),1);
    H = zeros(length(Est_old),length(Est_old));
    coste = 0.0;
    R2 = Est_old(1);
    R1 = Rt - R2;
    C = Est_old(2);
    for k = 1:length(t)
        e1 = exp(-(t(k)-t0)/(R2*C));
        e2 = exp(-t(k)/(R2*C));
        p = P(k,1); q = Q(k,1);
        AA = zeros(k,1); BB = zeros(k,1); CC = zeros(k,1);
        aa = 0; bb = 0; cc = 0;
        for m = 1:k
            AA(m,1) = Q(m)*exp(t(m)/(R2*C));
            BB(m,1) = t(m)*Q(m)*exp(t(m)/(R2*C));
            CC(m,1) = t(m)^2*Q(m)*exp(t(m)/(R2*C));
        end
        if k ~= 1
            aa = trapz(t(1:k),AA);
            bb = trapz(t(1:k),BB);
            cc = trapz(t(1:k),CC);
        end
        y     = p;
        f     = (p0 - R1*q0 -Pout)*e1 + R1*q + Pout + e2/C*aa;
        dfdR  = q0*e1 + (p0 - R1*q0 -Pout)*(t(k)-t0)/(R2^2*C)*e1 - q +  t(k)/(R2^2*C^2)*e2*aa - e2/(R2^2*C^2)*bb;
        dfdC  = (p0 - R1*q0 -Pout)*(t(k)-t0)/(R2*C^2)*e1 + (t(k)-R2*C)/(R2*C^3)*e2*aa - e2/(R2*C^3)*bb;
        dfdRR = 2*q0*(t(k)-t0)/(R2^2*C)*e1 + (p0 - R1*q0 -Pout)*((t(k)-t0)^2-2*(t(k)-t0)*R2*C)/(R2^4*C^2)*e1 + (t(k)^2-2*t(k)*R2*C)/(R2^4*C^3)*e2*aa - (2*t(k)+2*R2*C)/(R2^4*C^3)*e2*bb + e2/(R2^4*C^3)*cc;
        dfdRC = (t(k)-t0)*q0/(R2*C^2)*e1 +(p0 - R1*q0 -Pout)*((t(k)-t0)^2-(t(k)-t0)*R2*C)/(R2^3*C^3)*e1 + (t(k)*(t(k)-R2*C)-t(k)*R2*C)/(R2^3*C^4)*e2*aa + 2*(R2*C-t(k))/(R2^3*C^4)*e2*bb - 1/(R2^3*C^4)*e2*cc;
        dfdCC = (p0 - R1*q0 -Pout)*((t(k)-t0)^2-2*(t(k)-t0)*R2*C)/(R2^2*C^4)*e1 + (t(k)*(t(k)-R2*C)-3*(t(k)-R2*C)*R2*C-R2^2*C^2)/(R2^2*C^5)*e2*aa + (4*R2*C-2*t(k))/(R2^2*C^5)*e2*bb + 1/(R2^2*C^5)*e2*cc;
        
        grad(1,1) = grad(1,1) + (y-f)*dfdR;
        grad(2,1) = grad(2,1) + (y-f)*dfdC;

        H(1,1) = H(1,1) - dfdR*dfdR + (y-f)*dfdRR;
        H(1,2) = H(1,2) - dfdC*dfdR + (y-f)*dfdRC;
        H(2,1) = H(2,1) - dfdR*dfdC + (y-f)*dfdRC;
        H(2,2) = H(2,2) - dfdC*dfdC + (y-f)*dfdCC;

        coste = coste + (y-f)^2;
    end
    grad = -2*grad;
    H = -2*H;
    b = H\grad;
    %----------------Esto lo tuve que a鎙dir---------------------------
%     if Est_old(1) == 0 || Est_old(1) == Rt
%         b(1) = 0;
%     elseif Est_old(2) == 0
%         b(2) = 0;
%     end
    %------------------------------------------------------------------
    Est_new = Est_old - b;
    %---------------Esto lo tuve que a鎙dir----------------------------
%     if (Est_new(1) < 0 || Est_new(2) < 0 || Est_new(1) > Rt)
%        m = (Est_new(2)-Est_old(2))/(Est_new(1)-Est_old(1));
%        n = (Est_new(2)*(Est_new(1)-Est_old(1))-(Est_new(2)-Est_old(2))*Est_new(1))/(Est_new(1)-Est_old(1));
%        if Est_new(1) < 0
%            if n > 0
%                 Est_new(1) = 0;
%                 Est_new(2) = n;
%            else
%                 Est_new(2) = 0;
%                 Est_new(1) = -n/m;
%            end
%        elseif Est_new(2) < 0
%            Est_new(2) = 0;
%            Est_new(1) = -n/m;
%        elseif Est_new(1) > Rt
%            if m*Rt+n > 0
%                 Est_new(1) = Rt;
%                 Est_new(2) = m*Rt+n;
%            else
%                 Est_new(2) = 0;
%                 Est_new(1) = -n/m;
%            end
%        end
%     end
    %------------------------------------------------------------------
    Est_hist(it+1,1:2) = Est_new;
    Coste_hist(it,1) = coste;
end
R2 = Est_new(1);
R1 = Rt - Est_new(1);
C = Est_new(2);
FinalCost = Coste_hist(end,1);

if plots == 1
	eval(['figure(' num2str(fig_num) ')'])
	eval('hold on')
	eval(['plot(Est_hist(:,1), Est_hist(:,2), '':rs'', ''LineWidth'', ' num2str(linewidth) ', ''MarkerSize'', ' num2str(markersize) ')'])
	
	fig_num  = fig_num + 1;
	eval(['figure(' num2str(fig_num) ')'])
	subplot(3,1,1)
	eval(['plot(Est_hist(1:end-1,1), ''Color'', ''r'', ''LineStyle'', ''-'', ''LineWidth'', ' num2str(linewidth) ')']) %(mmHg)
	hold on
	xlabel('Iteration')
	ylabel('R_2 (mmHg新/mL)')
	title(['R_2 at each iteration'])
	subplot(3,1,2)
	eval(['plot(Est_hist(1:end-1,2), ''Color'', ''r'', ''LineStyle'', ''-'', ''LineWidth'', ' num2str(linewidth) ')']) %(mmHg)
	xlabel('Iteration')
	ylabel('C (mL/mmHg)')
	title(['C at each iteration'])
	subplot(3,1,3)
	eval(['plot(Coste_hist(1:end,1), ''Color'', ''r'', ''LineStyle'', ''-'', ''LineWidth'', ' num2str(linewidth) ')']) %(mmHg)
	xlabel('Iteration')
	ylabel('Cost (mmHg^{2})')
	title(['Cost at each iteration'])
end

end
