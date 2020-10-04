% The same as Plot_Estimated_Windkessel, but this one is for depicting the
% measured P and the estimated P.
function [mean_P, mean_Pest, fig_num] = Plot_Estimated_WindkesselB(t, P, Q, R1, R2, C, Pout, fig_num)
P0 = P(1);
t0 = t(1);
T  = t(end);
dt = t(2)-t(1);
Q0 = Q(1);
dt = t(2)-t(1);

tt(1) = t0;
for i = 2:12*length(t)
   tt(i) = tt(i-1) + dt; 
end
PP = [P; P; P; P; P; P; P; P; P; P; P; P];
QQ = [Q; Q; Q; Q; Q; Q; Q; Q; Q; Q; Q; Q];

for i = 1:length(tt)
   %Term 1
   term1 = (P0 - R1*Q0 - Pout)*exp(-(tt(i)-t0)/(R2*C));
   TERM1(i,1) = term1;
   %Term 2
   term2 = R1*QQ(i);
   TERM2(i,1) = term2;
   %Term 3
   term3 = Pout;
   TERM3(i,1) = term3;
   %Term 4
   AA = zeros(i,1);
   aa = 0;
   for j = 1:i
       AA(j,1) = QQ(j)*exp(tt(j)/(R2*C)); 
   end
   if i~=1
      aa = trapz(tt(1:i),AA); 
   end
   term4 = exp(-tt(i)/(R2*C))/C*aa;
   TERM4(i,1) = term4;
   %Sum
   P_est(i,1) = term1 + term2 + term3 + term4;
end

mean_P = mean(PP);
mean_Pest = mean(P_est);

%Plot
fig_num  = fig_num + 1;
eval(['figure(' num2str(fig_num) ')'])
plot(tt(1:1*T/dt+1),PP(8*T/dt:9*T/dt),tt(1:1*T/dt+1),P_est(8*T/dt:9*T/dt), 'linewidth', 2.5)
legend('Measured','Simulated')
xlabel('t (s)')
ylabel('P (mmHg)')

% fig_num  = fig_num + 1;
% eval(['figure(' num2str(fig_num) ')'])
% plot(tt,TERM1, tt, TERM2, tt, TERM3, tt, TERM4)
% legend('Initial value','R_{1}·Q', 'P_{out}', 'Integral')
end