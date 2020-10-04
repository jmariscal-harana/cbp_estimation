% P = REF.Pressure;
% Q = REF.Q_in;
% t = REF.t_in;
%- Modulus
Pfft = fft(P);
Pz = abs(Pfft(1:10));
Qfft = fft(Q);
Qz = abs(Qfft(1:10));
Modulus = Pz./Qz;
subplot(2,1,1)
plot(Modulus,'-ob','MarkerFaceColor','b')
axis([0 10 0 2])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'})
set(gca,'fontsize',14)
ylabel('Modulus Z [mmHg s/mL]')
box off
%- Phase
Pp = angle(Pfft(1:10));
Qp = angle(Qfft(1:10));
Phase = Pp - Qp;
subplot(2,1,2)
plot(Phase,'-ob','MarkerFaceColor','b')
axis([0 10 -10 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'})
line([0 10],[0 0],'linestyle','- -','color','k')
set(gca,'fontsize',14)
xlabel('Frequency [Hz]','fontsize',14)
ylabel('Phase \Phi [rad]')
box off