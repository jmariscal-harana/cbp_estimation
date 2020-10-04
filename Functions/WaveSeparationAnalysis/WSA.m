function [forward backward rhoCa] = WSA(P,U)
[rhoCa] = getRhoC(U, P);
pe = P;
ue = U;
dpe = diff(pe);
due = diff(ue);
dpep=(dpe + rhoCa * due)/2;        % forward pressure difference
dpem=(dpe - rhoCa * due)/2;        % backward pressure difference
forward=cumtrapz(dpep)+pe(1);         % forward pressure waveform
forward=[forward forward(end)];
forward = forward - min(forward);
backward=cumtrapz(dpem);               % backward pressure waveform
backward=[backward backward(end)];

function [rhoC] = getRhoC(flowPulse, pressPulse)

% dFlowPulse = diff(flowPulse, 1);
dFlowPulse = fsg(flowPulse, 1, 2, 5);
% dPressPulse = diff(pressPulse, 1);
dPressPulse = fsg(pressPulse, 1, 2, 5);
rhoC = sqrt(sum(dPressPulse .^ 2) / sum(dFlowPulse .^ 2));

function Z=fsg(Y,D,P,W)

[b,g]=sgolay(P,W);
HW=((W+1)/2)-1;
N=length(Y);
Z=zeros(1,N);
for n=HW+1:N-HW
    Z(n)=dot(g(:,D+1), Y(n-HW:n+HW));  
end

% fill in start and end values
if (D==0)
% use linear interpolation for D=0
    dZstart=(Z(HW+1)-Y(1))/HW;
    dZend=(Y(N)-Z(N-HW))/HW;
    for i=1:HW-1
        Z(i+1)=Y(1)+i*dZstart;
        Z(N-HW+i)=Z(N-HW)+i*dZend;
    end
else
% use constant value for derivatives
    Z(1:HW)=Z(HW+1)*ones(1,HW);
    Z(N-HW+1:N)=Z(N-HW)*ones(1,HW);
end
