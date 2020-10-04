function Z = Calculate_Cost_Matrix(t,X,Y,P,Q,Pout)
    t = t + 14*t(end);
    t0 = t(1,1);
    p0 = P(1,1);
    q0 = Q(1,1);
    Rt = (mean(P)-Pout)/mean(Q);
    for i = 1:length(X)
        for j = 1:length(Y)
            coste = 0.0;
            R2 = X(i);
            R1 = Rt - R2;
            C = Y(j);
            for k = 1:length(t)
                e1 = exp(-(t(k)-t0)/(R2*C));
                e2 = exp(-t(k)/(R2*C));
                p = P(k,1); q = Q(k,1);
                AA = zeros(k,1); %BB = zeros(k,1); CC = zeros(k,1);
                aa = 0; %bb = 0; cc = 0;
                for m = 1:k
                    AA(m,1) = Q(m)*exp(t(m)/(R2*C));
                    %BB(m,1) = t(m)*Q(m)*exp(t(m)/(R2*C));
                    %CC(m,1) = t(m)^2*Q(m)*exp(t(m)/(R2*C));
                end
                if k ~= 1
                    aa = trapz(t(1:k),AA);
                    %bb = trapz(t(1:k),BB);
                    %cc = trapz(t(1:k),CC);
                end
                y     = p;
                f     = (p0 - R1*q0 -Pout)*e1 + R1*q + Pout + e2/C*aa;
                %dfdR  = q0*e1 + (p0 - R1*q0 -Pout)*(t(k)-t0)/(R2^2*C)*e1 - q +  t(k)/(R2^2*C^2)*e2*aa - e2/(R2^2*C^2)*bb;
                %dfdC  = (p0 - R1*q0 -Pout)*(t(k)-t0)/(R2*C^2)*e1 + (t(k)-R2*C)/(R2*C^3)*e2*aa - e2/(R2*C^3)*bb;
                %dfdRR = 2*q0*(t(k)-t0)/(R2^2*C)*e1 + (p0 - R1*q0 -Pout)*((t(k)-t0)^2-2*(t(k)-t0)*R2*C)/(R2^4*C^2)*e1 + (t(k)^2-2*t(k)*R2*C)/(R2^4*C^3)*e2*aa - (2*t(k)+2*R2*C)/(R2^4*C^3)*e2*bb + e2/(R2^4*C^3)*cc;
                %dfdRC = (t(k)-t0)*q0/(R2*C^2)*e1 +(p0 - R1*q0 -Pout)*((t(k)-t0)^2-(t(k)-t0)*R2*C)/(R2^3*C^3)*e1 + (t(k)*(t(k)-R2*C)-t(k)*R2*C)/(R2^3*C^4)*e2*aa + 2*(R2*C-t(k))/(R2^3*C^4)*e2*bb - 1/(R2^3*C^4)*e2*cc;
                %dfdCC = (p0 - R1*q0 -Pout)*((t(k)-t0)^2-2*(t(k)-t0)*R2*C)/(R2^2*C^4)*e1 + (t(k)*(t(k)-R2*C)-3*(t(k)-R2*C)*R2*C-R2^2*C^2)/(R2^2*C^5)*e2*aa + (4*R2*C-2*t(k))/(R2^2*C^5)*e2*bb + 1/(R2^2*C^5)*e2*cc;

                %grad(1,1) = grad(1,1) + (y-f)*dfdR;
                %grad(2,1) = grad(2,1) + (y-f)*dfdC;

                %H(1,1) = H(1,1) - dfdR*dfdR + (y-f)*dfdRR;
                %H(1,2) = H(1,2) - dfdC*dfdR + (y-f)*dfdRC;
                %H(2,1) = H(2,1) - dfdR*dfdC + (y-f)*dfdRC;
                %H(2,2) = H(2,2) - dfdC*dfdC + (y-f)*dfdCC;

                coste = coste + (y-f)^2;
            end
            Z(j,i) = coste;
        end
    end
end