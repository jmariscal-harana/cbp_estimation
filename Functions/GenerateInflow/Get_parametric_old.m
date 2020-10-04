function Qtot_old = Get_parametric_old(dt,Qi,Qmax,Qmin,Qi2,ti,tmax,tmin,ti2,ts,td,T)


%% Q2 = fourth order (for old wave); no horizontal slope at t=0

% 15 unknowns
%LHS_old
% [P14 P13 P12 P11 P10     P24 P23 P22 P21 P20    P34 P33 P32 P31 P30]
LHS_O = [
    tmax^4      tmax^3      tmax^2      tmax^1  1       0 0 0 0 0    0 0 0 0 0;  %Q1(tmax) = Qmax
    0           0           0           0       1       0 0 0 0 0    0 0 0 0 0;  %Q1(0) = 0
%    0           0           0           1       0      0 0 0 0 0    0 0 0 0 0;  %Q1'(0) = 0
    4*tmax^3    3*tmax^2    2*tmax^1    1       0       0 0 0 0 0    0 0 0 0 0;  %Q1'(tmax) = 0
    4*3*ti^2    3*2*ti      2           0       0       0 0 0 0 0    0 0 0 0 0;  %Q1''(ti) =0
    ti^4        ti^3        ti^2        ti      1       0 0 0 0 0    0 0 0 0 0;  %Q1(ti) = Qi
    
    0           0           0           0       0       tmax^4      tmax^3     tmax^2   tmax    1       0 0 0 0 0; %Q2(tmax)=Qmax
    0           0           0           0       0       4*tmax^3    3*tmax^2   2*tmax   1       0       0 0 0 0 0; %Q2'(tmax)=0
    0           0           0           0       0       ti2^4       ti2^3      ti2^2    ti2     1       0 0 0 0 0; %Q2(ti2) = Qi2
    0           0           0           0       0       4*3*ti2^2   3*2*ti2    2        0       0       0 0 0 0 0; %Q2'(tmax)=0
    
    0           0           0           0       0       ts^4        ts^3       ts^2     ts      1       -ts^4      -ts^3     -ts^2    -ts  -1; %Q2(ts)=Q3(ts)
    0           0           0           0       0       4*ts^3      3*ts^2     2*ts     1       0       -4*ts^3    -3*ts^2   -2*ts    -1   0;  %Q2'(ts)=Q3'(ts)  
    0 0 0 0 0   0 0 0 0 0     td^4        td^3        td^2    td      1; %Q3(td) = 0
    0 0 0 0 0   0 0 0 0 0     4*td^3      3*td^2      2*td    1       0; %Q3'(td)=0
    0 0 0 0 0   0 0 0 0 0     tmin^4      tmin^3      tmin^2  tmin    1; %Q3(tmin) = Qmin
    0 0 0 0 0   0 0 0 0 0     4*tmin^3    3*tmin^2    2*tmin  1       0  %Q3'(tmin) = 0
    ];

%RHS
RHS_O = [
    Qmax;
    0;
%     0; %Q1'(0) = 0
    0;
    0;
    Qi;
    Qmax;
    0;
    Qi2;
    0;
    0;
    0;
    0;
    0;
    Qmin;
    0
];

%Solve system A*X=B: x = A\B;
clear X;
X = LHS_O\RHS_O;


%% Equation of flow curve solution - OLD
% dt=1e-3;
%Q1, t1
t1=[0:dt:tmax];
Q1 = X(1).*t1.^4 + X(2).*t1.^3 + X(3).*t1.^2 +X(4).*t1 + X(5);

%Q2,t2
t2=[tmax+dt:dt:ts];
Q2 = X(6).*t2.^4 + X(7).*t2.^3 + X(8).*t2.^2 + X(9).*t2 + X(10);

%Q3,t3
t3 = [ts+dt:dt:td];
Q3 = X(11).*t3.^4 + X(12).*t3.^3 + X(13).*t3.^2 + X(14).*t3 + X(15);

%Q4, t4
t4 = [td+dt:dt:T];
Q4 = t4.*0;

Qtot_old = [Q1, Q2, Q3, Q4]; 