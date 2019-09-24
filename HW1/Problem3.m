clc;
clear all;
%% Problem 3
syms q1 q2 q3 a1 a2 a3;
syms dq1 dq2 dq3;
syms m1 m2 m3 g;

% calculate frame transforms for position and jacobian
T01 = [cos(q1)  -sin(q1) 0       a1*cos(q1); ...
       sin(q1)  cos(q1) 0       a1*sin(q1); ...
       0        0       1       0; ...
       0        0       0       1];
   

T12 = [cos(q2)  -sin(q2) 0       a1*cos(q2); ...
       sin(q2)  cos(q2) 0       a1*sin(q2); ...
       0        0       1       0; ...
       0        0       0       1];
   
T23 = [cos(q3)  -sin(q3) 0       a1 * cos(q3); ...
       sin(q3)  cos(q3) 0       a1 * sin(q3); ...
       0        0       1       0; ...
       0        0       0       1];
   
T02 = T01*T12;
T03 = T02*T23

% position vector
P01 = T01(1:3,4)
P02 = T02(1:3,4)
P03 = T03(1:3,4)

% Potential ENergies
P1 = m1 * g * subs(P01(3), a1, a1/2);
P2 = m2 * g * subs(P02(3), a2, a2/2);
P3 = m3 * g * subs(P03(3), a3, a3/2);


% Kinetic Energies
jacob1 = [diff(P01(1:3),q1)];
jacob2 = [diff(P02(1:3),q1), diff(P02(1:3),q2)];
jacob3 = [diff(P03(1:3),q1), diff(P03(1:3),q2), diff(P03(1:3),q3)];

Vm1 = subs(jacob1, a1, a1/2) * dq1;
Vm2 = subs(jacob2, a2, a2/2) * [dq1; dq2];
Vm3 = subs(jacob3, a3, a3/2) * [dq1; dq2; dq3];

K1 = 0.5*m1*(Vm1.' * Vm1);
K2 = 0.5*m2*(Vm2.' * Vm2);
K3 = 0.5*m3*(Vm3.' * Vm3);
% Lagrange 
P = P1+P2+P3;
K = K1+K2+K3;
L = K-P

syms Q1 Q2 Q3 Q1(t) Q2(t) Q3(t) ddq1 ddq2 ddq3
% Lagrange equation to calculate Tau 1
diffq1dot = diff(L,dq1);

diffq1dot = subs(diffq1dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq1dot = diff(diffq1dot, t);

diffq1dot = subs(diffq1dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq1 = diff(L,q1);

Tau1 = simplify(diffq1dot - diffq1)
% Lagrange equation to calculate Tau2
diffq2dot = diff(L,dq2);

diffq2dot = subs(diffq2dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq2dot = diff(diffq2dot, t);

diffq2dot = subs(diffq2dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq2 = diff(L,q2);

Tau2 = simplify(diffq2dot - diffq2)
% Lagrange equation to calculate Tau3
diffq3dot = diff(L,dq3);

diffq3dot = subs(diffq3dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq3dot = diff(diffq3dot, t);

diffq3dot = subs(diffq3dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq3 = diff(L,q3);

Tau3= simplify(diffq3dot - diffq3)