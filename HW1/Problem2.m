clc;
clear all;
%% Problem 2
syms q1 q2 q3 a b c;
syms dq1 dq2 dq3;
syms m1 mt g;

% calculate frame transforms for position and jacobian
T01 = [cos(q1)  0       sin(q1) 0; ...
       sin(q1)  0       cos(q1) 0; ...
       0        -1      0       a; ...
       0        0       0       1];
   

T12 = [1        0       0       0; ...
       0        0       1       0; ...
       0        -1      0       b+q2; ...
       0        0       0       1];
   
T23 = [0        1       0       0; ...
       -1       0       0       0; ...
       0        0       1       c+q3; ...
       0        0       0       1];
   
T02 = T01*T12;
T03 = T02*T23

% position vector
P03 = T03(1:3,4)

% jacobian
jacob = [diff(P03(1:3),q1), diff(P03(1:3),q2), diff(P03(1:3),q3)]

% K and P for mass 1
K1 = 0;
P1 = m1*g*a;

% K and P for mass 2
Vmt = jacob * [dq1; dq2; dq3];
K2 = 0.5*mt*(Vmt.' * Vmt);
P2 = mt*g*P03(3);

K = K1+K2;
P = P1+P2;

% Lagrange 
L = K-P

syms Q1 Q2 Q3 Q1(t) Q2(t) Q3(t) ddq1 ddq2 ddq3
% Lagrange equation to calculate Tau 1
diffq1dot = diff(L,dq1);

diffq1dot = subs(diffq1dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq1dot = diff(diffq1dot, t);

diffq1dot = subs(diffq1dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq1 = diff(L,q1);

Tau1 = simplify(diffq1dot - diffq1)
% Lagrange equation to calculate F2
diffq2dot = diff(L,dq2);

diffq2dot = subs(diffq2dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq2dot = diff(diffq2dot, t);

diffq2dot = subs(diffq2dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq2 = diff(L,q2);

F2 = simplify(diffq2dot - diffq2)
% Lagrange equation to calculate F3
diffq3dot = diff(L,dq3);

diffq3dot = subs(diffq3dot, [q1 q2 q3 dq1 dq2 dq3], [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t)]);

diffq3dot = diff(diffq3dot, t);

diffq3dot = subs(diffq3dot, [Q1 Q2 Q3 diff(Q1(t),t) diff(Q2(t),t) diff(Q3(t),t) diff(Q1(t),t,t) diff(Q2(t),t,t) diff(Q3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

diffq3 = diff(L,q3);

F3= simplify(diffq3dot - diffq3)
