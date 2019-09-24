clc;
clear all;
%% Problem 2
syms q1 q2 q3 a1 a2 a3;
syms dq1 dq2 dq3;
syms m1 m2 m3 g;

% calculate frame transforms for position and jacobian
T01 = [cos(q1)  -sin(q1)0       a1*cos(q1); ...
       sin(q1)  cos(q1) 0       a1*sin(q1); ...
       0        0       1       0; ...
       0        0       0       1];