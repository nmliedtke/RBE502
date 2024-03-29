%% Question 2
%% Establish symbolic variables
syms a b c d
%% Part A
A1 = [   -2  1   0   0   0;...
        0   -1  1   0   0;...
        0   0   -2  1   0;...
        0   0   0   -2  1;...
        0   0   0   0   -1]
B1 = [0 0; 0 0; a b ; c 0; 0 d]

testM = [B1 A1*B1 (A1^2)*B1 (A1^3)*B1 (A1^4)*B1]
rank1 = rank(testM)
% In order for A1 and B1 to allow reachability, the matrix testM must be Full Row Rank
% In this case testM has 5 rows, so rank(testM) must return 5
% If a, b, c, d are all diferent values then the pair will always be
% reachable. However if some of their values are identical they
% pair may or may not lead to reachability. 
% If b and d or c and d equal zero then testM loses rank and becomes
% not reachable.
%% Part B
A2 = [   -1  1   0   0   0;...
        0   0  1   0   0;...
        0   0   -1  0   0;...
        0   0   0   -1  1;...
        0   0   0   0   0];
B2 = [0 0; 0 0;a b; 0 0; c d];

testM2 = [B2 A2*B2 (A2^2)*B2 (A2^3)*B2 (A2^4)*B2 ]
rank2 = rank(testM2)

% In order for A2 and B2 to allow reachability, the matrix testM must be Full Row Rank
% In this case testM2 has 5 rows, so rank(testM2) must return 5
% If a, b, c, d are all diferent values then the pair will always be
% reachable. However if some of their values are identical they
% pair may or may not lead to reachability. 
% If a=c and b=d pair is not reachable. If a and b or a and c or b and d or c and d equal zero then testM loses rank and becomes
% not reachable.