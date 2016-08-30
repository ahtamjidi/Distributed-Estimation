function [A,x0,B,C] = create_sys_gold()
%Source and receptors are not on the boundaries (except z = 0)
%10^5 system, Full State Measurements
global opt_dist
%% Simple 2*2 system
% % % % source.Q = [1 1];
% % % % A = [1 0;0 1];
% % % % %Receptors
% % % % 
% % % % B = [1 0;0 1];
% % % % % C = randn(9,2);%repmat(eye(2),4,1);[1 0]];
% % % % C = [repmat(eye(2),4,1);[1 0]];

% C = zeros(9,2);

%% 9*9 system
A = eye(9);
B = eye(9);
source.Q = ones(1,9);

% C = randn(9,2);%repmat(eye(2),4,1);[1 0]];
C = eye(9);


% x0 =zeros(size(A,1),1);
x0 = [1:9]';
opt_dist.order = 9;
opt_dist.A = A;
opt_dist.B = B;
opt_dist.C = C;
opt_dist.source = source;



