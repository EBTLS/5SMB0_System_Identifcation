%%%%%
%%% 5SMB0 System Identification
%%% Exercise 2 Problem 4
%%% Author: Jiaxuan Zhang
%%%%%
%% Question 1
% Answer: 
% The minimal order of excitation of u(t) is 2+4=6

%% Question 2
% Answer:
% We need to chose nb, nf, nk

%% Question 3
% Answer:
% Yes, we can because it is has higher order

%% Question 4
% load data
load('dataG0oe.mat')

figure
M1=oe(iddata(y,u),[2,4,3]);
present(M1);

% residue test
resid(iddata(y,u),M1);
grid on;

title('oe model')

% Answer: yes

%% Question 5
% Answer:
% No, because G0(q) and H0(q) are conflicted in ARX model

%% Question 6

figure
M2=arx(iddata(y,u),[4,2,3]);
present(M2);

% residue test
resid(iddata(y,u),M2);
grid on;

title('arx model')

% Answer: Yes

%% Question 7
% Answer:
% Yes

