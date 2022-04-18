%%%%%
%%% 5SMB0 System Identification
%%% Exercise 3 Problem 3
%%% Author: Jiaxuan Zhang
%%%%%

%% Question 1

% answer: Yes, it holds

%% Question 2

% answer: Yes, first, \epsilon meet convergence request, 
%         then if we choose correct model set, with input which has sufficient order, 
%         we can then meet consistency request, that is, the conclustion hold

%% Question 3
% answer: yes, S \in M
load('dataG0arx.mat')
data=iddata(y,u);
ARX.sys=arx(data,[4,2,3]);
ARX.resid=resid(data,ARX.sys);

figure
resid(data,ARX.sys)
grid on

%% Question 4
% answer: nop, because S \in M does not hold, so we can only have G=G0 but
%         we cannot have H=H0

%% Question 5
% answer: the result shows G is appropriate but H is not appropriate
OE.sys=oe(data,[2,4,3]);
OE.resid=resid(data,OE.sys);

figure
resid(data,OE.sys)
grid on;

%% Question 6
% answer: yes, confirm my previous conclusions



%% Question 7
% answer: the difference is caused by the difference of the process. In
%         Exercise 3, the system noise has different property. 


