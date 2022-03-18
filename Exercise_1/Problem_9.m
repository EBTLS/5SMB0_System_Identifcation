%%%%%
%%% 5SMB0 System Identification
%%% Exercise 1 Problem 9
%%% Author: Jiaxuan Zhang
%%%%%
%% load data
load('batch.mat')
figure
plot(y,'r')
hold on
plot(u,'b')
grid on


%% non-smoothed etfe
figure
G1 = etfe(iddata(y, u), [], 512)
bode(G1)

%% Judge window size
% Question: why we judge 30 as the turning point? why the increment of Ryu is weird?
figure
cra(iddata(y, u),500);
title('auto correlation')
grid on;

%% smoothed etfe
figure
G2 = etfe(iddata(y, u), 30, 512)
bode(G2)
