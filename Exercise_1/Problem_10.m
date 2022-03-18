%%%%%
%%% 5SMB0 System Identification
%%% Exercise 1 Problem 10
%%% Author: Jiaxuan Zhang
%%%%%

%% Generate Input-Output Data
G0 = tf([1 - 0.5 0.2 0.8 0.1], [1 -1.7 1.6 -0.8 0.25], 1);
H0 = tf([1], [1, -1.7, 1.6, -0.8, 0.25], 1);
N = 1024;
u = sign(randn(N, 1));
lambda = sqrt(0.1);
e = lambda * randn(N, 1);
y = lsim(G0, u) + lsim(H0, e);

%% plot signals
figure
plot(u)
hold on
plot(e)
hold on
plot(y)
legend('u', 'e', 'y')

%% non-smoothed etfe
% figure
G1 = etfe(iddata(y, u), [], 512)
% bode(G1)

%% judge smoothing window size
figure
cra(iddata(y, u), 500);
title('auto correlation')
grid on;

%% smoothed etfe
figure
G2 = etfe(iddata(y, u), 30, 512)
G3 = etfe(iddata(y, u), 50, 512)
G4 = etfe(iddata(y, u), 10, 512)
G5 = etfe(iddata(y, u), 200, 512)
bode(G1)
hold on
bode(G2)
bode(G3)
bode(G4)
bode(G5)
legend('0','30','50','10','250')
