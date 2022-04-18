%%%%%
%%% 5SMB0 System Identification
%%% Exercise 3 Problem 4
%%% Author: Jiaxuan Zhang
%%%%%
%% Data Generation
G0 = tf([1 - 0.5 0.2 0.8 0.1], [1 -1.7 1.6 -0.8 0.25], 1);
H0 = tf([1], [1, -1.7, 1.6, -0.8, 0.25], 1);
N = 500;
u = sign(randn(N, 1));
lambda = sqrt(0.1);
e = lambda * randn(N, 1);
y = lsim(G0, u) + lsim(H0, e);
data=iddata(y,u);
%% Question 1: ARX model
Q1.sys=arx(data,[4,5,0],arxOptions('Focus','prediction'));
Q1.res=resid(data,Q1.sys);

figure
resid(data,Q1.sys)
grid on

%% Question 2: FIR model
Q2.sys=arx(data,[0,10,0],arxOptions('Focus','prediction'))
Q2.res=resid(data,Q2.sys);


figure
subplot(1,2,1)
resid(data,Q2.sys)
subplot(1,2,2)
bode(Q2.sys)
grid on

%% Question 3: 
% answer: not all dynamics are captured by the FIR model, i.e., 
% there is still information about the dynamics present in the residual signal

%% Quesiont 4:
% answer: no, not an element
Q4.arx_sys=arx(data,[3,3,0],arxOptions('Focus','prediction'));
Q4.arx_res=resid(data,Q4.arx_sys);
Q4.oe_sys=oe(data,[3,3,0],arxOptions('Focus','prediction'));
Q4.oe_res=resid(data,Q4.oe_sys);

figure
subplot(1,2,1)
resid(data,Q4.arx_sys)
subplot(1,2,2)
resid(data,Q4.oe_sys)

%% Question 5:
% answer: here because our order chosen cannot guarantee s \in m, so here
% we need use approximate modelling to analyze it.
% Basaed on approximate modelling, in ARX model, becasue of H \neq 1, the e
% part has a different influence on the V (high-frequency misfit penalization in ARX identification)

figure 
bode(Q4.arx_sys,'r')
hold on
bode(Q4.oe_sys,'b')

%% Question 6:
% answer: yes, we can add low pass filter L, so that the low frequency part
% will get higher weight and smaller error

fdata = idfilt([y u], [0,1]);
Q6.sys = arx(fdata,[3,3,0],arxOptions('Focus','prediction'));

% original_sys=tf2ss([1 - 0.5 0.2 0.8 0.1], [1 -1.7 1.6 -0.8 0.25]);

figure 
bode(Q4.arx_sys)
hold on
bode(Q4.oe_sys)
hold on
bode(Q6.sys)
hold on
bode(G0)
legend('arx','oe','filtered arx','original')

