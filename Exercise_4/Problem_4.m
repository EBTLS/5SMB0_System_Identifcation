%%%%%
%%% 5SMB0 System Identification
%%% Exercise 4 Problem 4
%%% Author: Jiaxuan Zhang
%%%%%
%% load data
load('batch.mat')
figure
plot(y,'r')
hold on
plot(u,'b')
grid on

data=iddata(y,u);


%% Question 1: Select Full Order Model Structure

% NN=struc(1:20,1:20,1);
% v=arxstruc(data(1:1:floor(length(y)/2)),data(floor(length(y)/2)+1:1:end),NN);
% selstruc(v) % choose 4 2 1
% % order=[4,2,1];
% 
% sys=arx(data,[4,2,1]);
% resid(data,sys)

% use OE model
OE.sys=oe(data,[5,5,1]);
figure
resid(data,OE.sys) % G is perfect but H is not

% use BJ model
BJ.sys=bj(data,[5,5,5,5,1])
figure
resid(data,BJ.sys) % all are perfect

%% Question 2
BJ.sys

%% Question 3
% w=logspace(1e-3,10);
% [covG,G]=compcovGhat(BJ.sys,w);
% plot(w,covG)
% use systemIdentification Toolbox

%% Question 4: from solution
M = bj([y u],[5 5 5 5 1]); 
[A,B,C,D,F] = polydata(M);  % Retrieve polynomial coeff
H = tf(C,D,1); 
G = tf(B,F,1); 
[Phiu,w] = pwelch(u);  % Estimate input power spectrum
[Gmag,~] = bode(G,w);
[Hmag,~] = bode(H,w);  

sigma_ehat = M.Report.Fit.MSE; % Estimate variance of e(t)
Phiv = squeeze(Hmag).^2*sigma_ehat; 
Phiy = squeeze(Gmag).^2.*Phiu+Phiv;

figure;
loglog(w,Phiv,'linewidth',2); 
hold on; 
loglog(w,Phiv./Phiu,'linewidth',2); 
loglog(w,Phiv./Phiy,'linewidth',2);

