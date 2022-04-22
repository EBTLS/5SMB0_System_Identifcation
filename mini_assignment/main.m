%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TUE 5SMB0 System Identification
%%% Assignment 2022
%%% Authors: Jiaxuan Zhang, Yiting Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data Generation
% [ u , y ] = assignment_sys_36(r)
clear all
close all 
clc


%% Part 1: Understanding Saturation and Butterworth Filter

%% 1.1. Specification of the Butterworth Filter
F.num = [ 0.505 , 1.01 , 0.505 ];
F.denom = [ 1 , 0.7478 , 0.2722 ];
F.sys = tf(F.num, F.denom, -1, 'Variable', 'z^-1')

plotopts1 = bodeoptions;
plotopts1.FreqScale = 'linear';
plotopts1.XLim = {[0,4]};
plotopts1. Ylim = {[-20,20]};
plotopts1.Grid = 'on';
bode(F.sys, plotopts1)
% w = 2.1954 rad/s is the -3 db point, f = w/2/pi = 0.3494


%% 1.2. M Specification
% nu = 100;
% previous_u = -100; 
% for mag=0:1:100
%     r = mag * ones(1,nu);
%     [u,y] = assignment_sys_36(r);
%     u(1)
%     previous_u
%     if (u(1) == previous_u)
%         fprintf("M is equal to %d", u(1));
%         break;
%     end
%     previous_u = u(1);
% end
% 
% % M is equal to 3

N = 2000;
f = 10;
offset = 0;
A = 1;

t = (0 : 1/N : 1)';
r = A * sin(2*pi*f*t) + offset; 

% r = ones([1, N]);
[u,y] = assignment_sys_36(r);


figure;
subplot(1, 2, 1)
plot(u)
grid on;
xlabel("t");
ylabel("u");
title("fig");

subplot(1, 2, 2)
plot(y)
grid on;
xlabel("t");
ylabel("y");
title("fig");
% we now can know M is 3

%% Part 2: Nonparametric identification

%% 2.1. Nonparametric Identification
% system setup
N = 1024;
Ts = 1;
n_freq = 128;

% set the base frequency
w0 = 2/128;

% generate reference signal
t = [0: 1: N-1]/Ts;
r=0;
for i = 1: 1: n_freq
    r = r + 3 * sin(w0*i*t) + offset;
end

% analysis reference spectrum
w = fft(r);
w = circshift(w, N/2);
wf = [-N/2+1: N/2]*w0/2;
plot(wf, abs(w))
xlabel("w"); ylabel("amplitude"); title("FFT of reference signal");
grid on;

% plot reference signal
figure;
subplot(1, 2, 1)
plot(r)
grid on
xlabel("t"); ylabel("r"); title("reference signal");

% generate input-output
[u,y] = assignment_sys_36(r);
subplot(1, 2, 2)
plot(u)
grid on
xlabel("t"); ylabel("u"); title("input signal");

% ETFE identification
figure;
G_ETFE = etfe(iddata(y,u))
bode(G_ETFE)
grid on
title("bode diagram of the estimated model from ETFE");

% Check input signal spectrum
figure
w=fft(u);
plot(abs(w))
grid on;
title("DFT of input signal")

%% 2.3. Noise Spectrum Analysis
% Design a zero input, that is the output y will totally be noise signal
period = 128;
N = 1024;
r = 0*sin( [0: 1: N-1] * pi / period );

figure;
plot(r)
grid on
xlabel("t"); ylabel("r"); title("reference signal");

figure;
subplot(1, 2, 1)
[u,y] = assignment_sys_36(r);
plot(u)
grid on
xlabel("t"); ylabel("u"); title("input signal");

subplot(1, 2, 2)
plot(y)
grid on
xlabel("t"); ylabel("y"); title("output signal");

% analysis power specturm of noise signal
v = y;
figure;
[Pv, W] = cpsd(v, v, 1024);
cpsd(v, v, 1024)

%% Part 3: Expreiment Design

%% 3.1. Input Choice
% prefer PRBS

%% 3.2. Design the Input
[r,P,w] = Generate_PRBS(0.83, 3000, 3, true);
% I think 0.83 is a good threshold

%% Part 4: Parametric identification and validation

%% 4.1 OE Model

% % First Try OE Model
% % divide the dataset into training set and validation set
% [r,P,w] = Generate_PRBS(0.83, 3000, 3, false);
% [u,y] = assignment_sys_36(r);
% train_ratio = 0.7;
% u_train = u([1: floor(train_ratio * length(u))]);
% y_train = y([1: floor(train_ratio * length(y))]);
% u_test = u([floor(train_ratio * length(u))+1: end]);
% y_test = y([floor(train_ratio * length(y))+1: end]);
% 
% % prepare loop for order selection
% order_upper = 20;
% nb = 0; nf = 0; nk = 0;
% 
% OE.order.order_points = [];
% OE.order.order_config = [];
% OE.order.error = [];
% 
% for order_sum = 1: 1: order_upper
%     
%     for nb = 0: 1: order_sum
%         
%         for nf = 0: 1: order_sum - nb
%             
%             nk = order_sum - nb - nf;
%             % for OE model, nk and nb cannot be zero at the same time
%             if nk == 0 && nb == 0
%                 continue;
%             end
%             
%             % identify a model
%             temp_sys = oe(iddata(y_train, u_train), [nb, nf, nk]);
%             % cost validate error
%             temp_predict_error = pe(temp_sys, iddata(y_test, u_test), 1);
%             temp_error = temp_predict_error.y' * temp_predict_error.y / length(y_test);
%             
%             % abort too abnormal result
%             if temp_error > 50
%                 continue;
%             end
%             
%             OE.order.order_points = [OE.order.order_points, order_sum];
%             OE.order.order_config = [OE.order.order_config, [nb; nf; nk]];
%             OE.order.error = [OE.order.error, temp_error];
%             
%         end
%                   
%     end
% end
% scatter(OE.order.order_points, OE.order.error);


% % choose the optimal model
% [M, I] = min(OE.order.error);
% OE.model_configuration = OE.order.order_config(:, I);
% OE.model_order = OE.order.order_points(I);
% OE.model_configuration

% generate input-output signal
[r,P,w] = Generate_PRBS(0.83, 3000, 3, false);
[u,y] = assignment_sys_36(r);

% identify an OE model
OE.sys = oe(iddata(y, u), [3,13,1]);

% residual test
OE.resid = resid(iddata(y,u), OE.sys);
figure
resid(iddata(y,u), OE.sys)
grid on;

%% 4.1 ARX Model

% identify an ARX model
G_ARX = arx([y u], [13 3 1])

% do the residual test
resid([y u], G_ARX)


% Using the systemIdentification Toolbox 
% Using the r to generate u and y, then import these two time domian data
% into the Toolbox
% Using the polynomial model estimation with ARX model
% In the Polynomial Model screen

%% Part 5: Experimental verification of variance estimates

%% 5.1. Variance Esitmate
% using the same r to generate the u and y, set the loop number as 100
Loop = 100;

[r,P,w] = Generate_PRBS(0.83, 1000, 3, false);

% using the parameter obtained from previous section
nb = 4; nf = 14; nk = 1;
result = zeros([Loop, nb+nk+nf-1], 'double');

for i = 1:Loop
    [u,y] = assignment_sys_36(r);
    G_OE = oe([y u], [nb nf nk]);
    F = G_OE.Structure.F.Value(2:1:end);
    B = G_OE.Structure.B.Value(2:1:end);
    result(i,:) = [B F];
end
clear B F i
% the storing fromat is B(1x4) and F(1x13) in one raw
% result

%% 5.2. Inspect the Parameters obtained from the Monte Carlo Simulations
% due to the existance of the noise, each realization is different from
% each other

%% 5.3. Variance Comparison
cov = getcov(G_OE);

cov_est = diag(cov)';

cov_sim = zeros([2, (nb+nk+nf-1)], 'double');
for k = 1:(nb+nk+nf-1)
     
    cov_sim(1, k) = mean(result(:, k));
    cov_sim(2, k) = var(result(:, k));
end

figure;
plot(cov_sim(2,:))
hold on
plot(cov_est)
legend('simulation','estimation')
xlabel("t"); ylabel("y"); title("Variance comparison");

clear k

%% 5.4. Monte Carlo with Median Starting Point
result_54 = zeros([Loop, nb+nk+nf-1], 'double');

% for i = 1:Loop
%     [u,y] = assignment_sys_36(r);
%     G_OE = oe([y u], [nb nf nk]);
%     F_init = G_OE.Structure.F.Value;
%     B_init = G_OE.Structure.B.Value;
%     M_init = idpoly([], B_init, [], [], F_init);
%     M_oe = oe([y u], M_init);
%     F = G_OE.Structure.F.Value(2:1:end);
%     B = G_OE.Structure.B.Value(2:1:end);
%     result_54(i,:) = [B F];
% end

OEopt = oeOptions;
OEopt.SearchOptions.Tolerance = 0.0001;

for i = 1:Loop
    
    % Generate new data
    [u,y] = assignment_sys_36(r);
    
    if i ==1
        
        % generate an initial OE model
        M_oe = oe([y,u], [nb, nf, nk], OEopt);
        
    else
        % update starting point
        B_init = median(result_54([1: 1: i], [1: 1: nb]),1);
        F_init = median(result_54([1: 1: i], [nb+1: 1: end]),1);
        M_init = idpoly([], [1, B_init], [], [], [1, F_init]);
        
        % start OE model identification from the starting point
        M_oe = oe([y u], M_init, OEopt);
    end
    
    % record the parameter
    F = G_OE.Structure.F.Value(2:1:end);
    B = G_OE.Structure.B.Value(2:1:end);
    result_54(i,:) = [B F];
           
end
    

clear B F i
% the storing fromat is B(1x4) and F(1x13) in one raw
% result_54

cov_sim_54 = zeros([2, (nb+nk+nf-1)], 'double');
for k = 1:(nb+nk+nf-1)
     
    cov_sim_54(1, k) = mean(result_54(:, k));
    cov_sim_54(2, k) = var(result_54(:, k));
end

figure;
plot(cov_sim_54(2,:))
hold on
plot(cov_est)
legend('simulation','estimation')
xlabel("t"); ylabel("y"); title("Variance Comparison (median mode)");

clear k Loop

%% Part 6: Estimation of a Box Jenkins model for minimum variance

%% 6.1. Identify an BJ model
nb = 4;
nf = 14;
nc = 3;
nd = 5;
nk = 1;

% identify an BJ model
BJ.sys = bj(iddata(y,u), [nb, nc, nd, nf, nk]);

% resid test
figure;
resid(iddata(y,u), BJ.sys)
grid on

% bode graph of ETFE and identifed BJ model
figure;
bode(BJ.sys, G_ETFE)
legend('BJ','ETFE')
grid on

% bode graph of the noise signal power spectrum and identified H in the

% bode plot options
plotopts2 = bodeoptions;
plotopts2.FreqScale = 'linear';
plotopts2.XLim = {[0,4]};
% plotopts2. Ylim = {[-50,50]};
plotopts2.Grid = 'on';

% bode graph of identified H
figure
H = tf(BJ.sys.C, BJ.sys.D, -1, 'Variable', 'z^-1');
bode(H,plotopts2)
grid on;
title('bode graph of H1');

% noise spectrum 
figure
cpsd(v, v, 1024);
title('noise sepctrum')

%% 6.2. compute the variance
cov = getcov(BJ.sys);
figure 
plot(cov(2,:))

