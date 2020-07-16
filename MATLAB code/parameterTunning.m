
clear
close all

load('cumCasesCG.mat','cumCasesCG')

CRdata = cumCasesCG';
t = [1:30];

% Set Options
options = optimset('GradObj', 'off', 'MaxIter', 800);

CRdata_ = CRdata(1:30);
initial_theta = [0.01 0.01 1];
% Optimize
[theta, J, exit_flag] = fminunc(@(X)(calcMSE(CRdata_, t, X)), initial_theta, options);

X = theta;
CRtest = X(1)*exp(X(2)*t)-X(3);

figure, plot(CRdata,'-o');
figure(1),hold on, plot(CRtest,'-x')


%% Test SIRU model

clear 
close all
clc

load('cumCasesCG.mat','cumCasesCG')

CRdata = cumCasesCG';
% Set Options
options = optimset('GradObj', 'off','MaxFunEvals',800, 'MaxIter', 800,'Display','on');

t = [1:117];
CRdata_ = CRdata(1:max(t));
initial_theta = [7 0.15 28 0.0001];
% Optimize
[theta, J, exit_flag] = fminunc(@(param)(funSIRUmodel(CRdata_, t, param)), initial_theta, options);

[~, CRtest] = funSIRUmodel(CRdata_, t, theta);

figure, plot(CRdata,'-o');
figure(1),hold on, plot(CRtest,'-x')

figure, plot(diff(CRdata),'-o')