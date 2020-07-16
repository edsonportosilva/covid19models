
close all 
clear all
clc

% Referências: 
% [1] R.M. Cotta, C.P. Naveira-Cotta, and P. Magal, Parametric identification 
%     and public health measures influence on the COVID-19 epidemic evolution 
%     in Brazil, 2020. 

% Sumário de variáveis utilizadas no modelo SIRU:
%
% S(t): número de indivíduos passíveis de seren infectados no dia t
% I(t): número de indivíduos infecciosos e assintomáticos no dia t
% R(t): número de indivíduos infecciosos e notificados no dia t
% U(t): número de indivíduos infecciosos e não-notificados no dia t
% CR(t): cumulativo de indivíduos infectados notificados até o dia t
% CU(t): cumulativo de indivíduos infectados não-notificados até o dia t
% f(t): fração de indivíduos assintomáticos se tornarão casos notificados
% 1-f(t): fração de indivíduos assintomáticos se tornarão casos não-notificados
% DR(t): número diário de indivíduos notificados

% Premissas:
%
% 1. Indíviduos infecciosos notificados R(t) são removidos (ou isolados) da
% população e não causam novas infecções.
% 2. Individuos assintomáticos (I) são infecciosos por um período de Ti dias.
% 3. Indivíduos sintomáticos e notificados, ou não,(R ou U) são infecciosos por um período
% de Tr dias;
% 4. Todas as infecções acontecem via indivíduos dos grupos I ou U.

% Parâmetros do modelo:
t    = 0:150;   % Intervalo de tempo em dias
Ti   = 7.024;      % Tempo médio em que o indivíduo em I permanece infeccioso
Tr   = 7.024;      % Tempo médio em que o indivíduo em R permanece infeccioso
f0   = 0.1257;    % Taxa inicial de notificação de casos
X1   = 2.9445;  % Parâmetro de ajuste da curva de casos acumulados
X2   = 0.1408;  % Parâmetro de ajuste da curva de casos acumulados
X3   = 3.5587;  % Parâmetro de ajuste da curva de casos acumulados
N    = 28.05;     % Dias sem medidas de intervenção na dinâmica da epidemia
Nf   = length(t)+1; % Dias sem medidas de intervenção na contagem de casos
mu   = 0.0119; % Parâmetro de ajuste da curva de taxa de transmissão
muf  = 0;     % Parâmetro de ajuste da curva de taxa de notificação
fmax = 0;     % Parâmetro de ajuste da curva de taxa de notificação

% Parâmetros derivados utilizados no modelo:
nu  = 1/Ti; 
eta = 1/Tr; 
f   = f0*ones(1, length(t)); 
nu1 = nu*f;
nu2 = nu*(1-f);

S = 3.944e6;                  % Número inicial de indivíduos passíveis de infecção
I = X2*X3/(nu*f(1));          % Número inicial de indivíduos infectados assintomáticos
U = (1-f(1))*(nu/(eta+X2))*I; % Número inicial de indivíduos infectados, sintomáticos e não notificados
R = 0;                        % Número inicial de indivíduos infectados, sintomáticos e notificados
tau0 = (X2+nu)/S*((eta+X2)/((1-f(1))*nu+eta+X2)); % Taxa de infecção inicial
tau  = tau0*ones(1, length(t));

CR = 0;    % Número cumulativo de casos notificados
CU = 0;    % Número cumulativo de casos não-notificados
DR = 0;    % Número de notificações diárias

R0 = (tau0*S(1)/nu)*(1+(1-f(1))*nu/eta); % Taxa básica de reprodução

% SIRU model
for indT = 1:length(t)-1 
    
    if indT > N
    % A partir do dia N, a taxa de transmissão tau é afetada por medidas de
    % contenção da epidemia:
        tau(indT) = tau0*exp(-mu*(indT-N));
    end
    if indT > Nf
    % A partir do dia Nf, a taxa de notificação f é afetada por mudanças na
    % metodologia ou quantidade de testes:
        f(indT)=(fmax-f0)*(1-exp(-muf*(indT-Nf)))+f(1);
    end
    
    % Integração numérica:
    dS  = -tau(indT)*S(indT)*(I(indT)+U(indT));
    dI  = -dS - nu*I(indT);
    dR  = nu1(indT)*I(indT)-eta*R(indT);
    dU  = nu2(indT)*I(indT)-eta*U(indT);
    dDR = nu1(indT)*f(indT)*I(indT);    
  
    S(indT+1)  = S(indT) + dS;
    I(indT+1)  = I(indT) + dI; 
    R(indT+1)  = R(indT) + dR; 
    U(indT+1)  = U(indT) + dU; 
    CR(indT+1) = CR(indT) + nu1(indT)*I(indT);
    CU(indT+1) = CU(indT) + nu2(indT)*I(indT);
    DR(indT+1) = CR(indT+1) - CR(indT);
end

figure, 
plotScale = 1000;
hold on, plot(t, I/plotScale,'x','markerSize',4,'linewidth',1)
         plot(t, R/plotScale,'-.','markerSize',4,'linewidth',1)
         plot(t, U/plotScale,'--','markerSize',4,'linewidth',1)
         plot(t, CR/plotScale,'-d','markerSize',4,'linewidth',1)
         plot(t, CU/plotScale,'-o','markerSize',4,'linewidth',1)
         plot(t, DR/plotScale,'-sq','markerSize',4,'linewidth',1)
         xlabel('Tempo (dias)')
         ylabel(['Casos \times' num2str(plotScale)])
         title(['Taxa básica de reprodução inicial: R_0 = ' num2str(R0) ' Notificação inicial de casos: f_0 = ' num2str(100*f0) '%'])

load('cumCasesCG.mat','cumCasesCG')
CRdata = cumCasesCG';
plot(CRdata/plotScale,'k-o');

legend('População de infectados pre-sint./assint. I(t)', 'Casos sintomáticos ativos notificados R(t)', 'Casos sintomáticos ativos não-notificados U(t)',...
       'Cumulativo de casos notificados CR(t)', 'Cumulativo de casos não-notificados CU(t)', 'Notificações de casos no dia DR(t)','Cumulativo de casos PB','Location','NorthWest')