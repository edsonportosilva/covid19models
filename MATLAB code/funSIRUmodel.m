function [mse, CR] = funSIRUmodel(CRdata, t, param)

% Parâmetros do modelo:
%t    = 0:120;   % Intervalo de tempo em dias
Ti   = param(1); % Tempo médio em que o indivíduo em I permanece infeccioso
Tr   = param(1); % Tempo médio em que o indivíduo em R permanece infeccioso
f0   = param(2); % Taxa inicial de notificação de casos
X1   = 2.9445;   % Parâmetro de ajuste da curva de casos acumulados
X2   = 0.1408;   % Parâmetro de ajuste da curva de casos acumulados
X3   = 3.5587;   % Parâmetro de ajuste da curva de casos acumulados
N    = param(3);     % Dias sem medidas de intervenção na dinâmica da epidemia
Nf   = length(t)+1;  % Dias sem medidas de intervenção na contagem de casos
mu   = param(4);     % Parâmetro de ajuste da curva de taxa de transmissão
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

% Calcula MSE entre o curva do modelo e os dados de referência:
mse  = mean(abs(CRdata-CR).^2);

end