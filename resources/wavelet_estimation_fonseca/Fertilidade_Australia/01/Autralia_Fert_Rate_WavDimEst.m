
% Taxas anuais de fertilidade na Australia de 1921 a 2000
FertRate = csvread('Australia_Fert_Rate.csv');
% logaritmo das taxas
logFerRate = [FertRate(1,:);FertRate(2:10,1) log(FertRate(2:10,2:81))];

% log das taxas de fertilidade vistos como séries temporais para cada faixa
% de idade
plot(logFerRate(1,2:81),logFerRate(3:9,2:81))
%plot(1921:2000,logFerRate(3,2:81),'-',1921:2000,logFerRate(4,2:81),'--',...
%1921:2000,logFerRate(5,2:81),':',1921:2000,logFerRate(6,2:81),'-.',...
%1921:2000,logFerRate(7,2:81),'-',1921:2000,logFerRate(8,2:81),'--',...
%1921:2000,logFerRate(9,2:81),':','LineWidth',2)

set(gca,'LineStyleOrder',{'-','--',':','-.','-','--',':'})
legend('<20','20-24','25-29','30-34','35-39','40-44','>44')
legend('Location','southwest')
legend('boxoff')

n = 80;  % número de funcionais
nt = 64;  % número de pontos usados para cada ano
u = linspace(13,52,nt);  % nt valores igualmente espaçados de 13 a 52
% matriz para armazenar as funções suavizadas do log das taxas anuais
% avaliadas nos pontos de u referente às idades das mães
SmoLogFertRate = zeros(nt,n);

%fspl = fit(logFerRate(2:10,1),logFerRate(2:10,2),'smoothingspline');
%plot(fspl,logFerRate(2:10,1),logFerRate(2:10,2))
%plot(u,fspl(u))

% laço para obter a função suavizada da log-taxa em cada ano, sendo
% armazenada a função avaliada nos pontos do vetor u
for k=2:81
    fspl = fit(logFerRate(2:10,1),logFerRate(2:10,k),'smoothingspline');
    SmoLogFertRate(:,(k-1)) = fspl(u);
end

% log-taxas em função das idades para alguns anos
%plot(u,SmoLogFertRate(:,[1,20,40,60,80]),'LineWidth',2)
plot(u,SmoLogFertRate(:,1),'-',u,SmoLogFertRate(:,20),'--',...
    u,SmoLogFertRate(:,40),':',u,SmoLogFertRate(:,60),'-.',...
    u,SmoLogFertRate(:,80),'--','LineWidth',2)
xlabel('Age')
ylabel('Log fertility rate')
legend('1921','1940','1960','1980','2000','Position',[120 80 1 1])
%legend('Location','southwest')
legend('boxoff')


%*******************************************************************%
%*******************************************************************%

% Nessa etapa será feita uma análise de componentes principais funcional
% para as log-taxas em função das idades, onde serão usadas os funcionais
% discretizados (avaliados em pontos específicos), tomando como base o PCA
% descrito em Ramsay e Silverman (2005, p. 162 e 163). A base usada na
% expansão dos funcionais será a base de ondaletas.

% média cruzada (ao longo dos anos) das log-taxas
mu = mean(SmoLogFertRate,2);
% matriz com funcionais discretizados e centrados (com média zero)
X = (SmoLogFertRate - mu*ones(1,n))';

N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
% número de coefientes da base de ondaletas usada para o nt utilizado
J = 84;
% matriz de cujas linhas terão a decomposição em ondaletas
% de cada funcional observado (linha da matriz X)
Xdec = zeros(n,J);

% decomposição de cada funcional dos n anos
for ii = 1:n
    % As linhas de Xdec vão armazenar as decomposições das linhas de X,
    % i.e., dos funcionais. Lw vai armazenar os númeroes de coeficientes em
    % cada nível das decomposições, sendo seu último elemente o número de
    % pontos das colunas
    [Xdec(ii,:),Lw] = wavedec(X(ii,:),N,wname);
end

% versão discretizada do operador de variância-covariância
mV = (1/n)*(Xdec')*Xdec;

% Análise de autovalores e autovetores. As colunas da matriz B serão os
% coeficientes da decomposição das autofunções discretizadas, enquanto que
% os valores na diagonal de L serão seus autovalores correspondentes
[B,L] = eig(mV);

plot(diag(L(1:5,1:5)))  % gráfico dos maiores autovalores
sprintf('%8.2f, ',diag(L(1:5,1:4)))  % maiores autovalores
% porcentagem da variabilidade explicada por cada autofunção
100*diag(L(1:4,1:4))/norm(diag(L),1)

% recuperando autofunções avaliadas nos pontos do vetor u
h1_hat =  waverec(B(:,1),Lw,wname);
h2_hat =  waverec(B(:,2),Lw,wname);
h3_hat =  waverec(B(:,3),Lw,wname);

subplot(1,4,1); plot(u,mu);xlim([12 53])
subplot(1,4,2); plot(u,h1_hat);xlim([12 53])
subplot(1,4,3); plot(u,h2_hat);xlim([12 53])
subplot(1,4,4); plot(u,h3_hat);xlim([12 53])

% coeficientes da decomposição dos funcionais nas três primeiras
% autofunções
vbeta1 = X*h1_hat;
vbeta2 = X*h2_hat;
vbeta3 = X*h3_hat;
% série temporal formada pelos coeficientes da decomposição
subplot(1,3,1); plot(1921:2000,vbeta1)
subplot(1,3,2); plot(1921:2000,vbeta2)
subplot(1,3,3); plot(1921:2000,vbeta3)

%*******************************************************************%
%*******************************************************************%


p = 5;  % lag máximo utilizado
N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
% número de coefientes da base de ondaletas usada para o nt utilizado
J = 84;
% matriz de cujas linhas terão a decomposição em ondaletas
% de cada funcional observado (linha da matriz X)
Xdec = zeros(J,n);

% decomposição de cada funcional dos n anos
for ii = 1:n
    % As linhas de Xdec vão armazenar as decomposições das linhas de X,
    % i.e., dos funcionais. Lw vai armazenar os númeroes de coeficientes em
    % cada nível das decomposições, sendo seu último elemente o número de
    % pontos das colunas
    [Xdec(:,ii),Lw] = wavedec(SmoLogFertRate(:,ii),N,wname);
end

% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_dec = mean(Xdec,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = Xdec - mu_dec*ones(1,n);

% agora será obtida a matriz formada pelos coeficientes das
% decomposições das observações dos funcionais. 
C1 = C(:,1:(n-p));
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);

% após obter D, calculamos os autovetores e autovalores da mesma.
% as colunas de B terão os autovetores de D, que são os coeficientes de
% ondaletas das autofunções associadas aos autovalores na diagonal de L
[B,L] = eig(D);

plot(diag(L(1:5,1:5)))
sprintf('%8.4f, ',diag(L(1:10,1:10)))
diag(L(1:5,1:5))'
round(100*real(diag(L(1:5,1:5)))/norm(diag(L)),4)'

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

% Usando uma semente para o teste bootstrap
rng(2018);
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap( Xdec, Nboot, B(:,[1:d0]), p);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel

mu =  waverec(mu_dec,Lw,wname);
% recuperando autofunções avaliadas nos pontos do vetor u
h1_hat =  waverec(B(:,1),Lw,wname);
h2_hat =  waverec(B(:,2),Lw,wname);
h3_hat =  waverec(B(:,3),Lw,wname);
h4_hat =  waverec(B(:,4),Lw,wname);
%h5_hat =  waverec(B(:,5),Lw,wname);

subplot(2,3,1); plot(u,mu);xlim([12 53]); 
ylabel('$\mu(x)$','Interpreter','latex');xlabel('Age');
subplot(2,3,2); plot(u,h1_hat);xlim([12 53])
ylabel('$h_1(x)$','Interpreter','latex');xlabel('Age');
subplot(2,3,3); plot(u,h2_hat);xlim([12 53])
ylabel('$h_2(x)$','Interpreter','latex');xlabel('Age');
subplot(2,3,4); plot(u,h3_hat);xlim([12 53])
ylabel('$h_3(x)$','Interpreter','latex');xlabel('Age');
subplot(2,3,5); plot(u,h4_hat);xlim([12 53])
ylabel('$h_4(x)$','Interpreter','latex');xlabel('Age');


%*******************************************************************%
%*******************************************************************%

% Nessa simulação usamos o terceiro tipo de bootstrap, em que uma
% limiarização é aplicada na função média e nas autofunções para cirar um
% resíduo que será usado para formar as amostras bootstrap

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

% Usando uma semente para o teste bootstrap
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap_thresh(Xdec, Nboot, B(:,[1:d0]), p, N, wname);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel


%*******************************************************************%
%*******************************************************************%

% Nessa simulação usamos o segundo tipo de bootstrap, em que uma
% limiarização é aplicada nos funcionais observados antes de realizar um
% procedimento de bootstrap de Bathia et al (2010)

p = 5;  % lag máximo utilizado
N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
% número de coefientes da base de ondaletas usada para o nt utilizado
J = 84;
% matriz de cujas linhas terão a decomposição em ondaletas
% de cada funcional observado (linha da matriz X)
Xdec = zeros(J,n);

% decomposição de cada funcional dos n anos
for ii = 1:n
    % As linhas de Xdec vão armazenar as decomposições das linhas de X,
    % i.e., dos funcionais. Lw vai armazenar os númeroes de coeficientes em
    % cada nível das decomposições, sendo seu último elemente o número de
    % pontos das colunas
    [~,Xdec(:,ii),Lw] = wden(SmoLogFertRate(:,ii),'sqtwolog','h','mln',N,wname);
end

% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_dec = mean(Xdec,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = Xdec - mu_dec*ones(1,n);

% agora será obtida a matriz formada pelos coeficientes das
% decomposições das observações dos funcionais. 
C1 = C(:,1:(n-p));
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);

% após obter D, calculamos os autovetores e autovalores da mesma.
% as colunas de B terão os autovetores de D, que são os coeficientes de
% ondaletas das autofunções associadas aos autovalores na diagonal de L
[B,L] = eig(D);

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap( Xdec, Nboot, B(:,[1:d0]), p);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel

