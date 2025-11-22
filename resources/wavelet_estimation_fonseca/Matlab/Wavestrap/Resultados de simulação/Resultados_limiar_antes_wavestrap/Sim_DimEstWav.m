
n = 600; % tamanho amostral
d = 6; % parâmetro de dimensão
nt = 256; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
% N usado na função wavedec. Como temos 7 (2^7=256) níveis, ao selecionar
% N=3 teremos armazenados os níveis de detalhes 7, 6 e 5.
N = 3;
wname = 'db2';  % base de ondaletas utilizada
p = 5;  % lag máximo utilizado

J = 263;  % número de coefientes da base de ondaletas usada
NREP = 100;  % número de réplicas da simulação

Nboot = 100;  % número de réplicas bootstrap
alpha = .1;  % nível de significância do teste bootstrap
mPvalues = ones(NREP,8);  % valores-p dos testes dos 8 maiores autovalores
mEigval = zeros(NREP,8);  % 8 maiores autovalores de D
vDimSel = 9*ones(NREP,1);  % dimensões selecionadas em cada réplica

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4, quando d=4 são -.775,.65,-.525 e .4,
% e quado d=6 são -.8167,.7333,-.65,.5667,-.4833,.4
%model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',1.5);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);

tic;
for kk=1:NREP

% gerando os coeficientes do funcional não observado
xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
xi5 = simulate(model_xi5,n);xi6 = simulate(model_xi6,n);

% gerarndo o funcional não observado
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
       + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u)...
       + xi5(ii)*sqrt(2)*cos(5*pi*u) + xi6(ii)*sqrt(2)*cos(6*pi*u);
end

% gerando o ruído a ser adicionado ao funcional não observadao
mEps = zeros(nt,n);
for ii=1:n
    for jj=1:10
        mEps(:,ii) = mEps(:,ii) +...
            normrnd(0,1)*sqrt(2)*sin(pi*u*jj)/(2^(jj-1));
    end
end

% funcional observado
Y = X + mEps;

% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e
    % os coeficientes de detalhes
%    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,A(:,ii),~] = wden(Y(:,ii),'sqtwolog','h','mln',N,wname);
end
% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_A = mean(A,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = A - mu_A*ones(1,n);

% agora será obtida a matriz formada pelos coeficientes das
% decomposições das observações dos funcionais. 
C1 = C(:,[1:(n-p)]);
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,[(k+1):(n-p+k)])'*C(:,[(k+1):(n-p+k)]);
end
D = C1*D1*C1'/((n-p)^2);

% após obter D, calculamos os autovetores e autovalores da mesma.
% as colunas de B terão os autovetores de D, que são os coeficientes de
% ondaletas das autofunções associadas aos autovalores na diagonal de L
[B,L] = eig(D);
% 8 maiores autovalores de D
mEigval(kk,:) = diag(L([1:8],[1:8]));

% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_boot
    d_boot = DimEst_wavestrap( A, Nboot, B(:,[1:d0]), p);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(kk,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(kk,d0)>alpha)
        % se a hipótese for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel(kk) = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

end
toc;

vDimSel_n100 = vDimSel;mPvalues_n100 = mPvalues;mEigval_n100 = mEigval;
vDimSel_n300 = vDimSel;mPvalues_n300 = mPvalues;mEigval_n300 = mEigval;
vDimSel_n600 = vDimSel;mPvalues_n600 = mPvalues;mEigval_n600 = mEigval;
% salvando a matriz com os autovalores
save('file_Sim_comThresh_wavestrap_d6.mat','vDimSel_n100','mPvalues_n100','mEigval_n100',...
    'vDimSel_n300','mPvalues_n300','mEigval_n300',...
    'vDimSel_n600','mPvalues_n600','mEigval_n600');  
%file = matfile('file_Sim_comThresh_wavestrap_d6.mat');
%mDim_sim = example.mDim_sim;

tabulate(vDimSel)
mean(mPvalues)
mean(mEigval)
hist(d_boot)
L(d0+1,d0+1)



