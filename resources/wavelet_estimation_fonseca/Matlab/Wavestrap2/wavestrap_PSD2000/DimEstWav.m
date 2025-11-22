% Neste código simulamos uma amostra de uma série temporal de funções como
% descrito na estudo numérico de Bathia, para analisar as funções
% DimEst_boot e DimEst_wavestrap

n = 100; % tamanho amostral
d = 2; % parâmetro de dimensão
nt = 256; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
N = 6;  % nível usado na função wavedec
wname = 'db4';
p = 5;

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4, quando d=4 são -.775,.65,-.525 e .4,
% e quado d=6 são -.8167,.7333,-.65,.5667,-.4833,.4
model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',1.5);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%var_ar = 2.5;
%model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',var_ar);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',var_ar);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',var_ar);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',var_ar);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',var_ar);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',var_ar);

%rng('default')
xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
%xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
%xi5 = simulate(model_xi5,n);xi6 = simulate(model_xi6,n);

% gerarndo a função não observada
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u);
%       + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u);
%       + xi5(ii)*sqrt(2)*cos(5*pi*u) + xi6(ii)*sqrt(2)*cos(6*pi*u);
end

% gerando o ruído a ser adicionado na função não observada
mEps = zeros(nt,n);
for ii=1:n
    for jj=1:10
        mEps(:,ii) = mEps(:,ii) +...
            normrnd(0,1)*sqrt(2)*sin(pi*u*jj)/(2^(jj-1));
    end
end

% funcional observado
Y = X + mEps;

% número de coefientes da base de ondaletas usada
J = 294;
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);
%Aths = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e L
    % os coeficientes de detalhes
    %[D,L] = dwt(Y(:,ii),'db2');
    % os coeficientes são empilhados na i-ésima coluna de A
    %A(:,ii) = [D',L']';
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    %[thr,sorh,keepapp] = ddencmp('den','wv',Y(:,ii));
    %[~,Aths(:,ii),Lwths,~,~] = wdencmp('gbl',A(:,ii),Lw,'db1',1,thr,sorh,keepapp);
%    [~,A(:,ii),Lw] = wden(Y(:,ii),'sqtwolog','h','mln',N,'db2');
end
%A = Aths;
% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_A = mean(A,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = A - mu_A*ones(1,n);

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

NREP = 100; alpha = .05; vDimSel = 9;
d0 = 1;
mPvalues = ones(1,8);
tic;
while d0<=8
    d_boot = DimEst_wavestrap( A, NREP, B(:,[1:d0]), p);
%    d_boot = DimEst_wavestrap2( Y, NREP, B(:,[1:d0]), p, N, 'db2');
%    d_boot = DimEst_wavestrap2Residuals( Y, NREP, B(:,[1:d0]), p, N, 'db2');
%    d_boot = DimEst_wavestrap2SomeResiduals( Y, NREP, B(:,[1:d0]), p, N, 'db2');
%    d_boot = DimEst_wavestrap2Residuals_teste( Y, NREP, B(:,[1:d0]), p, N, 'db2');
%    d_boot = DimEst_wavestrap_PSD2000( A, Lw, NREP, B(:,[1:d0]), p);
%    d_boot = DimEst_wavestrap_TrueResiduals( A, Lw, NREP, B(:,[1:d0]), p, nt, N, wname);
%    d_boot = DimEst_wavestrap_0MeanResiduals( A, NREP, B(:,[1:d0]), p);
%    d_boot = DimEst_wavestrap_ResidualsDimPlus1( A, NREP, B(:,[1:(d0+1)]), p);
    mPvalues(d0) = sum(d_boot>L(d0+1,d0+1))/NREP;
    if (mPvalues(d0)>alpha)
        vDimSel = d0;
        d0 = 8;
    end
    d0 = d0 + 1;
end
toc;
vDimSel
mPvalues



ind_boot = zeros(263,1);
ind_boot(1:34) = 1:34;
for k=2:4
    cs_i = cumsum(Lw(1:(k-1))); cs_i = cs_i(end);
    cs_f = cumsum(Lw(1:k)); cs_f = cs_f(end);
%    ind_boot((cont_aux + 1):(cont_aux + Lw(k)))
%    ind_boot((cumsum(Lw(1:(k-1)))+1):cumsum(Lw(1:k))) = unidrnd(Lw(k),Lw(k),1) + cumsum(Lw(1:(k-1)));
    ind_boot((cs_i+1):cs_f) = unidrnd(Lw(k),Lw(k),1) + cs_i;
end

A_boot = A(ind_boot,:);
mu_A_boot = mean(A_boot,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = A_boot - mu_A_boot*ones(1,n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d0 = 2;
tic;
d_boot = DimEst_wavestrap2( Y, NREP, B(:,[1:d0]), p, N, 'db2');
toc;
sum(d_boot>L(d0+1,d0+1))/NREP
sum(d_boot>L(d0+1,d0+1)) <= (alpha*NREP)

figure
hist(d_boot)
sprintf('%8.2f',L(d0+1,d0+1))


NREP = 50; alpha = .1; vDimSel = 1;
d0 = 1;
mPvalues = ones(1,8);
tic;
while d0<=8
    %d_boot = DimEst_boot( Y, NREP, B(:,[1:d0]), p, N, 'db2');
    %d_boot = DimEst_wavestrap( A, NREP, B(:,[1:d0]), p);
    %d_boot = DimEst_wavestrap2( Y, NREP, B(:,[1:d0]), p, N, 'db2');
    d_boot = DimEst_wavestrap3( A, NREP, B(:,[1:d0]), p, N, 'db2');
    %d_boot = DimEst_wavestrap4( A, Lw, NREP, B(:,[1:d0]), p, N, 'db2');
    mPvalues(d0) = sum(d_boot>L(d0+1,d0+1))/NREP;
    if (mPvalues(d0)>alpha)
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end
toc;
vDimSel
mPvalues


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Realizando a análise de dimensão fazendo thresholding

J = 263;
Aths = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    %[A(:,ii),Lw] = wavedec(Y(:,ii),N,'db2');
    [~,Aths(:,ii),~] = wden(Y(:,ii),'sqtwolog','h','mln',N,'db2');
end
mu_Aths = mean(Aths,2);

C = Aths - mu_Aths*ones(1,n);
C1 = C(:,[1:(n-p)]);
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,[(k+1):(n-p+k)])'*C(:,[(k+1):(n-p+k)]);
end
D = C1*D1*C1'/((n-p)^2);

[B,L] = eig(D);

NREP = 50; alpha = .1; vDimSel = 1;
d0 = 1;
mPvalues = ones(1,8);
tic;
while d0<=8
    %d_boot = DimEst_wavestrap( A, NREP, B(:,[1:d0]), p);
    d_boot = DimEst_wavestrap3( A, NREP, B(:,[1:d0]), p, N, 'db2');
    mPvalues(d0) = sum(d_boot>L(d0+1,d0+1))/NREP;
    if (mPvalues(d0)>alpha)
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end
toc;
vDimSel
mPvalues


