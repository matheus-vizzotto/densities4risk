n = 600; % tamanho amostral
d = 6; % parâmetro de dimensão
nt = 500; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
du = u(2) - u(1);
p = 5;
vw = [.1,.1,.8]'; % vetor de pesos
vw = [1/3,1/3,1/3]'; % vetor de pesos
%vw = vw/norm(vw);
delta = length(vw); % número de pesos usados

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4, quando d=4 são -.775,.65,-.525 e .4,
% e quado d=6 são -.8167,.7333,-.65,.5667,-.4833,.4
%model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',1.5);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%var_ar = 2.5;
%model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',var_ar);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',var_ar);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',var_ar);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',var_ar);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',var_ar);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',var_ar);

%rng('default')
xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
xi5 = simulate(model_xi5,n);xi6 = simulate(model_xi6,n);

% gerarndo a função não observada
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
       + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u)...
       + xi5(ii)*sqrt(2)*cos(5*pi*u) + xi6(ii)*sqrt(2)*cos(6*pi*u);
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

mu_Y = mean(Y,2);

C = Y - mu_Y*ones(1,n);

% funcional observado centrado
mat_G = du*(C'*C);

% matriz com os produtos internos dos funcionais centrados multiplicados
% pelos pesos
Zk = zeros(n-p-delta+1,n-p-delta+1,p-delta+1);
for k=delta:(p-delta+1)
    mat_Z = FuncMatProdYcentAggregVet( C, vw, p, k);
    Zk(:,:,k) = du*mat_Z;
end

Kstar = sum(Zk,3)*mat_G(1:(n-p-delta+1),1:(n-p-delta+1))/((n-p-delta+1)^2);

[B,L] = eig(Kstar);
sprintf('%8.4f, ',diag(L(1:10,1:10)))


% calculando as autofunções como Bathia et al (2010) descreve na equação
% (2.13)
mEigFunc = zeros(nt,10);
for ii=1:10
    mEigFunc(:,ii) = du*C(:,1:(n-p-delta+1))*B(:,ii);
end
mEigFunc = MatGramSchm(mEigFunc);


NREP = 100; alpha = .05; vDimSel = 9;
d0 = 1;
mPvalues = ones(1,8);
tic;
while d0<=8
    d_boot = DimEst_Aggreg_du( Y, NREP, mEigFunc(:,1:d0), p, vw, du);
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

