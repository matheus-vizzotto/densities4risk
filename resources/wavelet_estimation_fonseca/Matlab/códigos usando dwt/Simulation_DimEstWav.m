
n = 100; % tamanho amostral
d = 4; % parâmetro de dimensão
nt = 256; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados

% processos para gerar a dependência temporal da função não observada
model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',.01);
model_xi2 = arima('Constant',0,'AR',{.65},'Variance',.01);
model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',.01);
model_xi4 = arima('Constant',0,'AR',{.4},'Variance',.01);

%rng('default')
xi1 = simulate(model_xi1,n);
xi2 = simulate(model_xi2,n);
xi3 = simulate(model_xi3,n);
xi4 = simulate(model_xi4,n);

% gerarndo a função não observada
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
       + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u);
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
J = 261;
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e L
    % os coeficientes de detalhes
    %[D,L] = dwt(Y(:,ii),'db2');
    % os coeficientes são empilhados na i-ésima coluna de A
    %A(:,ii) = [D',L']';
    A(:,ii) = wavedec(Y(:,ii),2,'db2');
end

% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_A = mean(A,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = A - mu_A*ones(1,n);

% agora será obtida a matriz formada pelos coeficientes das
% decomposições das observações dos funcionais. Aqui, tomamos
% p=5.
p = 5;

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
plot(diag(L(1:10,1:10)))
sprintf('%8.2f',diag(L(1:10,1:10)))

norm(diag(L(1:4,1:4)),1)/norm(diag(L),1)

% tomando a transformada inversa da função média
mu_hat =  idwt(mu_A(1:(J/2)),mu_A((J/2 + 1):end),'db2');

subplot(2,1,1); plot(Y(:,1));
subplot(2,1,2); plot(X(:,1));
subplot(2,1,2); plot(mu_hat);


NREP = 100;
alpha = .05;
d0 = 5;

tic;
d_boot = DimEst_boot( Y, NREP, B(:,1:d0), p, 'db2');
toc;
sum(d_boot>L(d0+1,d0+1)) <= (alpha*NREP)

hist(d_boot)
L(d0+1,d0+1)

