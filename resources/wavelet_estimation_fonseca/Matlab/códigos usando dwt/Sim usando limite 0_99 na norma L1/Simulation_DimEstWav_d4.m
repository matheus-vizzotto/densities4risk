
n = 600; % tamanho amostral
d = 4; % parâmetro de dimensão
nt = 256; % número de pontos para cada dia
J = 258;  % número de coefientes da base de ondaletas usada
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados

NREP = 200;
% essa matriz armazenará nas linhas os 10 maiores autovalores da matriz
% B em cada iteração
mDim_sim = zeros(NREP,J);
% esse vetor armazenará qual foi a dimensão selecionada em cada iteração
vDim_esc = zeros(NREP,1);

for kk=1:NREP
    
% processos para gerar a dependência temporal da função não observada
model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',.01);
model_xi2 = arima('Constant',0,'AR',{.65},'Variance',.01);
model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',.01);
model_xi4 = arima('Constant',0,'AR',{.4},'Variance',.01);

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

% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e L
    % os coeficientes de detalhes
    [D,L] = dwt(Y(:,ii),'db2');
    % os coeficientes são empilhados na i-ésima coluna de A
    A(:,ii) = [D',L']';
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

% essa matriz armazenará nas linhas os 10 maiores autovalores da matriz
% B em cada iteração
mDim_sim(kk,:) = diag(L)';
end

% salvando a matriz com os autovalores
save('mDim_sim_d4_n600.mat','mDim_sim');  
%example = matfile('mDim_sim_d4.mat');
%mDim_sim = example.mDim_sim;

for kk = 1:NREP
ii = 1;
while ii<11
    if norm(mDim_sim(kk,1:ii))/norm(mDim_sim(kk,:)) > .99;
        vDim_esc(kk) = ii;
        ii = 11;
    else
        ii = ii + 1;
    end
end
end

sprintf('%8.2f',mean(mDim_sim(:,[1:10])))
plot(mean(mDim_sim(:,[1:10])))
tabulate(vDim_esc)
