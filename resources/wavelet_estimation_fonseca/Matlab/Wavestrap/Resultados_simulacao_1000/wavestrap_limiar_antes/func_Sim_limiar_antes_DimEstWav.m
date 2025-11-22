function [ vDimSel, mEQMI_n, mEQMI_nt, mEigval, mPvalues ] = func_Sim_limiar_antes_DimEstWav(n,d,NREP,Nboot)
% Essa função é utilizada para simular a estimação da dimensão de uma série
% temporal de curvas seguindo como fizeram Bathia et al (2010). As entradas
% são o tamanho amostral (número de curvas observadas) n, a dimensão do
% processo d, o número de réplicas da simulação NREP, o número de réplicas
% bootstrap Nboot no teste da dimensão.

% Parâmetros já definidos da simulação
nt = 256; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
% N usado na função wavedec. Como temos 7 (2^7=256) níveis, ao selecionar
% N=3 teremos armazenados os níveis de detalhes 7, 6 e 5.
N = 3;
wname = 'db4';  % base de ondaletas utilizada
p = 5;  % lag máximo utilizado
J = 276;  % número de coefientes da base de ondaletas usada
alpha = .05;  % nível de significância do teste bootstrap

% Outputs da função
mPvalues = ones(NREP,8);  % valores-p dos testes dos 8 maiores autovalores
mEigval = zeros(NREP,8);  % 8 maiores autovalores de D
vDimSel = 9*ones(NREP,1);  % dimensões selecionadas em cada réplica
mEQMI_n = zeros(NREP,n);  % EMQI das curvas estimadas ao longo dos tempos
% EQMI das curvas estimadas ao longo do pontos avaliados
mEQMI_nt = zeros(NREP,nt);


% processos para gerar a dependência temporal da função não observada
X = zeros(nt,n);
if d==2
    model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);
else
    if d==4
        model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
    else
        % deve ser o caso d==6
        model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',1.5);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);
    end
end

for kk=1:NREP

if d==2
    % gerando os coeficientes do funcional não observado
    xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
    % gerarndo o funcional não observado
    for ii=1:n
        X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u);
    end
else
    if d==4
        xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
        for ii=1:n
            X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
               + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u);
        end
    else
        % deve ser o caso d==6
        xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);xi5 = simulate(model_xi5,n);xi6 = simulate(model_xi6,n);
        for ii=1:n
            X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
               + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u)...
               + xi5(ii)*sqrt(2)*cos(5*pi*u) + xi6(ii)*sqrt(2)*cos(6*pi*u);
        end
    end
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
    %[A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,A(:,ii),Lw] = wden(Y(:,ii),'sqtwolog','h','mln',N,wname);
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

H = zeros(nt,vDimSel(kk));  % matriz para armazenar as autofunções
Yhat = zeros(nt,n);  % matriz para armazenar os funcionais estimados
% tomando a transformada inversa da função média
mu_hat =  waverec(mu_A,Lw,wname);
% recuperando as autofunções a partir dos seus coeficientes de ondaletas
for ii=1:(vDimSel(kk))
    H(:,ii) = waverec(B(:,ii),Lw,wname);
end
% recuperando os funcionais estimados
for ii=1:n
    Yhat(:,ii) = mu_hat;
    for k=1:(vDimSel(kk))
        Yhat(:,ii) = Yhat(:,ii) + (C(:,ii)'*B(:,k))*H(:,k);
    end
end

mEQMI_n(kk,:) = mean((X - Yhat).^2);
mEQMI_nt(kk,:) = mean((X - Yhat).^2,2)';

end


end

