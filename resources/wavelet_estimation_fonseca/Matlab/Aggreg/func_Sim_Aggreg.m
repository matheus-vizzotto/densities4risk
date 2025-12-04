function [ vDimSel, mEQMI_n, mEQMI_nt, mEigval, mPvalues ] = func_Sim_Aggreg(n,d,NREP,Nboot,nt)
% Essa função é utilizada para simular a estimação da dimensão de uma série
% temporal de curvas seguindo como fizeram Bathia et al (2010). As entradas
% são o tamanho amostral (número de curvas observadas) n, a dimensão do
% processo d, o número de réplicas da simulação NREP, o número de réplicas
% bootstrap Nboot no teste da dimensão, o número de pontos nt para cada dia
% Aqui usamos o método de Bathia et al (2010)

% Parâmetros já definidos da simulação
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
du = u(2) - u(1);
p = 5;  % lag máximo utilizado
alpha = .05;  % nível de significância do teste bootstrap
vw = [0.1,0.1,0.8]'; % vetor de pesos
%vw = vw/norm(vw);
delta = length(vw); % número de pesos usados

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
[L,~] = sort(diag(L),'descend');
L = diag(L);

% 8 maiores autovalores de D
mEigval(kk,:) = diag(L([1:8],[1:8]));

% calculando as autofunções como Bathia et al (2010) descreve na equação
% (2.13)
mEigFunc = zeros(nt,10);
for ii=1:10
    mEigFunc(:,ii) = du*C(:,1:(n-p-delta+1))*B(:,ii);
end
mEigFunc = MatGramSchm(mEigFunc);

% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_Bathia
    d_boot = DimEst_Aggreg_du( Y, Nboot, mEigFunc(:,1:d0), p, vw, du);
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

Yhat = zeros(nt,n);  % matriz para armazenar os funcionais estimados
% recuperando os funcionais estimados
for ii=1:n
    Yhat(:,ii) = mu_Y;
    for k=1:(vDimSel(kk))
        Yhat(:,ii) = Yhat(:,ii) + (du*C(:,ii)'*mEigFunc(:,k))*mEigFunc(:,k);
    end
end

mEQMI_n(kk,:) = mean((X - Yhat).^2);
mEQMI_nt(kk,:) = mean((X - Yhat).^2,2)';

end


end

