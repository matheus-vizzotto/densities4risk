% Neste código simulamos uma amostra de uma série temporal de funções como
% descrito na estudo numérico de Bathia. Depois, geramos versões bootstrap
% da amostra observada Y usando o seguinte procedimento: os coeficientes
% das decomposições de cada dia são limiarizados com hard thresholding, onde
% identificamos os coef. que são ruído; esses coeficientes são repadronizados
% em cada nível para terem média zero e desvio padrão um; a reamostragem é
% feita sobre esses rezíduos padronizados, que são reescalados pela média e std
% dos níveis e colocados no lugar do coef. limiarizados. Com essas amostras
% bootstrap obtemos vetores de autovalores calculados, e então contruímos
% intervalos de confianças, que tem seus comprimentos comparados

n = 100; % tamanho amostral
d = 4; % parâmetro de dimensão
nt = 256; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
N = 3;  % nível usado na função wavedec
wname = 'db2';
Nboot = 100; % número de réplicas bootstrap
p = 5;  % lag máximo utilizado

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4, quando d=4 são -.775,.65,-.525 e .4,
% e quado d=6 são -.8167,.7333,-.65,.5667,-.4833,.4
model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);
model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);
model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);
model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
%model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);
%model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);

%rng('default')
xi1 = simulate(model_xi1,n);
xi2 = simulate(model_xi2,n);
xi3 = simulate(model_xi3,n);
xi4 = simulate(model_xi4,n);
%xi5 = simulate(model_xi5,n);
%xi6 = simulate(model_xi6,n);

% gerarndo a função não observada
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u)...
       + xi3(ii)*sqrt(2)*cos(3*pi*u) + xi4(ii)*sqrt(2)*cos(4*pi*u);
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
J = 263;
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);
Aths = zeros(J,n);
Aboot = zeros(J,n,Nboot);
ind_ths = zeros(J,n);

medj = zeros(N,n); sigj = zeros(N,n); numj = zeros(N,n);
r = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e Lw
    % o número de coeficientes em cada nível da decomposição
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    % obtendo valores defaut para realizar a limiarização
    %[thr,sorh,keepapp] = ddencmp('den','wv',Y(:,ii));
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,Aths(:,ii),~] = wden(Y(:,ii),'sqtwolog','h','mln',N,wname);
    % armazenando os indicadores dos coeficientes limiarizados
    ind_ths(:,ii) = (A(:,ii) ~= Aths(:,ii));
    
    % soma acumulada do número de coeficientes em cada nível
    cumsum_Lw = cumsum(Lw);
    cont = 1;  % contador
    for k=1:N
        % média dos coeficientes limiarizados no nível k
        medj(k,ii) = mean(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        % desvio padrão dos coeficientes limiarizados no nível k
        sigj(k,ii) = std(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        % número de coeficientes limiarizados no nível k
        numj(k,ii) = sum(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii));
        % padronizando os coeficientes limiarizados do nível k, para que
        % eles tenham média zero e variância um
        r(cont:sum(numj(:,ii)), ii) = (A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii) - medj(k, ii))/sigj(k, ii);
        cont = sum(numj(:,ii)) + 1;  % contador
    end
    
    % soma acumulada do número de coeficientes limiarizados em cada nível
    cumsum_numj = cumsum(numj(:,ii));
    for bb=1:Nboot
        % sorteando os coeficientes limiarizados padronizados
        rs = r(unidrnd(sum(ind_ths(:,ii)),sum(ind_ths(:,ii)),1), ii);
        % vetor com coeficientes limiarizados cujos zeros serão
        % substituídos pelos coeficientes re-padronizados
        Aboot(:,ii,bb) = Aths(:,ii);
        cont = 1;  % contador
        for k=1:N
            % substituindo os zeros pelos coeficientes re-padronizados pela
            % média e desvio padrão dos coeficientes limiarizados do nível
            % k
            Aboot(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k),ii,bb) = rs(cont:cumsum_numj(k))*sigj(k,ii) + medj(k,ii);
            cont = cumsum_numj(k) + 1;
        end
    end

end

% vetor para armazenar nas linhas os autovalores da matriz D formada com os
% coeficientes da amostra bootstrap formados anteriormente
vd_boot = zeros(Nboot,J);
for ii=1:Nboot
    mu_A = mean(Aboot(:,:,ii),2);
    C = Aboot(:,:,ii) - mu_A*ones(1,n);
    C1 = C(:,1:(n-p));
    D1 = zeros(n-p,n-p);
    for k=1:p
        D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
    end
    D = C1*D1*C1'/((n-p)^2);
    [~,L] = eig(D);
    vd_boot(ii,:) = diag(L)';
end

%mu_A = mean(Aths,2);
%C = Aths - mu_A*ones(1,n);

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
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);

% após obter D, calculamos os autovetores e autovalores da mesma.
% as colunas de B terão os autovetores de D, que são os coeficientes de
% ondaletas das autofunções associadas aos autovalores na diagonal de L
[B,L] = eig(D);

vIC = quantile(vd_boot(:,:),.975) - quantile(vd_boot(:,:),.025);
plot(vIC(1:10)./vIC(2:11))

