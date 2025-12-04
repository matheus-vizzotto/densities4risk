function [ vd0p1_boot ] = DimEst_wavestrap2( Y, Nboot, B, p, N, wname)
% Nessa função é realizado um procedimento bootstrap para estimar a 
% dimensão do processo. Obtemos os coeficientes de ondaleta para Y sob H0 e
% realizamos uma limiarização rígida, identificando quais destes coef. são
% ruído. A reamostragem é realizada sobre os coef. limiarizados
% (padronizados anteriormete), que são sorteados para substituir os zeros
% da matriz de coef. limiarizados, nos respectivos dias e níveis.
% Y é a matriz de funcionais observados, B é a matriz cujas linhas são
% observações das funções que geram Y. N e wname são argumentos da função
% waavdec. Nboot é o número de réplicas bootstrap. p é lag máximo. 

[J,d0] = size(B);
[nt,n] = size(Y);
vd0p1_boot = zeros(Nboot,1);
H = zeros(nt,d0);
Yhat = zeros(nt,n);
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);

for ii = 1:n
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
end
mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% tomando a transformada inversa da função média
mu_hat =  mean(Y,2);

for ii=1:d0
    H(:,ii) = waverec(B(:,ii),Lw,wname);
end

for ii=1:n
    Yhat(:,ii) = mu_hat;
    for k=1:d0
        Yhat(:,ii) = Yhat(:,ii) + (C(:,ii)'*B(:,k))*H(:,k);
    end
end


Aths = zeros(J,n);
Aboot = zeros(J,n,Nboot);
ind_ths = zeros(J,n);

medj = zeros(N,n); sigj = zeros(N,n); numj = zeros(N,n);
r = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e Lw
    % o número de coeficientes em cada nível da decomposição
    [A(:,ii),Lw] = wavedec(Yhat(:,ii),N,wname);
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,Aths(:,ii),~] = wden(Yhat(:,ii),'sqtwolog','h','mln',N,wname);
%    [~,Aths(:,ii),~] = wden(Yhat(:,ii),'sqtwolog','h','sln',N,wname);
    % armazenando os indicadores dos coeficientes limiarizados
    ind_ths(:,ii) = (A(:,ii) ~= Aths(:,ii));
    
    % soma acumulada do número de coeficientes em cada nível
    cumsum_Lw = cumsum(Lw);
    cont = 1;  % contador
    for k=1:N
        % média dos coeficientes limiarizados no nível k
        medj(k,ii) = mean( A(find( ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii) )+cumsum_Lw(k), ii) );
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

% vetor para armazenar os (d0+1)-ésimos maiores autovalores das matrizes
% D formadas com os coeficientes bootstraps
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
    vd0p1_boot(ii) = L(d0+1,d0+1);
end


end

