function [ vd0p1_boot ] = DimEst_wavestrap2Residuals( Y, Nboot, B, p, N, wname)
% Nessa função é realizado um procedimento bootstrap para estimar a 
% dimensão do processo. Obtemos os coeficientes de ondaleta para Y e
% realizamos uma limiarização rígida, identificando quais destes coef. são
% ruído. A reamostragem é realizada sobre os coef. limiarizados do ruído
% (padronizados anteriormete), que são sorteados para substituir os zeros
% da matriz de coef. Aths limiarizados, nos respectivos dias e níveis.
% Depois esse ruído é somado à matriz de coef. do funcional sob H0.
% Y é a matriz de funcionais observados, B é a matriz cujas linhas são
% observações das funções que geram Y. N e wname são argumentos da função
% waavdec. Nboot é o número de réplicas bootstrap. p é lag máximo. 

[J,d0] = size(B);
[nt,n] = size(Y);
vd0p1_boot = zeros(Nboot,1);
H = zeros(nt,d0);
Yhat = zeros(nt,n);
Epshat = zeros(nt,n);
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);
A_H0 = zeros(J,n);
A_res = zeros(J,n);

for ii = 1:n
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
%    [~,A(:,ii),Lw] = wden(Y(:,ii),'sqtwolog','h','mln',N,'db2');
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
    Epshat(:,ii) = Y(:,ii) - Yhat(:,ii);
end

% matriz cujos vetores armazenarão os coeficientes limiarizados cujos, de
% valor zero usando limiar rígido. Esses zeros serão substituídos pelos
% coeficientes detectados como ruído re-padronizados.
Aths = zeros(J,n);
ind_ths = zeros(J,n);

medj = zeros(N,n); sigj = zeros(N,n); numj = zeros(N,n);
r = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e Lw
    % o número de coeficientes em cada nível da decomposição
    [A_res(:,ii),~] = wavedec(Epshat(:,ii),N,wname);
    [A_H0(:,ii),Lw] = wavedec(Yhat(:,ii),N,wname);
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,Aths(:,ii),~] = wden(Epshat(:,ii),'sqtwolog','h','mln',N,wname);
    % armazenando os indicadores dos coeficientes limiarizados
    ind_ths(:,ii) = (A_res(:,ii) ~= Aths(:,ii));
    
    % soma acumulada do número de coeficientes em cada nível
    cumsum_Lw = cumsum(Lw);
    cont = 1;  % contador
    for k=1:N
        % média dos coeficientes limiarizados no nível k
        medj(k,ii) = mean( A_res(find( ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii) )+cumsum_Lw(k), ii) );
        % desvio padrão dos coeficientes limiarizados no nível k
        sigj(k,ii) = std(A_res(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        % número de coeficientes limiarizados no nível k
        numj(k,ii) = sum(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii));
        % padronizando os coeficientes limiarizados do nível k, para que
        % eles tenham média zero e variância um
        r(cont:sum(numj(:,ii)), ii) = (A_res(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii) - medj(k, ii))/sigj(k, ii);
        cont = sum(numj(:,ii)) + 1;  % contador
    end

end


% vetor para armazenar os (d0+1)-ésimos maiores autovalores das matrizes
% D formadas com os coeficientes bootstraps
for ii=1:Nboot
    
    Aboot = zeros(J,n);
    
    for jj=1:n
        ts = unidrnd(n);
        %ts = jj;
        % soma acumulada do número de coeficientes limiarizados em cada nível
        cumsum_numj = cumsum(numj(:,ts));
        % sorteando os coeficientes limiarizados padronizados
        rs = r(unidrnd(sum(ind_ths(:,ts)),sum(ind_ths(:,ts)),1), ts);
        % se o índice corresponder a um coeficiente que foi limiarizado,
        % Aths vai ter zero nele, daí será somado uma versão re-padronizada
        % do resíduo sorteado, em rs, multiplicado do desv. pad. do nível
        % somada da média do nível, respectivamente. 
        Aboot(:,jj) = A_H0(:,jj) + Aths(:,ts);
        cont = 1;  % contador
        for k=1:N
            % substituindo os zeros pelos coeficientes re-padronizados pela
            % média e desvio padrão dos coeficientes limiarizados do nível
            % k
            Aboot(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ts))+cumsum_Lw(k),jj) = Aboot(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ts))+cumsum_Lw(k),jj) + rs(cont:cumsum_numj(k))*sigj(k,ts) + medj(k,ts);
            cont = cumsum_numj(k) + 1;
        end 
    end
    
    mu_A = mean(Aboot,2);
    C = Aboot - mu_A*ones(1,n);
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

