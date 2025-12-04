function [ vd0p1_boot ] = DimEst_wavestrap2Residuals_teste( Y, Nboot, B, p, N, wname)
% Esse código faz a mesma coisa que é feita no código DimEst_wavestrap.m, a
% diferença que nesse código é usado o escopo do código
% DimEst_wavestrap2Residuals.m.

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
Aboot = zeros(J,n);

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

Aths = zeros(J,n);
%Aboot = zeros(J,n,Nboot);
ind_ths = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e Lw
    % o número de coeficientes em cada nível da decomposição
    [A_res(:,ii),~] = wavedec(Epshat(:,ii),N,wname);
    [A_H0(:,ii),Lw] = wavedec(Yhat(:,ii),N,wname);
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
end


% vetor para armazenar os (d0+1)-ésimos maiores autovalores das matrizes
% D formadas com os coeficientes bootstraps
for ii=1:Nboot
    
    Aboot = zeros(J,n);
    aux = zeros(J,n);
    
    for jj=1:n
        ts = unidrnd(n);
        %ts = jj;
        % vetor com coeficientes limiarizados cujos zeros serão
        % substituídos pelos coeficientes re-padronizados
        Aboot(:,jj) = A_H0(:,jj) + A_res(:,ts);
        aux(:,jj) = A_res(:,ts);        
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

