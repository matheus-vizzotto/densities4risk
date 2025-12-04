function [ vd0p1_boot ] = DimEst_wavestrap_ResidualsDimPlus1( A, NREP, B, p)
% Código que faz a mesma coisa que o DimEst_wavestrap.m, com a única
% diferença que os resíduos são centrados aqui para terem média zero

[J,d0] = size(B);
d0 = d0 - 1;
[~,n] = size(A);
vd0p1_boot = zeros(NREP,1);
A_H0 = zeros(J,n);
A_res = zeros(J,n);
Aboot = zeros(J,n);

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B(:,1:d0);
mEta2res = C'*B;

% obtendo diretamente os coeficientes do funcional sob H0 bootstrap e do
% resíduo gerado
for ii=1:n
    A_H0(:,ii) = mu_A + B(:,1:d0)*(mEta(ii,:))';
    A_res(:,ii) = A(:,ii) - mu_A - B*(mEta2res(ii,:))';
end

% deixando os resíduos com média zero ao longo dos dias
mu_A_res = mean(A_res,2);
A_res = A_res - mu_A_res*ones(1,n);


for jj=1:NREP
    
    % obtendo diretamente os coeficientes do funcional observado bootstrap
    for ii=1:n
        ts = unidrnd(n);
        Aboot(:,ii) = A_H0(:,ii) + A_res(:,ts);
    end

    mu_Aboot = mean(Aboot,2);
    C = Aboot - mu_Aboot*ones(1,n);

    C1 = C(:,[1:(n-p)]);
    D1 = zeros(n-p,n-p);
    for k=1:p
        D1 = D1 + C(:,[(k+1):(n-p+k)])'*C(:,[(k+1):(n-p+k)]);
    end
    Dboot = C1*D1*C1'/((n-p)^2);

    [~,Lboot] = eig(Dboot);
    vd0p1_boot(jj) = Lboot(d0+1,d0+1);
end

end

