function [ vd0p1_boot ] = DimEst_Bathia_du( Y, NREP, mEigFunc, p, du)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo usando o método proposto por Bathia et al (2010).
% Y é a matriz de funções observadas, mEigFunc é a matriz de autofunções
% calculadas usando a equação (2.13) de Bathia et al (2010),
% NREP é o número de réplicas bootstrap e p é lag máximo. 

[~,d0] = size(mEigFunc);
[nt,n] = size(Y);
vd0p1_boot = zeros(NREP,1);
Yboot = zeros(nt,n);

mu_Y = mean(Y,2);
C = Y - mu_Y*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = du*C'*mEigFunc;

for jj=1:NREP
    % obtendo diretamente os funcionais bootstrap
    for ii=1:n
        ts = unidrnd(n);
        Yboot(:,ii) = Y(:,ts) + mEigFunc*(mEta(ii,:) - mEta(ts,:))';
    end

    mu_Yboot = mean(Yboot,2);
    C = Yboot - mu_Yboot*ones(1,n);

    C1 = C(:,1:(n-p));
    D1 = zeros(n-p,n-p);
    for k=1:p
        D1 = D1 + du*C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
    end
    Dboot = D1*(du*(C1'*C1))/((n-p)^2);

    [~,Lboot] = eig(Dboot);
    % reordenando os autovalores em ordem decrescente
    [Lboot,~] = sort(diag(Lboot),'descend');
    Lboot = diag(Lboot);
    vd0p1_boot(jj) = Lboot(d0+1,d0+1);

end

end

