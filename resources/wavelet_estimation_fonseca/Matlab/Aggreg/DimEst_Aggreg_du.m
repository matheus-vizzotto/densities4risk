function [ vd0p1_boot ] = DimEst_Aggreg_du( Y, NREP, mEigFunc, p, vw, du)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo usando o método proposto por Bathia et al (2010).
% Y é a matriz de funções observadas, mEigFunc é a matriz de autofunções
% calculadas usando a equação (2.13) de Bathia et al (2010),
% NREP é o número de réplicas bootstrap e p é lag máximo. 

[~,d0] = size(mEigFunc);
[nt,n] = size(Y);
vd0p1_boot = zeros(NREP,1);
Yboot = zeros(nt,n);
delta = length(vw); % número de pesos usados

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
    
    % funcional observado centrado
    mat_G = du*(C'*C);
    % matriz com os produtos internos dos funcionais centrados multiplicados
    % pelos pesos
    Zk = zeros(n-p-delta+1,n-p-delta+1,p-delta+1);
    for k=delta:(p-delta+1)
        mat_Z = FuncMatProdYcentAggregVet( C, vw, p, k);
        Zk(:,:,k) = du*mat_Z;
    end
    
    Dboot = sum(Zk,3)*mat_G(1:(n-p-delta+1),1:(n-p-delta+1))/((n-p-delta+1)^2);
    
    [~,Lboot] = eig(Dboot);
    % reordenando os autovalores em ordem decrescente
    [Lboot,~] = sort(diag(Lboot),'descend');
    Lboot = diag(Lboot);
    vd0p1_boot(jj) = Lboot(d0+1,d0+1);

end

end

