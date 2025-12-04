function [ vd0p1_boot ] = DimEst_wavestrap3( A, Lw, Nboot, B, p, N, wname)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, as decomposições obtidas são aproveitadas na
% reamostragem, que teoricamente é feita sobre o resíduo obtido após ser 
% feita uma limiarização sobre os coeficientes das decomposições da média
% do processo e das autofunções.
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. 

[J,d0] = size(B);
[~,n] = size(A);
vd0p1_boot = zeros(Nboot,1);
Aboot = zeros(J,n);

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B;

Aeps = zeros(J,n);
Ahat = zeros(J,n);
for ii=1:n
    Ahat(:,ii) = mu_A + B*mEta(ii,:)';
    Aeps(:,ii) = A(:,ii) - Ahat(:,ii);
end
Aeps_boot = ResampWavestrap(Aeps, Lw, Nboot, N, wname);

for jj=1:Nboot
    Aboot = Ahat + Aeps_boot(:,:,jj);    
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

