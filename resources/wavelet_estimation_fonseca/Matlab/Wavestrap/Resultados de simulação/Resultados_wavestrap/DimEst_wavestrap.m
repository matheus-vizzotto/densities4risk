function [ vd0p1_boot ] = DimEst_wavestrap( A, NREP, B, p)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, mas diferentemente do que é feito em DimEst_boot.m,
% aqui aproveitamos as decomposições obtidas para realizar a reamostragem
% nos termos formados com os coeficientes, o que deixa o código mais rápido.
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. 

[J,d0] = size(B);
[~,n] = size(A);
vd0p1_boot = zeros(NREP,1);
Aboot = zeros(J,n);

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B;

for jj=1:NREP

% obtendo diretamente os coeficientes do funcional observado bootstrap
for ii=1:n
    ts = unidrnd(n);
    Aboot(:,ii) = A(:,ts) + B*(mEta(ii,:) - mEta(ts,:))';
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

