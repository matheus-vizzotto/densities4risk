function [ vd0p1_boot ] = DimEst_wavestrap_thresh( A, NREP, B, p, N, wname)
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
vd0p1_boot = zeros(NREP,1);
Aboot = zeros(J,n);

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B;

Bthr = B;
for ii=1:d0
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    Bthr(:,ii) = wden(B(:,ii),'sqtwolog','h','mln',N,wname);
end

for jj=1:NREP

% obtendo diretamente os coeficientes do funcional observado bootstrap
for ii=1:n
    ts = unidrnd(n);
    Aboot(:,ii) = A(:,ts) + Bthr*(mEta(ii,:) - mEta(ts,:))';
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
% reordenando os autovalores em ordem decrescente
[Lboot,~] = sort(diag(Lboot),'descend');
Lboot = diag(Lboot);

vd0p1_boot(jj) = Lboot(d0+1,d0+1);

end

end

