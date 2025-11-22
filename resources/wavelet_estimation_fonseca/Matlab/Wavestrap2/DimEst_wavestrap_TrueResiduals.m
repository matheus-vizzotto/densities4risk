function [ vd0p1_boot ] = DimEst_wavestrap_TrueResiduals( A, Lw, NREP, B, p, nt, N, wname)
% Nesse código realizamos o teste de hipóteses da dimensão utilizando o
% esquema verdadeiro de geração de resíduos ao invés do procedimento
% bootstrap, para verificar se o teste funciona com esses resíduos.

u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
[J,d0] = size(B);
[~,n] = size(A);
vd0p1_boot = zeros(NREP,1);
A_H0 = zeros(J,n);
A_res = zeros(J,n);

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B;

% matriz que armazenará os resíduos a serem adicionados na função
mEps = zeros(nt,n);

% obtendo diretamente os coeficientes do funcional observado bootstrap
for ii=1:n
    A_H0(:,ii) = mu_A + B*(mEta(ii,:))';
end

for jj=1:NREP
    
    mEps = zeros(nt,n);
    for ii=1:n
        % gerando os resíduos
        for k=1:10
            mEps(:,ii) = mEps(:,ii) +...
                normrnd(0,1)*sqrt(2)*sin(pi*u*k)/(2^(k-1));
        end
        % obtendo a decomposição dos resíduos
        [A_res(:,ii),~] = wavedec(mEps(:,ii),N,wname);
    end
    
    % decomposição da função sob H0 observada
    Aboot = A_H0 + A_res;

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

