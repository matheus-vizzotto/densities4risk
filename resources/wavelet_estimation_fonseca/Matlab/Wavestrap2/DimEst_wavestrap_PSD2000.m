function [ vd0p1_boot ] = DimEst_wavestrap_PSD2000( A, Lw, NREP, B, p)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, de maneira similar ao que é feito em
% DimEst_wavestrap.m, mas utilizando também a técnica de wavestrap proposta
% por Persival, Sardy e Davidson (2000).
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. Lw é o vetor
% retornado pela função wavedec, contendo o número de coeficientes em cada
% nível da decomposição.

[J,d0] = size(B);
[~,n] = size(A);
vd0p1_boot = zeros(NREP,1);
% matriz para armazenar os coeficientes da decomposição dos funcionais sob
% a hipótese nula
A_H0 = zeros(J,n);
% matriz para armazenar os coeficientes da decomposição dos resíduos
A_res = zeros(J,n);

% número de componentes da decomposição (níveis de detalhe + 1 de aproximação)
levDec = length(Lw) - 1;
% vetor com os índices sorteados no wavestrap
ind_boot = zeros(J,1);
% os primeiros coeficientes são de aproximação e não são sorteados
ind_boot(1:Lw(1)) = 1:Lw(1); 

mu_A = mean(A,2);
C = A - mu_A*ones(1,n);
% matriz de etas (as v.a.) para cada autofunção, nas colunas, para
% diferentes dias, nas linhas
mEta = C'*B;

% obtendo diretamente os coeficientes do funcional sob H0 bootstrap e do
% resíduo gerado
for ii=1:n
    A_H0(:,ii) = mu_A + B*(mEta(ii,:))';
    A_res(:,ii) = A(:,ii) - mu_A - B*(mEta(ii,:))';
end

for jj=1:NREP
    
    % coeficientes simulados com wavestrap
    Aboot = zeros(J,n);
    for ii=1:n
        
        % laço para sortear os níveis de detalhes no wavestrap
        for k=2:levDec
            % posição do vetor ind_boot com último valor já definido
            cs_i = cumsum(Lw(1:(k-1))); cs_i = cs_i(end);
            % última posição do vetor ind_boot que será definida agora
            cs_f = cumsum(Lw(1:k)); cs_f = cs_f(end);
            % sorteando os elementos de {cs_i + 1,...,cd_f}
            ind_boot((cs_i+1):cs_f) = unidrnd(Lw(k),Lw(k),1) + cs_i;
        end
        % tempo sorteado cujo resíduo vai ser usado
        ts = unidrnd(n);
        % coeficientes do functional observado bootstrap, onde A_H0(:,ii)
        % são os coeficientes de \tilde{f}, e A_res dos resíduos bootstrap.
        % Aqui, ts é referente ao resíduo sorteado, dentre os n possíveis,
        % e ind_boot embaralha seus coeficientes dos níveis de detalhe.
        Aboot(:,ii) = A_H0(:,ii) + A_res(ind_boot,ts);
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

