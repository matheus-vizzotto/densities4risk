function [Aboot] = ResampWavestrap(A, Lw, Nboot, N, wname)
% função para gerar os coeficientes de ondaletas de uma amostra
% bootstrap que foi gerada através do processo que usa limiarização
% para identificar quais coeficientes são ruído ou não na amostra
% original. O argumento A é a matriz dos coeficientes da amostra
% original, Lw é retornado da função wavedec enquanto N e wname
% são argumentos de wavedec.

[J,n] = size(A);

Aths = zeros(J,n);
Aboot = zeros(J,n,Nboot);
ind_ths = zeros(J,n);

medj = zeros(N,n); sigj = zeros(N,n); numj = zeros(N,n);
r = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % realizando a limiarização rígida ('h') com limiar universal
    % ('sqtwolog') e nível a nível ('mln')
    [~,Aths(:,ii),~] = wden(A(:,ii),Lw ,'sqtwolog','h','mln',N,wname);
    % armazenando os indicadores dos coeficientes limiarizados
    ind_ths(:,ii) = (A(:,ii) ~= Aths(:,ii));
    
    % soma acumulada do número de coeficientes em cada nível
    cumsum_Lw = cumsum(Lw);
    cont = 1;  % contador
    for k=1:N
        % média dos coeficientes limiarizados no nível k
        medj(k,ii) = mean(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        % desvio padrão dos coeficientes limiarizados no nível k
        sigj(k,ii) = std(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        % número de coeficientes limiarizados no nível k
        numj(k,ii) = sum(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii));
        % padronizando os coeficientes limiarizados do nível k, para que
        % eles tenham média zero e variância um
        r(cont:sum(numj(:,ii)), ii) = (A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii) - medj(k, ii))/sigj(k, ii);
        cont = sum(numj(:,ii)) + 1;  % contador
    end
    
    % soma acumulada do número de coeficientes limiarizados em cada nível
    cumsum_numj = cumsum(numj(:,ii));
    for bb=1:Nboot
        % sorteando os coeficientes limiarizados padronizados
        rs = r(unidrnd(sum(ind_ths(:,ii)),sum(ind_ths(:,ii)),1), ii);
        % vetor com coeficientes limiarizados cujos zeros serão
        % substituídos pelos coeficientes re-padronizados
        Aboot(:,ii,bb) = Aths(:,ii);
        cont = 1;  % contador
        for k=1:N
            % substituindo os zeros pelos coeficientes re-padronizados pela
            % média e desvio padrão dos coeficientes limiarizados do nível
            % k
            Aboot(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k),ii,bb) = rs(cont:cumsum_numj(k))*sigj(k,ii) + medj(k,ii);
            cont = cumsum_numj(k) + 1;
        end
    end

end


end

