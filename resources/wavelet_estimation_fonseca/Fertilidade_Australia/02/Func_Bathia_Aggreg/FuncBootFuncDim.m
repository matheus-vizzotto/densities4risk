function [ vecd0 ] = FuncBootFuncDim( fh_Y, fh_OrtNormEigFunc, mEta, d0, p, Nboot )
% Função para obter uma amostra bootstrap dos autovalores da matriz K* de
% produtos internos. Considerando que se dejesa testar se o (d0+1)-ésimo
% autovalor é zero, serão gerados funcionais bootstrap com dimensão d0,
% para os quais serão calculados o (d0+1)-ésimo autovalor da K*
% correspondente, que serão armazenados em vecd0.
% fh_Y é uma célula contendo os funcionais observados, fh_OrtNormEigFunc é uma
% célula com as autofunções do operador (as d0 primeiras), mEta 
% é a matriz com os coeficientes usados com as autofunções para formar as
% estimativas do funcional. d0 é a dimensão do funcional criado sob H0, p é
% a defasagem máxima e Nboot é o número de réplicas bootstrap.
% É preciso que fh_OrtNormEigFunc tenha d0 elementos e que mEta tenha d0
% colunas.

% tamanho da amostra
n = length(fh_Y);
% vetor que armazenará a amostra bootstrap do (d0+1)-ésimo autovalor
vecd0 = zeros(Nboot,1);

% gerando uma amostra bootstrap
for jj=1:Nboot
    
    % sorteando os índices do resíduo bootstrap reamostrado
    indboot = unidrnd(n,n,1);
    
    % matriz com o produto dos funcionais observados bootstrap centrados
    mat_G = integral(@(x) FuncMatProdYcentboot( fh_Y, fh_OrtNormEigFunc, mEta, indboot, x),0,1,'ArrayValued',true);
    
    % formando a matriz com mesmos autovalores que o operador
    Gk = zeros(n-p,n-p,p);
    for k=1:p
        Gk(:,:,k) = mat_G((1+k):(n-p+k),(1+k):(n-p+k));
    end
    Kstar = sum(Gk,3)*mat_G(1:(n-p),1:(n-p))/((n-p)^2);
    
    [~,EigValK] = eig(Kstar);
    vecd0(jj) = EigValK(d0+1,d0+1);
end

end

