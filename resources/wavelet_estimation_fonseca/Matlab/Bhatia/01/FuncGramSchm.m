function [ fh_OrtNormEigFunc ] = FuncGramSchm( fh_Ycent, EigVecK, d0)
% função para calcular as autofunções do operador analisado. fh_Ycent é uma
% célula armazenando os funcionais observados centrados, EigVecK é uma
% matriz contendo os autovetores (pelo menos os d0 primeiros) da matriz de
% produtos internos \hat{K}. As d0 primeiras autofunções serão retornadas.

% célula para armazenar autofunções do operador (ainda não ortonomalizadas)
fh_EigFunc = cell(d0,1);
% célula para armazenar as autofunções ortonormalizadas
fh_OrtNormEigFunc = cell(d0,1);
for ii=1:d0
    fh_EigFunc{ii} = @(z) LinCombFunc(fh_Ycent,EigVecK(:,ii),z);
    fh_OrtNormEigFunc{ii} = @(z) 0;
end

% aplicando o método de Gram-Schmidt nas autofunções
fh_OrtNormEigFunc{1} = @(z) fh_EigFunc{1}(z);
for ii=2:d0
    % etapa de ortogonalização
    vecGS = zeros(d0,1);
    for jj=1:(ii-1)
        vecGS(jj) = InnerProdFunc(fh_OrtNormEigFunc{jj},fh_EigFunc{ii},0,1)/InnerProdFunc(fh_OrtNormEigFunc{jj},fh_OrtNormEigFunc{jj},0,1);
    end
    fh_OrtNormEigFunc{ii} = @(z) fh_EigFunc{ii}(z) - LinCombFunc(fh_OrtNormEigFunc,vecGS,z);
end

for ii=1:d0
    % etapa que deixa as autofunções com norma 1
    innprod = InnerProdFunc(fh_OrtNormEigFunc{ii},fh_OrtNormEigFunc{ii},0,1);
    fh_OrtNormEigFunc{ii} = @(z) fh_OrtNormEigFunc{ii}(z)/sqrt(innprod);
end

end

