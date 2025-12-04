function [ y ] = FuncLinCombGramSchm( fh_Ycent, EigVecK, combVec, d0, x)
% Essa função calcula uma combinação linear das d0 autofunções
% ortonormalizadas do operador usando os coeficientes em combVec. O escopo
% desse código é similar à FuncGramSchm.m, mas retorna uma matriz dessa
% combinação avaliada em x.
% fh_Ycent é uma célula com os funcionais observados centrados e EigVecK é
% uma matriz com os autovetores da matriz de produtos internos

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

y = 0;
for ii=1:d0
    % etapa que deixa as autofunções com norma 1
    innprod = InnerProdFunc(fh_OrtNormEigFunc{ii},fh_OrtNormEigFunc{ii},0,1);
    fh_OrtNormEigFunc{ii} = @(z) fh_OrtNormEigFunc{ii}(z)/sqrt(innprod);
    % formando a combinação linear das autofunções
    y = y + combVec(ii)*fh_OrtNormEigFunc{ii}(x);
end

end

