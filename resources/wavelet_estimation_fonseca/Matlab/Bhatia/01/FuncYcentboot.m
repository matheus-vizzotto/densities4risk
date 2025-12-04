function [ mat_Ycentboot ] = FuncYcentboot( fh_Y, fh_OrtNormEigFunc, mEta, indboot, k, x)
% Essa função retorna o funcional centrado bootstrap avaliado em x. fh_Y é
% uma célula contendo os funcionais observados, fh_OrtNormEigFunc é uma
% célula com as autofunções, mEta é a matriz com os coeficientes usados com
% as autofunções para formar as estimativas do funcional.
% O indboot é um vetor contendo os índices sorteados para formar o resíduo
% reamostrado, enquanto k é o índice do funcional centrado bootstrap que se
% deseja avaliar em x.

% tamanho da amostra
n = length(fh_Y);
% número de autofunções usadas (dimensão do funcinal sob H0)
d0 = length(fh_OrtNormEigFunc);

% array para armazenar as d0 autofunções avaliadas em x. Essa etapa permite
% realizar as próximas etapas matricialmente.
mat_OrtNormEigFunc = zeros([size(x) d0]);
for ii=1:d0
    mat_OrtNormEigFunc(:,:,ii) = fh_OrtNormEigFunc{ii}(x);
end

% gerando uma amostra bootstrap
% a diferença dos vetores mEta é referente aos termos do funcional sob H0 e
% do resíduo bootstrap. Produto de Kronecker é usado para criar um array
% com os coeficientes nesse vetor, que é multiplicado termo a termo com as
% matrizes das autofunções. A combinação termina somando as matrizes do
% array. 
fh_Yboot = zeros([size(x) n]);
for ii=1:n
    fh_Yboot(:,:,ii) = sum(mat_OrtNormEigFunc.*reshape(kron((mEta(ii,:) - mEta(indboot(ii),:))',ones(size(x)))',[size(x) d0]),3) ...
        + fh_Y{indboot(ii)}(x);
end

% funcional observado centrado
mat_Ycentboot = fh_Yboot(:,:,k) - mean(fh_Yboot,3);

% isso era o que era feito antes, mas demorava por que as autofunções eram
% analisadas em cada réplicada, o que não era preciso
% fh_Yboot = zeros([size(x) n]);
% for ii=1:n
%     mat_Yboot(:,:,ii) = LinCombFunc(fh_OrtNormEigFunc,(mEta(ii,:) - mEta(indboot(ii),:))',x) ...
%         + fh_Y{indboot(ii)}(x);
% end
% fh_Ycentboot = zeros(size(x));
% for ii=1:n
%     mat_Ycentboot = fh_Yboot(:,:,ii) - mean(fh_Yboot,3);
% end

end

