function [ mat_Zboot ] = FuncMatProdYcentAggregBoot( fh_Y, fh_OrtNormEigFunc, mEta, indboot, w, p, k, x)
% Essa função retorna uma matriz (n-p-delta+1)x(n-p-delta+1) com o produto para diferentes tempos
% dos funcionais observados bootstrap centrados, multiplicados por pesos e avaliados em x, que deve
% ser um escalar. 
% fh_Y é uma célula contendo os funcionais observados, fh_OrtNormEigFunc é uma
% célula com as autofunções, mEta é a matriz com os coeficientes usados com
% as autofunções para formar as estimativas do funcional.
% O indboot é um vetor contendo os índices sorteados para formar o resíduo
% reamostrado
% w é o vetor de pesos usado na agregação dos dados, p é o lag máximo usado
% na análise  e k é o lag utilizado nos funcionais ponderados.

% tamanho da amostra
n = length(fh_Y);
% número de autofunções usadas (dimensão do funcinal sob H0)
d0 = length(fh_OrtNormEigFunc);
% número de pesos no vetor w
delta = length(w);

% vetor para armazenar as d0 autofunções avaliadas em x.
mat_OrtNormEigFunc = zeros(d0,1);
for ii=1:d0
    mat_OrtNormEigFunc(ii) = fh_OrtNormEigFunc{ii}(x);
end

% gerando uma amostra bootstrap
% a diferença dos vetores mEta é referente aos termos do funcional sob H0 e
% do resíduo bootstrap. Esse vetor diferença contém os coeficientes da
% combinação linear das autofunções avaliadas em x, mat_OrtNormEigFunc
fh_Yboot = zeros(n,1);
for ii=1:n
    fh_Yboot(ii) = (mEta(ii,:) - mEta(indboot(ii),:))*mat_OrtNormEigFunc ...
        + fh_Y{indboot(ii)}(x);
end

% vetor dos funcionais observados centrados
vec_Ycentboot = fh_Yboot - mean(fh_Yboot);

% obtendo os funcionais observados ponderados
mY = zeros(delta,delta);
vZ = zeros(n-p-delta+1,1);
for jj = 1:(n-p-delta+1)
    for ii = 0:(delta-1)
        mY(ii+1,:) = vec_Ycentboot((jj+k-ii):(jj+k+delta-1-ii))';
    end
    vZ(jj) = w'*mY*w;
end

% matriz com o produto 2 a 2 entre todos os elementos do vetor
mat_Zboot = vZ*vZ';

end

