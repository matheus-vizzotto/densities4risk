% Como criar criar vetor de funções e acessar somente uma delas
% essas funções são conhecidas como anonymous functions
y = @(x) [2*x; 3*x; 5*x];
% função que retorna a somente uma das funções avaliada no argumento
indexat = @(expr, index) expr(index);
indexat(y(1), 3)

% tentativa de criar um vetor de funções
a = [1,2,3];
temp=cell(1,size(a,2));
temp{1,1}=1; 
for k=2:size(a,2)
temp{k}=sprintf('@(x) x.^%d',a(k));
end
aux = @(x) temp;

% Criando um célula de funções. Cada elemento dessa célula vai uma função
fh{1}=@(x) x.*x;
fh{2}=@(x) x.*x.*x;
fh{3}=@(x) x.*x.*x.*x;
fh{4}=@(x) sin(x);
fh{5}=@(x) cos(x);

% o cellfun converte a célula para uma função, que quando é avalida em z
% vai soltar os valores de cada elemento em z. O false indica que os
% valores estarão armazenados em outra célula
foo = @(z) cellfun( @(f,x) f(x), fh, repmat({z},size(fh)), 'uni', false);
foo(2)
% Faz o mesmo que a função anterior, mas o true indica que os valores do
% output serão escalares, formando um vetor
foo = @(z) cellfun( @(f,x) f(x), fh, repmat({z},size(fh)), 'uni', true);
foo(2)

% tentativa de criar uma matriz de funções usando células
fh = cell(2,2);
fh{1,1}=@(x) x.*x;
fh{1,2}=@(x) x.*x.*x;
fh{2,1}=@(x) x.*x.*x.*x;
fh{2,2}=@(x) sin(x);
foo2 = @(z) cellfun( @(f,x) f(x), fh, repmat({z},size(fh)), 'uni', true);
foo2(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cenário de simulação usado por Bathia et al (2010)

n = 50; % tamanho amostral
d = 2; % parâmetro de dimensão
p = 5;

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4
model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);
model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);

xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);

% célula para armazenar as funções não observada
fh_X = cell(n,1);
for ii=1:n
%    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u);
    fh_X{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x);
end

% function_handle que retorna um vetor com cada funcional observado
% avaliado em z
X = @(z) cellfun( @(f,x) f(x), fh_X, repmat({z},size(fh_X)), 'uni', true);

% função que retorna o elemento (i,j) da função
indexat = @(expr, i, j) expr(i,j);

% gerando a função ruído a ser adicionado na função não observada
fh_Eps = cell(n,1);
for ii=1:n
    vZ = normrnd(0,1,10,1);
    fh_Eps{ii} = @(x) sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
        + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
        + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
        + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
        + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
end

% função que retorna um vetor com os valores do resíduo
Eps = @(z) cellfun( @(f,x) f(x), fh_Eps, repmat({z},size(fh_Eps)), 'uni', true);

% célula que armazenará o funcional observado
fh_Y = cell(n,1);
% aqui o funcional observado já será feito usando os termos da função não
% observada, sem usar X ou fh_X
for ii=1:n
    vZ = normrnd(0,1,10,1);
    fh_Y{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x)...
        + sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
        + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
        + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
        + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
        + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
end


%%%%%%%%%%%%%

% primeira tentativa de criar o funcional observado usando function_handle
Y = @(z) cellfun( @(f,x) f(x), fh_Y, repmat({z},size(fh_Y)), 'uni', true);
%Y = @(z) cellfun( @(f,x) f(x), fh_Y, repmat({z},size(fh_Y)),'UniformOutput',false);

% função média
Ybar = @(z) sum(Y(z))/n;

Kstar = zeros(n-p,n-p);
for k=1:p
    Gk = zeros(n-p,n-p);
    for tt=1:(n-p)
         for ss=tt:(n-p)
             % função para pegar somente o índice tt+k da função Y
             Yt = @(z) indexat(Y(z),tt+k,1) - Ybar(z);
             % vetorizando a função Yt para poder calcular a integral
             % numérica no InnerProdFunc
             Yt = @(z) VecFunc(Yt,z);
             Ys = @(z) indexat(Y(z),ss+k,1) - Ybar(z);
             Ys = @(z) VecFunc(Ys,z);
             innprodGk = InnerProdFunc(Yt,Ys,0,1);
             Gk(tt,ss) = innprodGk;
             Gk(ss,tt) = innprodGk;
         end
    end
    
    Kstar = Kstar + Gk;
end

G0 = zeros(n-p,n-p);
for tt=1:(n-p)
     for ss=tt:(n-p)
         Yt = @(z) indexat(Y(z),tt,1) - Ybar(z);
         Yt = @(z) VecFunc(Yt,z);
         Ys = @(z) indexat(Y(z),ss,1) - Ybar(z);
         Ys = @(z) VecFunc(Ys,z);
         innprodGk = InnerProdFunc(Yt,Ys,0.001,.999);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
     end
end

Kstar = Kstar*G0/((n-p)^2);

%%%%%%%%%%%%%

tt = 1; ss = 1;
tic;
         Yt = @(z) indexat(Y(z),tt,1) - Ybar(z);
         Yt = @(z) VecFunc(Yt,z);
         Ys = @(z) indexat(Y(z),ss,1) - Ybar(z);
         Ys = @(z) VecFunc(Ys,z);
         innprodGk = InnerProdFunc(Yt,Ys,0.001,.999);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
toc;

% tentativa usando as anonymous functions armazenadas na célula fh_Y.
% Dentro das chaves {.} está o índice da função usada, avaliada em (.)
% A função média é calculada com a função FuncYbar
tic;
         Yt = @(z) fh_Y{tt}(z) - FuncYbar(fh_Y,n,z);
         Ys = @(z) fh_Y{ss}(z) - FuncYbar(fh_Y,n,z);
         innprodGk = InnerProdFunc(Yt,Ys,0.001,.999);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
toc;

%%%%%%%%%%%%%%%%

% célula para o funcional observado centrado
fh_Ycent = cell(n,1);
for ii=1:n
    fh_Ycent{ii} = @(z) fh_Y{ii}(z) - FuncYbar(fh_Y,n,z);
end

tt=1; ss=1;
G0 = zeros(n-p,n-p);
tic;
    % agora o produto interno é feito direto com o funcional centrado
    innprodGk = InnerProdFunc(fh_Ycent{tt},fh_Ycent{ss},0,1);
    G0(tt,ss) = innprodGk;
    G0(ss,tt) = innprodGk;
toc;

Kstar = zeros(n-p,n-p);
for k=1:p
    Gk = zeros(n-p,n-p);
    for tt=1:(n-p)
         for ss=tt:(n-p)
             innprodGk = InnerProdFunc(fh_Ycent{tt+k},fh_Ycent{ss+k},0,1);
             Gk(tt,ss) = innprodGk;
             Gk(ss,tt) = innprodGk;
         end
    end
    
    Kstar = Kstar + Gk;
end

G0 = zeros(n-p,n-p);
for tt=1:(n-p)
     for ss=tt:(n-p)
         innprodGk = InnerProdFunc(fh_Ycent{tt},fh_Ycent{ss},0,1);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
     end
end

% matriz usada para obter os autovalores das autofunções
Kstar = Kstar*G0/((n-p)^2);
% obtendo autovalores e autovetores de Kstar
[EigVecK,EigValK] = eig(Kstar);

% calculando as d0 primeiras autofunções do operador analisado
d0 = 2;
fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, d0);

% matriz com as va's que são os coeficientes das autofunções
mEta = zeros(n,d0);
for ii=1:n
   for jj=1:d0
       mEta(ii,jj) = InnerProdFunc(fh_Ycent{ii},fh_OrtNormEigFunc{jj},0,1);
   end
end


% gerando uma amostra bootstrap
fh_Ystar = cell(n,1);
for ii=1:n
    k = unidrnd(n);
    % o primeiro termo termo é referente o funcional definido com dimensão
    % d0, enquanto os dois últimos termos são referentes ao resíduo
    % amostrado
    fh_Ystar{ii} = @(z) LinCombFunc(fh_OrtNormEigFunc,mEta(ii,:)',z) ...
        + fh_Y{k}(z) - LinCombFunc(fh_OrtNormEigFunc,mEta(k,:)',z);
end


%%%%%%%%%%%%%%%%%

% testando os tempos de exucução de algumas funções

%FuncYcent(fh_Y,1,n,.5)
%fh_Ycent{1}(.5)

auxtest = rand(100,100);
tic;
%fh_Y{1}(auxtest);
fh_Ycent{1}(auxtest);
%FuncYcent(fh_Y,1,n,auxtest);
toc;

auxtest = rand(100,100);
tic;
%fh_OrtNormEigFunc{1}(auxtest);
%LinCombFunc(fh_OrtNormEigFunc,mEta(1,:)',auxtest);
FuncLinCombGramSchm(fh_Ycent, EigVecK, mEta(1,:)', d0, auxtest);
toc;

tic;
fh_Ycentboot2{1}(aux)
%FuncYcentboot( fh_Y, fh_OrtNormEigFunc, mEta, indboot, 1, auxtest);
toc;


