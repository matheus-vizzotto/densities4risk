
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matriz de funções usando células
fh = cell(2,2);
fh{1,1}=@(x) x.*x;
fh{1,2}=@(x) x.*x.*x;
fh{2,1}=@(x) x.*x.*x.*x;
fh{2,2}=@(x) sin(x);
foo2 = @(z) cellfun( @(f,x) f(x), fh, repmat({z},size(fh)), 'uni', true);
foo2(2)
% integral das funções contidas na matriz. A opção 'ArrayValued' true
% permite realizar esse tipo de integral.
integral(foo2,0,1,'ArrayValued',true)

% comparando o tempo usado para formar a matriz Kstar da forma anterior e
% usando uma matriz de funções dos produtos dos funcionais centrados

tic;
G0 = zeros(n-p,n-p);
for tt=1:(n-p)
     for ss=tt:(n-p)
         innprodGk = InnerProdFunc(fh_Ycent{tt},fh_Ycent{ss},0,1);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
     end
end
toc;

tic;
G02 = integral(@(x) FuncMatProdYcent( fh_Y, x),0,1,'ArrayValued',true);
toc;

sum(sum(abs(G0 - G02(1:45,1:45)),2))

%%%%%%%

% Comparando o tempo de execução da forma anterior de obter a matriz mEta e
% uma forma usando uma matriz de funções do produto dos funcionais
% centrados e das autofunções

tic;
FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, .2 );
toc;

tic;
% matriz com as va's que são os coeficientes das autofunções
mEta = zeros(n,d0);
for ii=1:n
   for jj=1:d0
       mEta(ii,jj) = FuncInnerProd(fh_Ycent{ii},fh_OrtNormEigFunc{jj},0,1);
   end
end
toc;

tic;
mEta = integral(@(x) FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, x ),0,1,'ArrayValued',true);
toc;



