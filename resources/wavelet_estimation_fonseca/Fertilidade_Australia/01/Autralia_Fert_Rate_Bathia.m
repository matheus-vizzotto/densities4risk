
% Taxas anuais de fertilidade na Australia de 1921 a 2000
FertRate = csvread('Australia_Fert_Rate.csv');
% logaritmo das taxas
logFerRate = [FertRate(1,:);FertRate(2:10,1) log(FertRate(2:10,2:81))];


n = 80;  % número de funcionais
% variáveis para usar agregação de dados
vw = [.1,.3,.5]'; % vetor de pesos
%vw = vw/norm(vw);
delta = length(vw); % número de pesos usados
p = 5; % lag máximo usado

% adicionando caminho para chamar as funções necessárias
addpath(genpath('Func_Bathia_Aggreg'))

% célula com funcionais observados
fh_spl = cell(n,1);

% idades centradas no intervalo [0,1]
vAgeCent01 = (logFerRate(2:10,1) - 13)./(52 - 13);
% laço para obter a função suavizada da log-taxa em cada ano, sendo
% armazenada a função avaliada nos pontos do vetor u
for k=2:81
    fh_spl{k-1} = fit(logFerRate(2:10,1),logFerRate(2:10,k),'smoothingspline');
%    fh_spl{k-1} = fit(vAgeCent01,logFerRate(2:10,k),'smoothingspline');
    fh_spl{k-1} = @(x) vec2mat(fh_spl{k-1}(x),size(x,2));
end

% hold on
% fplot(fh_spl{1},[0 1])
% fplot(fh_spl{20},[0 1])
% fplot(fh_spl{40},[0 1])
% fplot(fh_spl{60},[0 1])
% fplot(fh_spl{80},[0 1])
% hold off

%%%%%%%%%%%%%%%%
% método com agregação de dados

tic;
% funcional observado centrado
mat_G = integral(@(x) FuncMatProdYcent( fh_spl, x),13,52,'ArrayValued',true);

% matriz com os produtos internos dos funcionais centrados multiplicados
% pelos pesos
Zk = zeros(n-p-delta+1,n-p-delta+1,p-delta+1);
for k=delta:(p-delta+1)
    mat_Z = integral(@(x) FuncMatProdYcentAggreg( fh_spl, vw, p, k, x),13,52,'ArrayValued',true);
    Zk(:,:,k) = mat_Z;
end
toc;

Kstar = sum(Zk,3)*mat_G(1:(n-p-delta+1),1:(n-p-delta+1))/((n-p-delta+1)^2);

[EigVecK,EigValK] = eig(Kstar);
% 10 maiores autovalores
diag(EigValK(1:10,1:10))'
% 5 maiores autovalores divididos pelo maior deles (x100)
round(100*real(diag(EigValK(1:5,1:5)))/norm(diag(EigValK)),4)'

%%%%%%%%%%%%%%%%%%%%%%
% método do Bathia et al (2010)

Gk_BYZ = zeros(n-p,n-p,p);
for k=1:p
    Gk_BYZ(:,:,k) = mat_G((1+k):(n-p+k),(1+k):(n-p+k));
end

Kstar_BYZ = sum(Gk_BYZ,3)*mat_G(1:(n-p),1:(n-p))/((n-p)^2);

[EigVecK_BYZ,EigValK_BYZ] = eig(Kstar_BYZ);
% 10 maiores autovalores
diag(EigValK_BYZ(1:10,1:10))'
% 5 maiores autovalores divididos pelo maior deles (x100)
round(100*real(diag(EigValK_BYZ(1:5,1:5)))/norm(diag(EigValK_BYZ)),4)'

%%%%%%%%%%%%%%%%%%%%%%

% célula para o funcional observado centrado
fh_Ycent = cell(n,1);
for ii=1:n
    fh_Ycent{ii} = @(z) fh_spl{ii}(z) - FuncYbar(fh_spl,n,z);
end

% calculando as d0 primeiras autofunções do operador analisado
dmax = 5;
fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK_BYZ(:,1:10), dmax, 13, 52);
mEta = integral(@(x) FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, x ),13,52,'ArrayValued',true);

rng(2018)
Nboot = 1; alpha = .05; vDimSel = 1;
checkDimSel = 0;
d0 = 1;
mPvalues = ones(1,dmax);
tic;
while d0<=dmax
    fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK_BYZ, d0, 13, 52);
    d_boot = FuncBootFuncDim(fh_Y, fh_OrtNormEigFunc, mEta(:,1:d0), d0, p, Nboot);
    mPvalues(d0) = sum(d_boot>EigValK(d0+1,d0+1))/Nboot;
    if (mPvalues(d0)>alpha)&&(checkDimSel==0)
        vDimSel = d0;
        checkDimSel = 1;
    end
    d0 = d0 + 1;
end
toc;
vDimSel
mPvalues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


