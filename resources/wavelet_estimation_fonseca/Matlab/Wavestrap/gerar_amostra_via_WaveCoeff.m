% Nesse código é feito o esquema para gerar amostras bootstrap do funcional
% observado Y, usando limiarização para identificar quais coeficientes são
% ruídos, e realizando a reamostragem nestes.


n = 1; % tamanho amostral
d = 2; % parâmetro de dimensão
nt = 64; % número de pontos para cada dia
u = linspace(0.01,.99,nt)'; % vetor ntx1 de pontos equiespaçados
N = 3;  % nível usado na função wavedec
wname = 'db2';

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4, quando d=4 são -.775,.65,-.525 e .4,
% e quado d=6 são -.8167,.7333,-.65,.5667,-.4833,.4
model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);
model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);

%rng('default')
xi1 = simulate(model_xi1,n);
xi2 = simulate(model_xi2,n);

% gerarndo a função não observada
X = zeros(nt,n);
X = xi1*sqrt(2)*cos(pi*u) + xi2*sqrt(2)*cos(2*pi*u);

% gerando o ruído a ser adicionado na função não observada
mEps = zeros(nt,n);
for jj=1:10
    mEps = mEps + normrnd(0,1)*sqrt(2)*sin(pi*u*jj)/(2^(jj-1));
end

% funcional observado
Y = X + mEps;

% número de coefientes da base de ondaletas usada
J = 263;
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);
Aths = zeros(J,n);

% D aqui irá receber os coeficiente de aproximação e L
% os coeficientes de detalhes
%[D,L] = dwt(Y(:,ii),'db2');
% os coeficientes são empilhados na i-ésima coluna de A
%A(:,ii) = [D',L']';
[A,Lw] = wavedec(Y,N,wname);
[thr,sorh,keepapp] = ddencmp('den','wv',Y);
[~,Aths,Lwths,~,~] = wdencmp('gbl',A,Lw,wname,1,thr,'h',keepapp);
ind_ths = (A~=Aths);

%%%%%%%%%%%%%%%%%%%%%%%%

medj = [mean(A(find(ind_ths(11:20))+10)), mean(A(find(ind_ths(21:38))+20)), mean(A(find(ind_ths(39:71))+38))];
sigj = [std(A(find(ind_ths(11:20))+10)), std(A(find(ind_ths(21:38))+20)), std(A(find(ind_ths(39:71))+38))];
numj = [sum(ind_ths(11:20)), sum(ind_ths(21:38)), sum(ind_ths(39:71))];

r = [(A(find(ind_ths(11:20))+10) - medj(1))/sigj(1); (A(find(ind_ths(21:38))+20) - medj(2))/sigj(2); ...
    (A(find(ind_ths(39:71))+38) - medj(3))/sigj(3)];
rs = r(unidrnd(length(r),length(r),1));

Aths2 = Aths;
cumsum_numj = cumsum(numj);
Aths2(find(ind_ths(11:20))+10) = rs(1:cumsum_numj(1))*sigj(1) + medj(1);
Aths2(find(ind_ths(21:38))+20) = rs((cumsum_numj(1)+1):cumsum_numj(2))*sigj(2) + medj(2);
Aths2(find(ind_ths(39:71))+38) = rs((cumsum_numj(2)+1):cumsum_numj(3))*sigj(3) + medj(3);

Yboot = waverec(Aths2,Lwths,wname);

%%%%%%%%%%%%%%%%%%%%%%%%

medj = zeros(N,1); sigj = zeros(N,1); numj = zeros(N,1);
cumsum_Lw = cumsum(Lw);
r = zeros(sum(ind_ths),1);
cont = 1;
for k=1:N
    medj(k) = mean(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1)))+cumsum_Lw(k)));
    sigj(k) = std(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1)))+cumsum_Lw(k)));
    numj(k) = sum(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1)));
    r(cont:sum(numj)) = (A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1)))+cumsum_Lw(k)) - medj(k))/sigj(k);
    cont = sum(numj) + 1;
end

rs = r(unidrnd(length(r),length(r),1));
Aths2 = Aths;
cumsum_numj = cumsum(numj);
cont = 1;
for k=1:N
    Aths2(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1)))+cumsum_Lw(k)) = rs(cont:cumsum_numj(k))*sigj(k) + medj(k);
    cont = cumsum_numj(k) + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 3; % tamanho amostral
nt = 64; % número de pontos para cada dia
Nboot = 3; % número de réplicas bootstrap

%rng('default')
xi1 = simulate(model_xi1,n);
xi2 = simulate(model_xi2,n);

% gerarndo a função não observada
X = zeros(nt,n);
for ii=1:n
    X(:,ii) = xi1(ii)*sqrt(2)*cos(pi*u) + xi2(ii)*sqrt(2)*cos(2*pi*u);
end

% gerando o ruído a ser adicionado na função não observada
mEps = zeros(nt,n);
for ii=1:n
    for jj=1:10
        mEps(:,ii) = mEps(:,ii) +...
            normrnd(0,1)*sqrt(2)*sin(pi*u*jj)/(2^(jj-1));
    end
end

% funcional observado
Y = X + mEps;

% número de coefientes da base de ondaletas usada
J = 71;
% matriz de cujas colunas terão a decomposição em ondaletas
% de cada funcional observado (coluna da matriz Y)
A = zeros(J,n);
Aths = zeros(J,n);
Aboot = zeros(J,n,Nboot);
ind_ths = zeros(J,n);

medj = zeros(N,n); sigj = zeros(N,n); numj = zeros(N,n);
r = zeros(J,n);

% decomposição de cada funcional dos n dias
for ii = 1:n
    % D aqui irá receber os coeficiente de aproximação e L
    % os coeficientes de detalhes
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    [thr,sorh,keepapp] = ddencmp('den','wv',Y(:,ii));
    [~,Aths(:,ii),~,~,~] = wdencmp('gbl',A(:,ii),Lw,wname,1,thr,'h',keepapp);
    [~,Aths(:,ii),~] = wden(Y(:,ii),'sqtwolog','h','mln',N,'db2');
    ind_ths(:,ii) = (A(:,ii) ~= Aths(:,ii));
    
    cumsum_Lw = cumsum(Lw);
    cont = 1;
    for k=1:N
        medj(k,ii) = mean(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        sigj(k,ii) = std(A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii));
        numj(k,ii) = sum(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii));
        r(cont:sum(numj(:,ii)), ii) = (A(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k), ii) - medj(k, ii))/sigj(k, ii);
        cont = sum(numj(:,ii)) + 1;
    end
    
    cumsum_numj = cumsum(numj(:,ii));
    for bb=1:Nboot
        rs = r(unidrnd(sum(ind_ths(:,ii)),sum(ind_ths(:,ii)),1), ii);
        Aboot(:,ii,bb) = Aths(:,ii);
        cont = 1;
        for k=1:N
            Aboot(find(ind_ths((cumsum_Lw(k)+1):cumsum_Lw(k+1),ii))+cumsum_Lw(k),ii,bb) = rs(cont:cumsum_numj(k))*sigj(k,ii) + medj(k,ii);
            cont = cumsum_numj(k) + 1;
        end
    end

end





