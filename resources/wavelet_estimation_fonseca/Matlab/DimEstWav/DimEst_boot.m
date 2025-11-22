function [ vd0p1_boot ] = DimEst_boot( Y, NREP, B, p, N, wname)
% o vetor vd0_boot conterá as dimensões bootstrap (d0+1) estimadas quando
% o processo tem dimensão d0 (número de colunas de B)
% Y é a matriz de funcionais observados, nas colunas, para diferentes
% dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. wname é o nome
% da base de ondaletas utilizada.

[J,d0] = size(B);
[nt,n] = size(Y);
vd0p1_boot = zeros(NREP,1);
Yhat = zeros(nt,n);
H = zeros(nt,d0);
A = zeros(J,n);
Yboot = zeros(nt,n);
Aboot = zeros(J,n);

for ii = 1:n
    [A(:,ii),Lw] = wavedec(Y(:,ii),N,wname);
    %[~,A(:,ii),Lw] = wden(Y(:,ii),'sqtwolog','h','mln',N,wname);
end
mu_A = mean(A,2);
C = A - mu_A*ones(1,n);

% tomando a transformada inversa da função média
mu_hat =  waverec(mu_A,Lw,wname);

for ii=1:d0
    H(:,ii) = waverec(B(:,ii),Lw,wname);
end

for ii=1:n
    Yhat(:,ii) = mu_hat;
    for k=1:d0
        Yhat(:,ii) = Yhat(:,ii) + (C(:,ii)'*B(:,k))*H(:,k);
    end
end

mEps_hat = Y - Yhat;

for jj=1:NREP

for ii=1:n
    Yboot(:,ii) = Yhat(:,ii) + mEps_hat(:,unidrnd(n));
    Aboot(:,ii) = wavedec(Yboot(:,ii),N,wname);
    %[~,Aboot(:,ii),Lw] = wden(Yboot(:,ii),'sqtwolog','h','mln',N,wname);
end

mu_Aboot = mean(Aboot,2);
C = Aboot - mu_Aboot*ones(1,n);

C1 = C(:,[1:(n-p)]);
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,[(k+1):(n-p+k)])'*C(:,[(k+1):(n-p+k)]);
end
Dboot = C1*D1*C1'/((n-p)^2);

[Bboot,Lboot] = eig(Dboot);
vd0p1_boot(jj) = Lboot(d0+1,d0+1);

end

end

