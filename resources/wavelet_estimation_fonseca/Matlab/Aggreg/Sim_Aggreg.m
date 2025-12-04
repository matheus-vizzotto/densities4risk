% Estudo de simulação da dimensão do processo de séries temporais de curvas
% usando ondaletas. Nesse código é utilizada a função func_Sim_BYZ2010 e
% cálculos em paralelo. As dimensões utilizadas variam em d=2,4,6 e os
% tamanhos amostrais usados são n=100,300,600. Ao final da simulação temos
% um arquivo para cada dimensão, com os dados gerados para cada cada valor
% de n.
% A técnica utilizada foi o método que utiliza agregação.
% Na simulaçao obtemos para cada replica os oito maiores autovalores
% gerados, valores-p obtidos, dimensoes estimadas
% os resultados sao armazenados nos arquivos file_Sim_BYZ2010_d2.mat,
% file_Sim_BYZ2010_d4.mat e file_Sim_BYZ2010_d6.mat

NREP = 200;  % número de réplicas da simulação
Nboot = 200;  % número de réplicas bootstrap
nt = 500;  % número de pontos para cada dia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulações para dimensão igual a 2
d = 2;

%%%%%%
n=100;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n100=vDimSel; mPvalues_n100=mPvalues; mEigval_n100=mEigval;
mEQMI_n_n100=mEQMI_n; mEQMI_nt_n100=mEQMI_nt;

%%%%%%
n=300;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n300=vDimSel; mPvalues_n300=mPvalues; mEigval_n300 = mEigval;
mEQMI_n_n300=mEQMI_n; mEQMI_nt_n300=mEQMI_nt;

save('file_Sim_Aggreg_d2_part1.mat','vDimSel_n100','mPvalues_n100','mEigval_n100','mEQMI_n_n100','mEQMI_nt_n100',...
    'vDimSel_n300','mPvalues_n300','mEigval_n300','mEQMI_n_n300','mEQMI_nt_n300');  

%%%%%%
n=600;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n600=vDimSel; mPvalues_n600=mPvalues; mEigval_n600 = mEigval;
mEQMI_n_n600=mEQMI_n; mEQMI_nt_n600=mEQMI_nt;
%%%%%%

% salvando as matrizes
save('file_Sim_Aggreg_d2_part2.mat','vDimSel_n600','mPvalues_n600','mEigval_n600','mEQMI_n_n600','mEQMI_nt_n600');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulações para dimensão igual a 4
d = 4;

%%%%%%
n=100;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n100=vDimSel; mPvalues_n100=mPvalues; mEigval_n100=mEigval;
mEQMI_n_n100=mEQMI_n; mEQMI_nt_n100=mEQMI_nt;

%%%%%%
n=300;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n300=vDimSel; mPvalues_n300=mPvalues; mEigval_n300 = mEigval;
mEQMI_n_n300=mEQMI_n; mEQMI_nt_n300=mEQMI_nt;

save('file_Sim_Aggreg_d4_part1.mat','vDimSel_n100','mPvalues_n100','mEigval_n100','mEQMI_n_n100','mEQMI_nt_n100',...
    'vDimSel_n300','mPvalues_n300','mEigval_n300','mEQMI_n_n300','mEQMI_nt_n300');  

%%%%%%
n=600;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n600=vDimSel; mPvalues_n600=mPvalues; mEigval_n600 = mEigval;
mEQMI_n_n600=mEQMI_n; mEQMI_nt_n600=mEQMI_nt;
%%%%%%

% salvando as matrizes
save('file_Sim_Aggreg_d4_part2.mat','vDimSel_n600','mPvalues_n600','mEigval_n600','mEQMI_n_n600','mEQMI_nt_n600');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulações para dimensão igual a 6
d = 6;

%%%%%%
n=100;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n100=vDimSel; mPvalues_n100=mPvalues; mEigval_n100=mEigval;
mEQMI_n_n100=mEQMI_n; mEQMI_nt_n100=mEQMI_nt;

%%%%%%
n=300;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n300=vDimSel; mPvalues_n300=mPvalues; mEigval_n300 = mEigval;
mEQMI_n_n300=mEQMI_n; mEQMI_nt_n300=mEQMI_nt;

save('file_Sim_Aggreg_d6_part1.mat','vDimSel_n100','mPvalues_n100','mEigval_n100','mEQMI_n_n100','mEQMI_nt_n100',...
    'vDimSel_n300','mPvalues_n300','mEigval_n300','mEQMI_n_n300','mEQMI_nt_n300');  

%%%%%%
n=600;
mPvalues = ones(NREP,8); mEigval = zeros(NREP,8); vDimSel = 9*ones(NREP,1);
mEQMI_n = zeros(NREP,n); mEQMI_nt = zeros(NREP,nt);

tic
parfor k=1:NREP
    [vDimSel(k),mEQMI_n(k,:),mEQMI_nt(k,:),mEigval(k,:),mPvalues(k,:)] = func_Sim_Aggreg(n,d,1,Nboot,nt);
end
toc
vDimSel_n600=vDimSel; mPvalues_n600=mPvalues; mEigval_n600 = mEigval;
mEQMI_n_n600=mEQMI_n; mEQMI_nt_n600=mEQMI_nt;
%%%%%%

% salvando as matrizes
save('file_Sim_Aggreg_d6_part2.mat','vDimSel_n600','mPvalues_n600','mEigval_n600','mEQMI_n_n600','mEQMI_nt_n600');  


