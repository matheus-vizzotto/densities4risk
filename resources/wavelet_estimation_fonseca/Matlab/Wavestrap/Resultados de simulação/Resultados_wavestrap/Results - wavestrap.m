% Arquivo para ler e analisar os dados simulados

fd2 = matfile('file_Sim_wavestrap_d2.mat');
fd4 = matfile('file_Sim_wavestrap_d4.mat');
fd6 = matfile('file_Sim_wavestrap_d6.mat');


mTbn100d2 = tabulate(fd2.vDimSel_n100);mTbn300d2 = tabulate(fd2.vDimSel_n300);mTbn600d2 = tabulate(fd2.vDimSel_n600);
mTbn100d4 = tabulate(fd4.vDimSel_n100);mTbn300d4 = tabulate(fd4.vDimSel_n300);mTbn600d4 = tabulate(fd4.vDimSel_n600);
mTbn100d6 = tabulate(fd6.vDimSel_n100);mTbn300d6 = tabulate(fd6.vDimSel_n300);mTbn600d6 = tabulate(fd6.vDimSel_n600);

mTbn100d2
mTbn300d2
mTbn600d2

mTbn100d4
mTbn300d4
mTbn600d4

mTbn100d6
mTbn300d6
mTbn600d6


% boxplot dos valores-p dos testes do d-ésimo e (d+1)-ésimo autovalor, para
% cada tamanho amostral e valor verdadeiro de d

mPvn100d2 = fd2.mPvalues_n100 ;mPvn300d2 = fd2.mPvalues_n300 ;mPvn600d2 = fd2.mPvalues_n600;
mPvn100d4 = fd4.mPvalues_n100 ;mPvn300d4 = fd4.mPvalues_n300 ;mPvn600d4 = fd4.mPvalues_n600;
mPvn100d6 = fd6.mPvalues_n100 ;mPvn300d6 = fd6.mPvalues_n300 ;mPvn600d6 = fd6.mPvalues_n600;

subplot(2,1,1)
aux = [mPvn100d2(:,1),mPvn300d2(:,1),mPvn600d2(:,1),...
    mPvn100d4(:,3),mPvn300d4(:,3),mPvn600d4(:,3),...
    mPvn100d6(:,5),mPvn300d6(:,5),mPvn600d6(:,5)];
xt = {'\begin{tabular}{c}d=2\\n=100\end{tabular}';'\begin{tabular}{c}d=2\\n=300\end{tabular}';'\begin{tabular}{c}d=2\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=4\\n=100\end{tabular}';'\begin{tabular}{c}d=4\\n=300\end{tabular}';'\begin{tabular}{c}d=4\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=6\\n=100\end{tabular}';'\begin{tabular}{c}d=6\\n=300\end{tabular}';'\begin{tabular}{c}d=6\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(a) p-values for the d-th largest eigenvalue');
subplot(2,1,2)
aux = [mPvn100d2(:,2),mPvn300d2(:,2),mPvn600d2(:,2),...
    mPvn100d4(:,4),mPvn300d4(:,4),mPvn600d4(:,4),...
    mPvn100d6(:,6),mPvn300d6(:,6),mPvn600d6(:,6)];
xt = {'\begin{tabular}{c}d=2\\n=100\end{tabular}';'\begin{tabular}{c}d=2\\n=300\end{tabular}';'\begin{tabular}{c}d=2\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=4\\n=100\end{tabular}';'\begin{tabular}{c}d=4\\n=300\end{tabular}';'\begin{tabular}{c}d=4\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=6\\n=100\end{tabular}';'\begin{tabular}{c}d=6\\n=300\end{tabular}';'\begin{tabular}{c}d=6\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(b) p-values for the (d+1)-th largest eigenvalue');

