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
line([0:10],.05*ones(11,1),'LineStyle','--');
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
line([0:10],.05*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(b) p-values for the (d+1)-th largest eigenvalue');


% Gráfico do EQM das funções estimadas ao longo das n observações

subplot(3,3,1)
plot(1:100,mean(fd2.mEQMI_n_n100));xlim([1 100]);ylim([.2 2]);
title('d=2, n=100');
%hline = refline([0 mean(mean(fd2.mEQMI_n_n100))]);hline.Color = 'r';
subplot(3,3,2)
plot(1:300,mean(fd2.mEQMI_n_n300));xlim([1 300]);ylim([.2 1]);
title('d=2, n=300');
subplot(3,3,3)
plot(1:600,mean(fd2.mEQMI_n_n600));xlim([1 600]);ylim([.2 1]);
title('d=2, n=600');

subplot(3,3,4)
plot(1:100,mean(fd4.mEQMI_n_n100));xlim([1 100]);ylim([1 3]);
title('d=4, n=100');
subplot(3,3,5)
plot(1:300,mean(fd4.mEQMI_n_n300));xlim([1 300]);ylim([.2 1]);
title('d=4, n=300');
subplot(3,3,6)
plot(1:600,mean(fd4.mEQMI_n_n600));xlim([1 600]);ylim([.2 1]);
title('d=6, n=600');

subplot(3,3,7)
plot(1:100,mean(fd6.mEQMI_n_n100));xlim([1 100]);ylim([1 5]);
title('d=6, n=100');
subplot(3,3,8)
plot(1:300,mean(fd6.mEQMI_n_n300));xlim([1 300]);ylim([.2 1]);
title('d=6, n=300');
subplot(3,3,9)
plot(1:600,mean(fd6.mEQMI_n_n600));xlim([1 600]);ylim([.2 1]);
title('d=6, n=600');


% Gráfico do EQM das funções estimadas ao longo das nt observações

subplot(3,3,1)
plot(1:256,mean(fd2.mEQMI_nt_n100));xlim([1 256]);ylim([0 4]);
title('d=2, n=100');
subplot(3,3,2)
plot(1:256,mean(fd2.mEQMI_nt_n300));xlim([1 256]);ylim([0 1]);
title('d=2, n=300');
subplot(3,3,3)
plot(1:256,mean(fd2.mEQMI_nt_n600));xlim([1 256]);ylim([0 1]);
title('d=2, n=600');

subplot(3,3,4)
plot(1:256,mean(fd4.mEQMI_nt_n100));xlim([1 256]);ylim([0 4]);
title('d=4, n=100');
subplot(3,3,5)
plot(1:256,mean(fd4.mEQMI_nt_n300));xlim([1 256]);ylim([0 1]);
title('d=4, n=300');
subplot(3,3,6)
plot(1:256,mean(fd4.mEQMI_nt_n600));xlim([1 256]);ylim([0 1]);
title('d=4, n=600');

subplot(3,3,7)
plot(1:256,mean(fd6.mEQMI_nt_n100));xlim([1 256]);ylim([0 6]);
title('d=6, n=100');
subplot(3,3,8)
plot(1:256,mean(fd6.mEQMI_nt_n300));xlim([1 256]);ylim([0 2]);
title('d=6, n=300');
subplot(3,3,9)
plot(1:256,mean(fd6.mEQMI_nt_n600));xlim([1 256]);ylim([0 2]);
title('d=6, n=600');




