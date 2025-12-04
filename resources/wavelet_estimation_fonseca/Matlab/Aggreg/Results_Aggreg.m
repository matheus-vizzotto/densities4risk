
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for case 1
% the weights used in this case are (0.1,0.1,0.8)

fd2_c1p1 = matfile('simulation_case1/file_Sim_Aggreg_d2_part1.mat');
fd2_c1p2 = matfile('simulation_case1/file_Sim_Aggreg_d2_part2.mat');
fd4_c1p1 = matfile('simulation_case1/file_Sim_Aggreg_d4_part1.mat');
fd4_c1p2 = matfile('simulation_case1/file_Sim_Aggreg_d4_part2.mat');
fd6_c1p1 = matfile('simulation_case1/file_Sim_Aggreg_d6_part1.mat');
fd6_c1p2 = matfile('simulation_case1/file_Sim_Aggreg_d6_part2.mat');


mTbn100d2_c1 = tabulate(fd2_c1p1.vDimSel_n100);
mTbn300d2_c1 = tabulate(fd2_c1p1.vDimSel_n300);
mTbn600d2_c1 = tabulate(fd2_c1p2.vDimSel_n600);
mTbn100d4_c1 = tabulate(fd4_c1p1.vDimSel_n100);
mTbn300d4_c1 = tabulate(fd4_c1p1.vDimSel_n300);
mTbn600d4_c1 = tabulate(fd4_c1p2.vDimSel_n600);
mTbn100d6_c1 = tabulate(fd6_c1p1.vDimSel_n100);
mTbn300d6_c1 = tabulate(fd6_c1p1.vDimSel_n300);
mTbn600d6_c1 = tabulate(fd6_c1p2.vDimSel_n600);

mTbn100d2_c1
mTbn300d2_c1
mTbn600d2_c1

mTbn100d4_c1
mTbn300d4_c1
mTbn600d4_c1

mTbn100d6_c1
mTbn300d6_c1
mTbn600d6_c1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for case 2
% the weights used in this case are (1/3,1/3,1/3)

fd2_c2p1 = matfile('simulation_case2/file_Sim_Aggreg_d2_part1.mat');
fd2_c2p2 = matfile('simulation_case2/file_Sim_Aggreg_d2_part2.mat');
fd4_c2p1 = matfile('simulation_case2/file_Sim_Aggreg_d4_part1.mat');
fd4_c2p2 = matfile('simulation_case2/file_Sim_Aggreg_d4_part2.mat');
fd6_c2p1 = matfile('simulation_case2/file_Sim_Aggreg_d6_part1.mat');
fd6_c2p2 = matfile('simulation_case2/file_Sim_Aggreg_d6_part2.mat');


mTbn100d2_c2 = tabulate(fd2_c2p1.vDimSel_n100);
mTbn300d2_c2 = tabulate(fd2_c2p1.vDimSel_n300);
mTbn600d2_c2 = tabulate(fd2_c2p2.vDimSel_n600);
mTbn100d4_c2 = tabulate(fd4_c2p1.vDimSel_n100);
mTbn300d4_c2 = tabulate(fd4_c2p1.vDimSel_n300);
mTbn600d4_c2 = tabulate(fd4_c2p2.vDimSel_n600);
mTbn100d6_c2 = tabulate(fd6_c2p1.vDimSel_n100);
mTbn300d6_c2 = tabulate(fd6_c2p1.vDimSel_n300);
mTbn600d6_c2 = tabulate(fd6_c2p2.vDimSel_n600);

mTbn100d2_c2
mTbn300d2_c2
mTbn600d2_c2

mTbn100d4_c2
mTbn300d4_c2
mTbn600d4_c2

mTbn100d6_c2
mTbn300d6_c2
mTbn600d6_c2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for case 3
% the weights used in this case are (0.1,0.3,0.5)

fd2_c3p1 = matfile('simulation_case3/file_Sim_Aggreg_d2_part1.mat');
fd2_c3p2 = matfile('simulation_case3/file_Sim_Aggreg_d2_part2.mat');
fd4_c3p1 = matfile('simulation_case3/file_Sim_Aggreg_d4_part1.mat');
fd4_c3p2 = matfile('simulation_case3/file_Sim_Aggreg_d4_part2.mat');
fd6_c3p1 = matfile('simulation_case3/file_Sim_Aggreg_d6_part1.mat');
fd6_c3p2 = matfile('simulation_case3/file_Sim_Aggreg_d6_part2.mat');


mTbn100d2_c3 = tabulate(fd2_c3p1.vDimSel_n100);
mTbn300d2_c3 = tabulate(fd2_c3p1.vDimSel_n300);
mTbn600d2_c3 = tabulate(fd2_c3p2.vDimSel_n600);
mTbn100d4_c3 = tabulate(fd4_c3p1.vDimSel_n100);
mTbn300d4_c3 = tabulate(fd4_c3p1.vDimSel_n300);
mTbn600d4_c3 = tabulate(fd4_c3p2.vDimSel_n600);
mTbn100d6_c3 = tabulate(fd6_c3p1.vDimSel_n100);
mTbn300d6_c3 = tabulate(fd6_c3p1.vDimSel_n300);
mTbn600d6_c3 = tabulate(fd6_c3p2.vDimSel_n600);

mTbn100d2_c3
mTbn300d2_c3
mTbn600d2_c3

mTbn100d4_c3
mTbn300d4_c3
mTbn600d4_c3

mTbn100d6_c3
mTbn300d6_c3
mTbn600d6_c3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for case 4
% the weight used in this case is (1)

fd2_c4p1 = matfile('simulation_case4/file_Sim_Aggreg_d2_part1.mat');
fd2_c4p2 = matfile('simulation_case4/file_Sim_Aggreg_d2_part2.mat');
fd4_c4p1 = matfile('simulation_case4/file_Sim_Aggreg_d4_part1.mat');
fd4_c4p2 = matfile('simulation_case4/file_Sim_Aggreg_d4_part2.mat');
fd6_c4p1 = matfile('simulation_case4/file_Sim_Aggreg_d6_part1.mat');
fd6_c4p2 = matfile('simulation_case4/file_Sim_Aggreg_d6_part2.mat');


mTbn100d2_c4 = tabulate(fd2_c4p1.vDimSel_n100);
mTbn300d2_c4 = tabulate(fd2_c4p1.vDimSel_n300);
mTbn600d2_c4 = tabulate(fd2_c4p2.vDimSel_n600);
mTbn100d4_c4 = tabulate(fd4_c4p1.vDimSel_n100);
mTbn300d4_c4 = tabulate(fd4_c4p1.vDimSel_n300);
mTbn600d4_c4 = tabulate(fd4_c4p2.vDimSel_n600);
mTbn100d6_c4 = tabulate(fd6_c4p1.vDimSel_n100);
mTbn300d6_c4 = tabulate(fd6_c4p1.vDimSel_n300);
mTbn600d6_c4 = tabulate(fd6_c4p2.vDimSel_n600);

mTbn100d2_c4
mTbn300d2_c4
mTbn600d2_c4

mTbn100d4_c4
mTbn300d4_c4
mTbn600d4_c4

mTbn100d6_c4
mTbn300d6_c4
mTbn600d6_c4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boxplot dos valores-p dos testes do d-ésimo e (d+1)-ésimo autovalor, para
% cada tamanho amostral e valor verdadeiro de d

mPvn100d4_c1 = fd4_c1p1.mPvalues_n100;
mPvn300d4_c1 = fd4_c1p1.mPvalues_n300;
mPvn600d4_c1 = fd4_c1p2.mPvalues_n600;

mPvn100d4_c4 = fd4_c4p1.mPvalues_n100;
mPvn300d4_c4 = fd4_c4p1.mPvalues_n300;
mPvn600d4_c4 = fd4_c4p2.mPvalues_n600;


subplot(2,1,1)
aux = [mPvn100d4_c1(1:200,3),mPvn100d4_c4(:,3),mPvn300d4_c1(1:200,3),...
    mPvn300d4_c4(:,3),mPvn600d4_c1(1:200,3),mPvn600d4_c4(:,3)];
xt = {'\begin{tabular}{c}case 1\\n=100\end{tabular}';'\begin{tabular}{c}case 4\\n=100\end{tabular}';...
    '\begin{tabular}{c}case 1\\n=300\end{tabular}';'\begin{tabular}{c}case 4\\n=300\end{tabular}';...
    '\begin{tabular}{c}case 1\\n=600\end{tabular}';'\begin{tabular}{c}case 4\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(a) p-values for the 4-th largest eigenvalue');
subplot(2,1,2)
aux = [mPvn100d4_c1(1:200,4),mPvn100d4_c4(:,4),mPvn300d4_c1(1:200,4),...
    mPvn300d4_c4(:,4),mPvn600d4_c1(1:200,4),mPvn600d4_c4(:,4)];
xt = {'\begin{tabular}{c}case 1\\n=100\end{tabular}';'\begin{tabular}{c}case 4\\n=100\end{tabular}';...
    '\begin{tabular}{c}case 1\\n=300\end{tabular}';'\begin{tabular}{c}case 4\\n=300\end{tabular}';...
    '\begin{tabular}{c}case 1\\n=600\end{tabular}';'\begin{tabular}{c}case 4\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(b) p-values for the 5-th largest eigenvalue');

