function [ cs_phi,cs_psi ] = wavcoef_rand_sqrtdens(j0,j1,wname)
% para os níveis j0 e j1 da decomposição de ondaletas, sendo wname
% alguma ondaleta de Daubechies, essa função gera aleatoriamente 
% coeficientes da decomposição de ondaletas de alguma densidade.

h = wfilters(wname);
N = length(h)/2;

n_ind_j0_phi = 2^j0 + 2*N;
cs_phi = normrnd(0,1,n_ind_j0_phi,1)/(2^j0);

cs_psi = zeros(2^(j1+1)-2^j0+2*N*(j1-j0+1),1);
cs_psi_dim = zeros(j1-j0+1,1);

k = 1;
for jj = j0:j1
    cs_psi_dim(jj-j0+1) = 2^jj + 2*N;
    cs_psi(k:sum(cs_psi_dim)) = normrnd(0,1,cs_psi_dim(jj-j0+1),1)/(2^jj);
    k = sum(cs_psi_dim) + 1;
end

cs_gs = [cs_phi',cs_psi']';
cs_phi = cs_phi/norm(cs_gs);
cs_psi = cs_psi/norm(cs_gs);

end

