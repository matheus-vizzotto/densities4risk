gs_coef = zeros(32,1);
gs_d = [1,1,2,4,8,16,32];

gs_coef(1) = normrnd(0,1,1,1);
for ii = 0:4
	gs_coef((2^ii + 1):(2^(ii+1))) = normrnd(0,1,2^ii,1)/(2^ii);
end
gs_coefnorm = gs_coef/(norm(gs_coef));

N = 2;
j0 = 0;
j1 = 5;
x = .49;
n = 15;
wname = 'db2';

n_ind_j0_phi = 2^j0 + 2*N;
cs_phi = normrnd(0,1,n_ind_j0_phi,1)/(2^j0);

cs_psi = zeros(2^(j1+1)-2^j0+2*N*(j1-j0+1),1);
cs_psi_dim = zeros(j1-j0+1,1);

cs_gs = [cs_phi',cs_psi']';
cs_gs = cs_gs/norm(cs_gs);

k = 1;
for jj = j0:j1
    cs_psi_dim(jj-j0+1) = 2^jj + 2*N;
    cs_psi(k:sum(cs_psi_dim)) = normrnd(0,1,cs_psi_dim(jj-j0+1),1)/(2^jj);
    k = sum(cs_psi_dim) + 1;
end

xv = linspace(0.01,.99,500);
gsv = zeros(500,1);

for ll=1:500
x = xv(ll);    
gs = 0;
for k = (1-2*N):(2^j0)
    gs = gs + cs_phi(k+2*N)*(2^(j0/2))*phi_n((2^j0)*x-k,n,wname);
end

ki = 1;
for jj = j0:j1
    for k = (-N):(2^jj + N-1)
        gs = gs + cs_psi(ki)*(2^(jj/2))*psi_n((2^jj)*x-k,n,wname);
        ki = ki + 1;
    end
end
gsv(ll) = gs;
end


[cs_phi,cs_psi] = wavcoef_rand_sqrtdens(j0,j1,wname);
rand_dens(x,j0,j1,n,wname,cs_phi,cs_psi)

xv = linspace(0,1,500);
gv = zeros(500,1);
for ii = 1:500
    gv(ii) = rand_dens(xv(ii),j0,j1,n,wname,cs_phi,cs_psi);
end
plot(xv,gv)

