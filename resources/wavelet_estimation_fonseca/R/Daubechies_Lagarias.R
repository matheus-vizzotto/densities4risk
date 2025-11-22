# exemplo 3.5.2 da pág. 90 do 'Statistical modeling by wavelets'

# filtro de ondaletas da DAUB2
h <- c(.482963,.836516,.224144,-.12941)
# aqui temos N=2
N <- length(h)/2
# matrizes usadas no produto com os diádicos, formadas através dos elementos
# do filtro de ondaletas
T0 <- matrix(0,2*N-1,2*N-1)
T1 <- matrix(0,2*N-1,2*N-1)
for(i in 1:(2*N-1)){
  for(j in 1:(2*N-1)){
    # caso seja tomado um índice do filtro fora de {0,...,2N-1}, é 
    # atribuído um valor zero ao elemento da matriz
    h_ind <- ifelse((2*i-j-1<0),0,ifelse((2*i-j-1>2*N-1),0,h[2*i-j]))
    T0[i,j] <- sqrt(2)*h_ind
    h_ind <- ifelse((2*i-j<0),0,ifelse((2*i-j>2*N-1),0,h[2*i-j+1]))
    T1[i,j] <- sqrt(2)*h_ind
  }
}
# matrizes obtidas
T0
T1

# essa função retornará os n primeiros valores da representação diádica de um
# número x em (0,1)
dyadn <- function(x,n){
  d <- rep(0,n)
  xd <- x
  for(i in 1:n){
    d[i] <- floor(2*xd)
    xd = 2*xd - floor(2*xd)
  }
  return(d)
}
# essa função retornará a representação diádica da parte inteira de
# um número x>0
dyad_int <- function(x){
  d = numeric(0)
  xi <- floor(x)
  i = 1
  d[i] <- xi%%2
  while(xi>1){
    d[i] <- xi%%2
    xi <- xi%/%2
    i = i + 1
  }
  d[i] <- xi
  return(d[length(d):1])
}

# essa função retornará os n primeiros valores da representação diádica de
# um número x>0
dyadn_pos <- function(x,n){
  if(x<1){
    d <- dyadn(x,n)
  }else{
    xi <- floor(x)
    xd <- x - floor(x)
    di <- dyad_int(xi)
    if(length(di)>=n){
      d <- di[1:n]
    }else{
      k <- n - length(di)
      d <- rep(0,n)
      d[1:length(di)] <- di
      d[(length(di)+1):n] <- dyadn(xd,k)
    }
  }
  return(d)
}

vd <- dyadn_pos(.45,20)

# agora realizamos o produto das matrizes T0 e T1 conforme a representação
# diádica de x
if(vd[1]==0) mT <- T0 else mT <- T1  
for(i in 2:length(vd)){
  if(vd[i]==0) mT <- mT%*%T0 else mT <- mT%*%T1
}
mT # resultado

# função para calcular a função escala avaliada em x com aproximação
# de ordem n nos diáticos e para um filtro de ondaletas h
phi.n <- function(x,n,h){
  if(x<=0 || x>=length(h)-1){
    return(0)
  }else{
    k <- 1
    y <- x
    while(y>=0){
      if(y<1){
        vd <- dyadn(y,n)
        if(vd[1]==0) mT <- T0 else mT <- T1  
        for(i in 2:length(vd)) if(vd[i]==0) mT <- mT%*%T0 else mT <- mT%*%T1
        return(mT[k,1])
      }else{
        y <- y - 1
        k <- k + 1
      }
    }
  }
}

phi.n(2.7389,20,h)

# função para calcular a função de ondaleta psi avaliada em x com aproximação
# de ordem n nos diáticos e para um filtro de ondaletas h
psi.n <- function(x,n,h){
  if(x<=1-length(h)/2 || x>=length(h)/2){
    return(0)
  }else{
    psiv = 0
    for(i in (-1):(length(h)-2)) psiv = psiv + 
        ((-1)^i)*sqrt(2)*h[2+i]*phi.n(2*x+i,n,h)
    return(psiv)
  }
}

phi.n(.45,20,h)
psi.n(.45,20,h)

############

# exemplo para o teorema 3.5.5, para avaliar a função de ondaleta (mãe)
# em um valor x

x <- 0.45
vu <- rep(0,2*N-1)
for(i in 0:(2*N-2)){
  # se para algum i o índice i+1-floor(2x) é negativo ou maior que 2N-1,
  # então o componente de vu é zero
  h_ind <- ifelse((i+1-floor(2*x)<0),0,ifelse((i+1-floor(2*x)>2*N-1),0,h[i+1-floor(2*x)+1]))
  vu[i+1] <- ((-1)^(1-floor(2*x)))*h_ind
}
# n primeiros valores da representação diádica de 2x
vd2x <- dyadn_pos(2*x,20)
# produto das matrizes T0 e T1 para 2x
if(vd2x[1]==0) mT2x <- T0 else mT2x <- T1  
for(i in 2:length(vd2x)){
  if(vd2x[i]==0) mT2x <- mT2x%*%T0 else mT2x <- mT2x%*%T1
}
# vetor v
vxn <- mT2x%*%rep(1,2*N-1)/(2*N-1)
# aproximação da função de ondaleta avaliada em x
vu%*%vxn


########
# outra forma

x <- 0.45
vu <- rep(0,2*N-1)
for(i in (-1):(2*N-2)) vu[i+2] <- ((-1)^i)*h[i+2]

# n primeiros valores da representação diádica de 2x
vd2x <- dyadn_pos(2*x,20)
if(vd2x[1]==0) mT2x <- T0 else mT2x <- T1  
for(i in 2:length(vd2x)){
  if(vd2x[i]==0) mT2x <- mT2x%*%T0 else mT2x <- mT2x%*%T1
}

vxn <- c(0,mT2x%*%rep(1,2*N-1)/(2*N-1))
sqrt(2)*vu%*%vxn

#########
x <- .45
psi_value = 0
for(i in (-1):(2*N-2)){
  if(2*x+i<0||2*x+i>2*N-1){
    psi_value = psi_value + 0
  }else{
    vd <- dyadn_pos(2*x+i,20)
    if(vd[1]==0) mT <- T0 else mT <- T1  
    for(j in 2:length(vd)){
      if(vd[j]==0) mT <- mT%*%T0 else mT <- mT%*%T1
    }
    psi_value = psi_value + ((-1)^(-i))*h[2+i]*sqrt(2)*mT[1,1]
  }
}
psi_value
  

