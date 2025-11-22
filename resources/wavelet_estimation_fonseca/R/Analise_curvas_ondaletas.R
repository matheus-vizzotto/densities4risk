# Código para simular a identificação da dimensionalidade de curvas via 
# ondaletas

# Essa será a densidade verdadeira, de dimensão 2, de uma variável aleatória
# que assume valores em (0,1)
f_t <- function(x){
  (2/3)*(pi/2)*cos(x*pi/2) + (1/3)*(pi*sqrt(2)/4)*cos(x*pi/4)
}

# o erro presente na densidade observada será gerado por sin((x-.5)/pi)

# A ideia para obter a densidade observada g_t em um dia é gerar várias 
# observações unif(0,1), avaliar em g_t e formar um histograma

set.seed(2017)
x <- runif(1000)
y <- f_t(x) + sin((x-.5)/pi)
plot(x,y)


############################################

require(wavethresh)
data <- rclaw(100)
datahr <- denproj(data, J=8, filter.number=2,family="DaubExPhase")
data.wd <- denwd(datahr)
hist(data)
plotdenwd(data.wd, top.level=(datahr$res$J-1))

