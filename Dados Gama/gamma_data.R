#Author: Eduardo Ensslin
#Script para a reprodução dos resultados obtidos no documento
#de defesa de qualificação de doutorado.

#carregando pacotes
library(dplyr)
library(scales)
require(ggplot2)
require(ggthemes)
theme_set(theme_pander())
require(reshape2)

#semente para a reprodução dos dados
set.seed(123456, kind="Mersenne-Twister")

###################################################
#     DADOS SIMULADOS COM DISTRIBUICAO GAMMA      #
###################################################

#Definindo o valor dos parâmetros de escala e tamanho amostral
l1 <- 1
l2 <- 10
n1 <- 50
n2 <- 50

#Gerando os dados simulados ~ gamma
z <- c(rgamma(n1, l1,2), rgamma(n2, l2,2))
z
ggplot(data.frame(z), aes(x=1:100, y=z)) +
  geom_line() +
  geom_vline(xintercept = 50, col="darkorange") +
  xlab("N") + ylab("f(z)")

#Estimação dos valores para uma distribuição gamma
Estim <- matrix(data=rep(0, 100), nrow=3, ncol=100)
Estim[1,] <- 1:100

for(el in 5:95) {
  Estim[2, el] <- (mean(z[1:el]))/2
  Estim[3, el] <- (mean(z[(el+1):100]))/2
}

Estim.df <- melt(data.frame(Position=Estim[1,], lambda1=Estim[2,], lambda2=Estim[3,]), 
                 measure.vars = 2:3, 
                 value.name = "Estimates")
#Estim.df
ggplot(Estim.df, aes(x=Position, y=Estimates, group=variable, col=variable)) +
  geom_line()


#Conta auxiliar para o cômputo da estatística de teste
N <- ((n1*n2)/(n1+n2))

#Definindo os valores fixados
k=2
beta=0.5

###################################################
#        SHANNON - DADOS SIMULADOS GAMMA          #
###################################################

#Cálculo estatística de teste
SH <- rep(0, 100)

for(el in 5:95) {
  lambda1 <- Estim[2,el]
  lambda2 <- Estim[3,el]
  SH[el] = N*(abs(sqrt(k))*(log(lambda2)-log(lambda1)))^2
}
edge <- which.max(SH)
#windows()
ggplot(data.frame(SH), aes(x=1:100, y=SH)) +
  geom_line() +
  xlab("n") + ylab(expression(T["(h,phi)"]))+
  geom_hline(yintercept = qchisq(c(.9, .95, .99), df=1),
             col="orange", linetype=2) +
  geom_vline(xintercept = edge, col="red")


###################################################
#         RÉNYI - DADOS SIMULADOS GAMMA           #                           
###################################################

#Cálculo estatística de teste
SH <- rep(0, 100)

#Contas auxiliares simplificadas
frac1 <- (1/(1-beta))
frac2 <- (1/beta)
gama1 <- gamma(beta*(k-1)+1)
gama2 <- gamma(beta*(k-1)+2)
gama3 <- gamma(beta*(k-1)+3)
renyi1 <- frac1*((gama2/gama1)^2-2*beta*k*(gama2/gama1)+beta^2*k^2)
renyi2 <- frac2*((gama3/gama1)-2*k*(gama2/gama1)+beta*k^2)  


for(el in 5:95) {
  lambda1 <- Estim[2,el]
  lambda2 <- Estim[3,el]
  SH[el] = N*(abs(sqrt(renyi1-renyi2)*(log(lambda2)-log(lambda1))))^2
}
edge <- which.max(SH)
#windows()
ggplot(data.frame(SH), aes(x=1:100, y=SH)) +
  geom_line() +
  xlab("n") + ylab(expression(T["(h,phi)"]))+
  geom_hline(yintercept = qchisq(c(.9, .95, .99), df=1),
             col="orange", linetype=2) +
  geom_vline(xintercept = edge, col="red")


###################################################
#        ARIMOTO - DADOS SIMULADOS GAMMA          #                       
###################################################

#Cálculo estatística de teste
SH <- rep(0, 100)

#Contas auxiliares simplificadas
potArmt1 <- (((((beta^2)+(beta*k)+1))/beta))
potArmt2 <- ((beta*(beta+k))/beta)
potArmt3 <- ((beta*(beta+k)-beta)/beta)
potArmt4 <- ((beta*(beta+k)-(2*beta))/beta)
gamaArmt1 <- gamma((beta+k-1)/beta)
gamaArmt2 <- gamma(((2*beta)+k-1)/beta)
gamaArmt3 <- gamma(((3*beta)+k-1)/beta)
fracArmt <- (1/gamma(k))
constArmt1 <- ((-beta^potArmt1*gamaArmt1^(beta-2)*gamaArmt2^2)
               +(2*k*beta^potArmt1*gamaArmt1^(beta-1)*gamaArmt2)
               -(k^2*beta^potArmt1*gamaArmt1^beta))
constArmt2 <- ((beta^potArmt2*gamaArmt1^(beta-1)*gamaArmt3)
               -(2*k*beta^potArmt3*gamaArmt1^(beta-1)*gamaArmt2)
               +(k^2*beta^potArmt4*gamaArmt1^beta))

for(el in 5:95) {
  lambda1 <- Estim[2,el]
  lambda2 <- Estim[3,el]
  SH[el] = N*(abs(sqrt(fracArmt*(constArmt1+constArmt2))
                  *(2/(beta-1))*(lambda2^((beta-1)/2)-lambda1^((beta-1)/2))))^2
}
edge <- which.max(SH)
#windows()
ggplot(data.frame(SH), aes(x=1:100, y=SH)) +
  geom_line() +
  xlab("n") + ylab(expression(T["(h,phi)"]))+
  geom_hline(yintercept = qchisq(c(.9, .95, .99), df=1),
             col="orange", linetype=2) +
  geom_vline(xintercept = edge, col="red")


###################################################
#      HAVRDA-CHARVAT - DADOS SIMULADOS GAMMA     #                 
###################################################

#Cálculo estatística de teste
SH <- rep(0, 100)

#Contas auxiliares simplificadas
supB <- beta*(1-k)
gamaB <- beta*(k-1)

for(el in 5:95) {
  lambda1 <- Estim[2,el]
  lambda2 <- Estim[3,el]
  SH[el] = N*(abs(sqrt((1/(gamma(k)^beta))*((beta^(supB-2)*gamma(gamaB+3))
                                            -(2*(beta^(supB-1))*k*gamma(gamaB+2))
                                            +(beta^(supB)*(k^2)*gamma(gamaB+1)))))
              *((2/(beta-1))*((lambda1^((1-beta)/2))-(lambda2^((1-beta)/2)))))^2
}
edge <- which.max(SH)
#windows()
ggplot(data.frame(SH), aes(x=1:100, y=SH)) +
  geom_line() +
  xlab("n") + ylab(expression(T["(h,phi)"]))+
  geom_hline(yintercept = qchisq(c(.9, .95, .99), df=1),
             col="orange", linetype=2) +
  geom_vline(xintercept = edge, col="red")


###################################################
#    SHARMA-MITTAL - DADOS SIMULADOS GAMMA        #                 
###################################################

#Cálculo estatística de teste
SH <- rep(0, 100)

#Contas auxiliares simplificadas
SMexp1 <- (k^2*(1-beta))*(exp((beta-1)*(k+log(gamma(k)))))
SMexp2 <- (k*(exp((beta-1)*(k+log(gamma(k)))))) 
kbeta <- (2/(k*(beta-1)))
SMpot <- ((k*(beta-1))/2)

for(el in 5:95) {
  lambda1 <- Estim[2,el]
  lambda2 <- Estim[3,el]
  SH[el] = N*(abs(sqrt(SMexp1+SMexp2)*(kbeta*((lambda2^SMpot)-(lambda1^SMpot)))))^2
}
edge <- which.max(SH)
#windows()
ggplot(data.frame(SH), aes(x=1:100, y=SH)) +
  geom_line() +
  xlab("n") + ylab(expression(T["(h,phi)"]))+
  geom_hline(yintercept = qchisq(c(.9, .95, .99), df=1),
             col="orange", linetype=2) +
  geom_vline(xintercept = edge, col="red")