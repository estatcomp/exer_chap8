

install.packages("FNN")

install.packages("data.table")

install.packages("mltools")
 
install.packages("latticeExtra")




library(boot)
library(FNN)
library(latticeExtra)
library(data.table)
library(mltools)



# ___________
# 
#   8.1
# ___________


#Implement the two-sample Cramer-von Mises test for equal distributions as a
#permutation test. Apply the test to the data in Examples 8.1 and 8.2.





# 1) Definir x e y


attach(chickwts)

x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))

detach(chickwts)


  #vetor unindo as amostras

z <- (c(x,y))



# 2) Fazer um teste T para x=y (name = "t0")


t0 <- t.test(x,y)$statistics




# 4) Gerar uma estatística alternativa "Cramer - Von - Mises"


  # tamanho das amostras

m = length(x)


n = length(y)




  # definind ecdf(G,F (x,y))


  g <- ecdf(y) #dist.empírica y (Gn)
  
  
  gx <- g(x)  #Gn(x)
  
  
  gy <- g(y)  #Gn(y)

  
  
  f <- ecdf(x) #dist.empírica de x
  
  fx <- f(x)  # Fn(x)
  
  fy <- f(y)  #Fn(y)
  
  

  
  #função

  for(n in 1: length (x)){
    
    
    
    for (m in 1:length(y)){
      
      
      
      w2=m*n/(m+n)^2*(sum(fx[n]-gx[n])^2 + sum(fy[m]-gy[m])^2)
      
      
    }

    w2   #valor exato de    
     
  }    
    




# 5) Gerar um p-valor = mean (c(t0,C-M)>= t0)




p <- mean(c(t0,w2))

p






# 8.1 Método alternativo


#definindo parâmetros
  
  
R <- 999 
z <- c(x, y) 
K <- 1:26
D <- numeric(R) 
options(warn = -1)


#tamanhos da amostra

m = length(x)

n = length(y)



#valor exato


for(n in 1: length (x)){
  
  
  
  for (m in 1:length(y)){
    
    
    
    w2=m*n/(m+n)^2*(sum(fx[n]-gx[n])^2 + sum(fy[m]-gy[m])^2)
    
    
  }
  
  w2   #valor exato da estatística N-C  
  
}    




# algorítimo de geração de nova varável com dist.conjunta

for (i in 1:R) {
  
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  w2_2[i] <- ks.test(x1, y1, exact = FALSE)$statistic
  
  w2_2
  
}



#cálculo do p-valor

p2 <- mean(c(w2_2,w2) >= w2)


p2







#________________________________

#             8.2
#________________________________



# Implement the bivariate Spearman rank correlation test for independence
# [255] as a permutation test. The Spearman rank correlation test statistic can be obtained from function cor with method = "spearman". Compare the
# achieved significance level of the permutation test with the p-value reported
# by cor.test on the same samples.



# p-valor empírico

R = 1000

n = 15

x <- runif(n,10,20)

y <- runif(n, 10,20)

cor.test(x,y,method = "spearman")


cor_empirica <-cor.test(x,y,method = "spearman")$statistic

cor_empirica


#p-valor da permutação

z <- c(x,y)

K <- 1:30


corpermut <- numeric(R)


  # loop

for (i in 1:R){
  
  k <- sample(K, size = 15, replace = FALSE)
  
  x1 <- z[k]
  
  y1 <- z[-k]
  
  corpermut[i] <- cor.test(x1, y1)$statistic

}



p <- mean(c(cor_empirica,corpermut) >= cor_empirica)

p # reduz-se muito o p-valor.




#________________

#       8.3

#________________


#definindo variáveis

x1<- rnorm(25,0,2)

x2 <- rnorm(20,0,1)



#estatística de ~variância

X1 <- x1-mean(x1)

X2 <- x2-mean(x2)



#função valor exato



func5 <- function(x1,x2){
  
  
  X1 <- x1-mean(x1)
  
  X2 <- x2-mean(x2)
  
  
  # contanto extremos
  
  out1 <- sum(X1 > max(X2)) + sum(X1 < min(X2))
  
  out2 <- sum(X2 > max(X1)) + sum(X2 < min(X1))
  
  return(max(c(out1, out2) > 5))
  
}



max5 <- func5(x1,x2) # valor exato



# variáveis do algoritimo de permutação


n <- 1000 # amostragens


z <- c(x1,x2) # vetor com as duas obs


nz <-length(z)




for (i in 1:n){
  
  cont5permut <- numeric(r) #vetor vazio
  
  k = sample(length(z), size=length(z),replace = FALSE) # amostragens
  
  cont5permut[i] = func5(z[k],z[-k]) # função nas reamostragens
  
  }



p5 <-mean( c(max5,cont5permut) > max5 ) # p- valor da permutação
  
p5  





#___________________

#_______8.4_________

# __________________




data("chickwts")

#puxando os dados
x <- with(chickwts, as.vector(weight[feed == "sunflower"]))
y <- with(chickwts, as.vector(weight[feed == "linseed"]))
z <- cbind( c(x, y), rep(0, length(c(x, y))) )



#função de teste knn

my.knn <- function(z, nn) {
  t_n.i <- function(z, nn, i) {
    n_g <- nrow(z) / 2 # length x and length y
    n <- nrow(z)
    z <- z[i, ]
    z <- cbind(z, rep(0, n))
    knn <- get.knn(z, k = nn) # Obtaining the k nearest neighbors
    n_1 <- knn$nn.index[1:n_g, ]       # Dividing x
    n_2 <- knn$nn.index[(n_g + 1):n, ] # and y
    i_1 <- sum(n_1 < n_g + .5) # Obtaining the sum of
    i_2 <- sum(n_2 > n_g + .5) # the indicator functions
    return( (i_1 + i_2) / (nn * n) ) # Test statistic
  }
  perm <- boot(data = z # 10000 permutation samples
               , statistic = t_n.i
               , sim = "permutation"
               , R = 10000
               , nn = 3)
  p <- with(perm, mean( c(t, t0) >= t0 )) # p-value
 
  
   return(
    
    
    
    histogram(with(perm, c(t, t0)) # Histogram
              , type = "density"
              , col = "#0080ff"
              , xlim = c(-.1, 1.1)
              , xlab = paste0("Replicates of T(n, ", nn, ") statistic")
              , ylab = list(rot = 0)
              , main = paste0("Permutation distribution of T(n, "
                              , nn, ") statistic")
              , sub = list(substitute(paste(hat(p), " = ", pvalue)
                                      , list(pvalue = p))
                           , col = 2)
              , panel = function(...){
                panel.histogram(...)
                panel.abline(v = p, col = 2, lwd = 2)
              })
  )
}
my.knn(z, 3) # k (or r) nearest neighbors = 3


  
  




