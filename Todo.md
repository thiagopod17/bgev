


# TO DO LIST FOR FUTURE VERSIONS


- find a more well behaved integrate function or put a warning in the output when some of it does not converge so that the user interpret correctly the test results. 
- start values for mle from estimation of unimodal GEV and choose delta = 0. 



## correcao limites bgev.mle 

os limites da função bgev.mle não estão sendo atualizados quando vc troca na função, acho que é porque a linha que define o lower e upper se repete dentro da função ai quando troca essa parte ela não atualiza e mantém os mesmos limites

eu deixei essas duas linhas como comentario e funcionou normal
# FUNCTION:
# get reasonable starting values using
a genetic algorithm

#lower = c(-3,0.1,-3,-0.9)
#upper = c(30,30,3,3)


A função likbgev está devolvendo o valor errado, deveria ser
logl <- sum(log(dbgevd))

Eu testei com uma amostra muito grande e alguns valores iam para menos infinito, então eu removi os valores que o log iam para infinito ou menos infinito da soma final
logl <- sum(log(dbgevd)*is.finite(log(dbgevd)),na.rm=TRUE) #inf times 0 = NaN
  
