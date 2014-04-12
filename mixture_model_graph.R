
#
# GENERATION DU GRAPHE 
#

library(igraph);
library(mclust)


# n=10 points  and K=2 classes
K <- 2
n <- 10 
probas_connections_vertices <- matrix(NA,K,K)

for(k in 1:K)
{
  for(l in k:K)
  {
    if(l==k)
    {
      # elements of same group are very likely connected
      probas_connections_vertices[k,l] <- 0.9 
    }
    else
    {
      # elements of different groups are unlikely connected
      probas_connections_vertices[k,l] <- runif(1)/10 
      # the graph is not oriented
      probas_connections_vertices[l,k] <- probas_connections_vertices[k,l]
    }
  }
}


groupe <- matrix(NA,n,2)
A <- matrix(0,n,n)

# each row of groupe represents a vertice, each vertice is labelled by its number i and the groupe gr the belongs to.
for(i in 1:n)
{
  gr <- trunc(runif(1)*K) + 1
  groupe[i,1] <- i
  groupe[i,2] <- gr
}

# A is the connectivity matrix. A[i,j]=1 means that vertices i and j are connected.
for(i in 1:(n-1))
{
  for(j in (i+1):n)
  {
    A[i,j] <- rbinom(1,1,probas_connections_vertices[groupe[i,2],groupe[j,2]])
    A[j,i]=A[i,j]
  }
}

# we use igraph library to plot the graph
igraphX = graph.adjacency(A,mode=c("undirected"))
mode(igraphX)
str(igraphX)
plot(igraphX)



#
# FONCTIONS 
#


mise_a_jour <- function(i,k,gamma,p,pi,A)
{ 
  temp <- 1
  for (j in 1:nrow(A))
  {
    if(i!=j)
    {
      for(l in 1:ncol(gamma))
      {
        temp <- temp*(p[k,l]^A[i,j]*(1-p[k,l])^(1-A[i,j]))^gamma[j,l]
      }
    }
  }  
  return (pi[k]*temp)
}

mise_a_jour_P <- function(gamma,A)
{
  K <- ncol(gamma)
  p <- matrix(0,K,K)
  for (k in 1:K)
  {
    for(l in 1:K)
    {
      temp1 <- 0
      temp2 <- 0
      for(i in 1:nrow(A))
      {
        for(j in 1:nrow(A))
        {
          temp1 <- temp1 + gamma[i,k]*gamma[j,l]*A[i,j]
          temp2 <- temp2 + gamma[i,k]*gamma[j,l]
        }
      }
      p[k,l]=temp1/temp2
    }
  }
  return(p)
}

norme <- function(matrix){return(max(abs(matrix)))}

normalisation <- function(matrix)
{
  for(i in 1:nrow(matrix))
  {
    matrix[i,] <- matrix[i,]/sum(matrix[i,])
  }
  return(matrix)
}


Initialisation <- function(A,K)
{
  km <- kmeans(A,centers=K)  
  p <- matrix(0,K,K); 
  pi <- rep(NA,K);gamma=matrix(0,nrow(A),K)
  for (k in 1:K)
  {
    pi[k] <- length(which(km$cluster==k))/nrow(A)
  }
  for(i in 1:nrow(A))
  {
    for(k in 1:K)
    {
      if(km$cluster[i]==k)
      {
        gamma[i,k] <- 1
      }
    }
  }
  for (k in 1:K)
  {
    for(l in 1:K)
    {
      if(k==l)
      {
        n_prime=nrow(A[which(km$cluster==k),which(km$cluster==l)])
        if(n_prime==0)
        {
          p[k,l]=0
        } 
        else
        {
          p[k,l] <- sum(A[which(km$cluster==k),which(km$cluster==l)])/(n_prime*(n_prime-1))
        }
        #p[k,l]=p[l,k]
      }
      else
      {
        if(ncol(A[which(km$cluster==k),which(km$cluster==l)])*nrow(A[which(km$cluster==k),which(km$cluster==l)])==0)
        {
          p[k,l] <- 0
        }
        else
        {
          p[k,l] <- sum(A[which(km$cluster==k),which(km$cluster==l)])/(ncol(A[which(km$cluster==k),which(km$cluster==l)])*nrow(A[which(km$cluster==k),which(km$cluster==l)]))
        }
      }
    }             
  }
  newList=list("pi"=pi , "p"=p, "gamma"=gamma)
  return(newList)
}



#
# ALGORITHME EM
#

compteur <- 0
compteurE <- 0
epsilon <- 0.01
epsilonE <- 0.01
condition1 <- 1
#attention pas forcemment le meme K que celui d'avant
K=2

#initialisation de l'algo avec kmeans
initialisation <- Initialisation(A,K)
pi <- initialisation$pi
gamma <- initialisation$gamma
p <- initialisation$p 

while(condition1>epsilon)
{
  compteur <- compteur+1
  
  # phase E (Excpectation)
  gamma_avant_M <- gamma
  #gamma_avant_E=gamma+1
  condition2 <- 1
  while(condition2>epsilonE)
  {
    gamma_avant_E <- gamma;
    for(i in 1:n)
    {
      for(k in 1:K)
      {
        gamma[i,k] <- mise_a_jour(i,k,gamma_avant_E,p,pi,A)
      }
    }
    gamma <- normalisation(gamma)
    condition2=norme(gamma-gamma_avant_E)}
    compteurE=compteurE+1
    
    #phase M (Maximization)
    pi_avant <- pi
    p_avant <- p
    for(k in 1:K)
    {
      pi[k] <- sum(gamma[,k])/n
    }
    p <- mise_a_jour_P(gamma,A)
    #stabilisation de pi (proportions) et p (les probas d'interconnexion)
    condition1 <- max(norme(pi-pi_avant),norme(p-p_avant)) 
}


#
# TRACE DU RESULTAT
#

# autant de couleurs que de classes K => g?n?ration al?atoire de couleurs
liste_couleurs = rep(NA,K)
liste_couleurs[1]="red";
liste_couleurs[2]="blue";
liste_couleurs[3]="green";
liste_couleurs[4]="yellow"; 

groupe_posteriori = matrix(NA,n,2)

for(i in 1:n)
{
  numero_classe_i <- which(max(gamma[i,])==gamma[i,])
  groupe_posteriori[i,1] <- i
  groupe_posteriori[i,2] <- liste_couleurs[numero_classe_i]
}

couleurs_sommets <- as.vector(groupe_posteriori[,2])
igraphX <- set.vertex.attribute(igraphX, "color", value=couleurs_sommets)
plot(igraphX)





#------------------------------------------------------------------------------------------------------------------------------------------------------------------

#
# LIBRAIRY MIXER APPLIED TO SIMPLE EXAMPLE
#

library(igraph)
library(mixer)
library(RCurl)
library(rjson)


PI <- matrix(NA,2,2);
PI[1,1] <- 0.85
PI[1,2] <- 0.05
PI[2,1] <- 0.05
PI[2,2] <- 0.9

#equiprobability to be in one of the two groups, here the two groups have the same size
alpha1 <- 0,5
alpha2 <- 0,5

n <- 30
groupe <- matrix(NA,n,2)
X <- matrix(0,n,n)

for(i in 1:n)
{
  if(runif(1)<0.5)
  {
    groupe[i,2] <- 1
    groupe[i,1] <- i
  } 
  else
  {
    groupe[i,2] <- 2
    groupe[i,1] <- i
  }
}

for (i in 1:(n-1))
{
  for (j in (i+1):n)
  {
    X[i,j] <- rbinom(1,1,PI[groupe[i,2],groupe[j,2]])
    X[j,i] <- X[i,j]
  }
}

igraphX = graph.adjacency(X,mode=c("undirected"))
mode(igraphX)
str(igraphX)
plot(igraphX) 
Mix <- mixer(X) 
taux = Mix$output[[1]]$Taus
groupePosteriori <- matrix(NA,n,2)

for (i in 1:n)
{
  if(taux[1,i]>taux[2,i])
  {
    groupePosteriori[i,2]= "red"
    groupePosteriori[i,1]= i
  } 
  else
  {
    groupePosteriori[i,2] <- "yellow"
    groupePosteriori[i,1] <- i
  }  
}

listecouleurs = as.vector(groupePosteriori[,2])
igraphX <- set.vertex.attribute(igraphX, "color", value=listecouleurs);
plot(igraphX)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

#
# FACEBOOK PROJECT  
#

#http://tuxette.nathalievilla.org/?p=953&lang=en
#http://blog.revolutionanalytics.com/2012/01/visualize-your-facebook-friends-network-with-r.html
   
library(igraph)
library(mixer)
library(RCurl)
library(rjson)

token <-  "yourtoken" # paste you 
data <-  getURL(sprintf( "https://graph.facebook.com/%s&access_token=%s", "764127769?fields=friends.fields(name,mutualfriends)", token ),ssl.verifypeer = FALSE)

friends  <- fromJSON(data);
n <- length(friends$friends$data);
nbr_id_name <- matrix(NA,n,3);
X <- matrix(0,n,n)

#friends.name = sapply(friends$friends$data, function(x) iconv(x$name,"UTF-8","ASCII//TRANSLIT"))
for (i in 1:n)
{
  nbr_id_name[i,1] <- i
  nbr_id_name[i,2] <- friends$friends$data[[i]]$id
  nbr_id_name[i,3] <- friends$friends$data[[i]]$name
}

for (i in 1:n)
{
  mutual_friends_i <- friends$friends$data[[i]]$mutualfriends$data
  nbr_mutual_friends_i <- length(friends$friends$data[[i]]$mutualfriends$data)
  for (j in 1:nbr_mutual_friends_i)
  {
    X[i,which(nbr_id_name[,2]==friends$friends$data[[i]]$mutualfriends$data[[j]]$id)]=1
  }
}  

igraphX <- graph.adjacency(X,mode=c("undirected"))
mode(igraphX);#str(igraphX)
plot(igraphX,vertex.size=3)
Mix <- mixer(X,qmin=1,qmax=10)
plot(Mix, frame=5)
nbr_classes <- 3
liste_couleurs <- rep(NA,nbr_classes)
liste_couleurs_Font <- rep(NA,nbr_classes)

# for(j in 1:nbr_classes){liste_couleurs[j]=col2rgb(j)}
#liste_couleurs[1]="red";liste_couleurs[2]="blue";liste_couleurs[3]="green";#liste_couleurs[4]="orange";liste_couleurs[5]="purple";liste_couleurs[6]="pink";liste_couleurs[7]="orange";liste_couleurs[8]="black";

liste_couleurs[1] <- "red"
liste_couleurs[2] <- "grey52"
liste_couleurs[3] <- "green"
liste_couleurs_Font[1]="red4";
liste_couleurs_Font[2]="black";
liste_couleurs_Font[3]="green4";

groupePosteriori <- matrix(NA,n,3)

for (i in 1:n)
{
  if(max(Mix$output[[nbr_classes]]$Taus[,i])==0)
  {
    groupePosteriori[i,1]=i;
    groupePosteriori[i,2]="grey";
    groupePosteriori[i,3]="grey"
  }
  else
  {
    numero_classe_i <- which(Mix$output[[nbr_classes]]$Taus[,i]==max(Mix$output[[nbr_classes]]$Taus[,i]))
    groupePosteriori[i,1] <- i
    groupePosteriori[i,2] <- liste_couleurs[numero_classe_i]
    groupePosteriori[i,3] <- liste_couleurs_Font[numero_classe_i]
  }
}

vertices_label = rep("",n) #nbr_id_name[,3]

for(i in 2:n)
{
  if(i%%2==0)
  {
    vertices_label[i]=nbr_id_name[i,3]
  }
}

couleurs_sommets <- as.vector(groupePosteriori[,2])
couleurs_Font <- as.vector(groupePosteriori[,3])
igraphX <- set.vertex.attribute(igraphX, "color", value=couleurs_sommets)
plot(igraphX,vertex.size=1,vertex.label=vertices_label,vertex.label.color=couleurs_Font, vertex.label.family=c("Calibri"),vertex.label.font=2,vertex.label.cex=0.5,edge.width=0.1,asp=0.5)




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
