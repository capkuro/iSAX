#########################################################################################################
# iSAX is an R package which provides access to iSA technology developed by 
# VOICES from the Blogs. It is released for academic use only and licensed 
# under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License 
# see http://creativecommons.org/licenses/by-nc-nd/4.0/
# Warning: Commercial use of iSA is protected under the U.S. provisional patent application No. 62/215264
#########################################################################################################

iSA <-
function(Strain, Stest, Dtrain, nboot=1000, predict=FALSE, ret.boot=FALSE, seqlen=5,
   sparse=FALSE, verbose=TRUE, tolerance=0){
    ptm <- proc.time()
    #  require(quadprog)
    #require(data.table)
    
	#Funcion para crear la tabla de probabilididades
    ptable <- function(x)
    {
        a <- data.table(x)
        b <-  a[, .N, by=x]
        tmp <- b$N
        names(tmp) <- b$x
        return(tmp/sum(tmp))
    }
    
    if(verbose)
    cat("\niSAX...\n")
    
    #S Vector de stems 
    S <- c(Strain,Stest)
    
	#D Posibles Categorias (Categorias training + Categorias NA testing)
    D <- c(Dtrain, rep(NA, length(Stest)))
    
	#Se obtiene el tamaño del feature string convertido a hexadecimal
    nc <- nchar(S[1])
    if(seqlen>0){
		#Asignacion de seqlen según las caracteristicas de feature string
        if(nc<=3)
            seqlen <- 1
        if(nc>10 & seqlen<3)
            seqlen <- 3
        if(seqlen>=nc)
            seqlen <- nc
		
		#Obtencion de numero de sequencias en total para hacer el split de iSAX
        nseq <- floor(nc/seqlen)
		#se genera un arreglo con la suma acumulada de seqlen, repetido nseq veces
        splits <- cumsum(rep(seqlen,nseq))
		#N° de splits
        nsplits <- length(splits)
		#concatenando nc a splits, o en su defecto asignandolo al ultimo elemento del arreglo splits
        if(splits[nsplits] > nc-2){
            splits[nsplits] <- nc
        } else {
            splits <- c(splits,nc)
        }
		#Nuevo N° de splits
        nsplits <- length(splits)
		#N° de textos a clasificar
        lS <- length(S)
		#nuevos stems seccionados
        newS <- character(nsplits * lS)
		#Si existen splits
        if(nsplits>2){
			#Secciona S desde el indice 1 hasta splits[1], añadiendole la letra 'a' al principio
            newS[1:lS] <- paste0(letters[1],substr(S,1,splits[1]))
			#Secciona el resto de s, desde splits[i] hasta splits[i+1], asignandole la letra numero i+1 ("b","c",...,etc)
            for(i in 1:(nsplits-1)){
                newS[(lS*i+1):(lS*(i+1))] <- paste0(letters[i+1],substr(S,splits[i]+1,splits[i+1]))
            }
			#Resamplea D, asignando a newD con las repeticiones (resampleo de newS se le asignan la misma categoria que el S original)
            newD <- rep(D, nsplits)
            D <- newD
            S <- newS
        }
    }
    #Indices de elementos de testing
    idd <- which(is.na(D))
	#Separacion de S y D para 
    Strain <- S[-idd]
    Dtrain <- D[-idd]
    
    # distribuzione di D nel training set
	#Se obtienen la distribucion de probabilidad de D (ver funcion en linea 17)
    pD.train <- ptable(Dtrain)

	
    # distribuzione complessiva di S
	#se obtiene la distribucion total de las stems S (seccionadas(newS))
    pS <- ptable(S)
    nS <- names(pS) # stem names
    lS <- length(pS) # number of unique stems
    
	#revisar el else de este if, debe de realizar lo mismo pero para matrices sparse
    if(sparse){
     prSD <- xtabs(~ Strain + Dtrain, sparse=sparse)
     sums <- apply(prSD,2,sum)
     for(i in 1:ncol(prSD))
      prSD[,i] <- prSD[,i]/sums[i]
    } else {
		#prop.table <- "Express Table Entries As Fraction Of Marginal Table"
     prSD <- prop.table(table(Strain,Dtrain),2)
    }
	
	#se obtiene el numero de columnas en conjunto el nombre de las filas y las columnas
    lD <- ncol(prSD)
    nD <- colnames(prSD)
    nSD <- rownames(prSD)
    
	#se genera una matriz de lS (revisar linea 93) x lD (linea 107)
    P <- matrix(0,  nrow=lS, ncol=lD) #, dimnames=list(nS,nD))
	#se obtiene los indices de los match entre nSD (linea 109) y nS(linea 92)
    idx <- match(nSD, nS )
	#Asigna a los indices idx, los valores obtenidos en prSD (linea 103)
    P[idx, ] <- as.numeric(prSD)
	#Matriz con resultados de Error estandar, z value y confidencia
    tb <- matrix(, lD, 4)
    rownames(tb) <- nD
    colnames(tb) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    #Matriz temporal
    tmp <- matrix(,lD,1)
    colnames(tmp) <- "iSA"
    rownames(tmp) <- nD
    
    
    # contrained iSA
	#Numero de columnas de P (linea 112)
    q <- ncol(P)
    ind <- matrix(0, q, 2)
    rownames(ind) <- nD
    ind[, 1] <- 1
    q0 <- ncol(P)
	#matrix defining the constraints under which we want to minimize the quadratic function.
    Amat <- matrix(0, q0, q0 * 2 + 1)
	#Asignar valores 1 a la primera columna
    Amat[, 1] <- rep(1, q0)
	#asigna para cada columna de 2 a (q0+1) el valor de una matriz identidad
    Amat[, 2:(q0 + 1)] <- diag(1, q0)
	#asignar para cada columna de (q0+2):(2*q0+1) el valor de una matriz con diagonal -1
    Amat[, (q0 + 2):(2 * q0 + 1)] <- diag(-1, q0)
	#vector lleno de 0 de largo de q(linea 129)
    lbd <- rep(0, q)
	#vector lleno de 1 de largo de q(linea 129)
    ubd <- rep(1, q)
    const <- 1
	
	#bvec vector holding the values of b_0 (defaults to zero). <- para resolver la programacion quadratica, ver linea 156
    bvec <- c(const, lbd[ind[, 1] == 1], -ubd[ind[, 1] == 1])
	# se calcula el absoluto del determinante del resultado de la multiplicacion la transpuesta de P (t(P)) y P 
    adetP <- abs(det(t(P)%*% P))
    if(verbose)
     cat(sprintf("\nNote: abs(det(P'*P))=%g\n", adetP))
	#checkeo de si se puede realizar la inversa de la matriz
    if(adetP < tolerance)
      stop("Matrix P'*P is not invertible")

	#Intentar encontrar solucion al problema de programacion cuadratica con:
	# t(P)%*% P <- matriz que aparece en la funcion cuadratica a minimizar
	# t(pS) <- vector que aparece en la funcion cuadratica a minimizar
	# Amat <- ver linea 137:141
	# bvec <- ver linea 149
	# meq  <- the first meq constraints are treated as equality constraints, all further as inequality constraints (defaults to 0).
    aa <- try(solve.QP(t(P)%*% P, t(pS) %*% P, Amat,bvec, meq = 1), TRUE)
	#si falla repite NA q0 veces
    if(class(aa) == "try-error"){
        b <- rep(NA, q0)
    } else {
	#sino, b se le asigna el vector solucion el cual contiene los valores de solucion para el problema de programacion cuadratica
        b <- aa$solution
    }
	#se calcula la varianza
    sigma2 <- sum((pS - as.numeric(P %*% b))^2)/(length(pS)-ncol(P)-1)
    qp <- tmp
	#a la primera columna de qp se le asigna los valores de la solucion
    qp[,1] <- b
    colnames(qp) <- "iSA"
	#tabc se le asigna la matriz de error estandar, z value y confidencia (linea 118)
    tabc <- tb
	# se asigna la primera columna de qp a tabc
    tabc[,1] <- qp[,1]
	#Error estandar: se calcula la raiz de la multiplicacion de la varianza por la diagonal (diag) de la inversa (solve) de la multiplicación entre la transpuesta (t) de P y P
    serr <- try(sqrt(sigma2 * diag(solve(t(P)%*% P))), TRUE)
    if(class(serr)=="try-error")
    serr <- rep(NA, q0)
    tabc[,2] <- serr # Error estandar
    tabc[,3] <- tabc[,1]/tabc[,2] #Valor Z
    tabc[,4] <- pnorm(tabc[,3],lower.tail=FALSE)*2 # Norma P
    boot <- NULL
    
    if(nboot>0){
		#bootstraping <- https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
        if(verbose)
        cat("\nbootstrapping...please wait")
        for(i in 1:nboot){
			#se obtiene length(pS) numeros al azar (con repeticion) entre 1 y length(pS), los cuales se utilizaran como indices para las filas
            idx <- sample(1:length(pS), length(pS), replace=TRUE)
            #se asigna a tP las filas de P los cuales esten contenido en P
			tP <- P[idx,]
			#vector de probabilidad de S (linea 91) para la nueva muestra
            tpS <- pS[idx]
			#Intentar encontrar solucion al problema de programacion cuadratica con:
			# t(tP)%*% tP <- matriz que aparece en la funcion cuadratica a minimizar
			# t(tpS)%*% tP <- vector que aparece en la funcion cuadratica a minimizar
			# Amat <- ver linea 137:141
			# bvec <- ver linea 149
			# meq  <- the first meq constraints are treated as equality constraints, all further as inequality constraints (defaults to 0).
            aa <- try(solve.QP(t(tP)%*% tP, t(tpS) %*% tP, Amat,bvec, meq = 1),TRUE)
			# se genera un vector el cual se les asignara las soluciones obtenidas en aa
            b <- rep(NA, q0)
            if(class(aa) != "try-error"){
                b <-aa$solution
            }
			#append a la matriz boot una fila nueva añadiendole b
            boot <- rbind(boot, b)
        }
        
        if(verbose)
        cat("\n")

        b.cf <- tmp
		#a la primera columna se le asigna las medias de las columnas
        b.cf[,1] <- colMeans(boot,na.rm=TRUE)
        colnames(b.cf) <- "iSAb"
		#aplica la funcion de desviacion de estandar a los extremos de boot 
		#apply() <- Returns a vector or array or list of values obtained by applying a function to margins of an array or matrix. 
        b.sd <- apply(boot,2, function(u) sd(u, na.rm=TRUE))
        #print(b.sd)
		#b.t corresponde al coeficiente de variacion (mu/sigma)
        b.t <- NA
        if(!is.na(sum(b.sd)))
        b.t <- b.cf/b.sd
        
        #se le asigna la matriz de error estandar, z value y confidencia (linea 118)
        tabb <- tb
		#se remplaza la primera columna con el vector de medias
        tabb[,1] <- b.cf
		#se remplaza la segunda columa con el vector de sigma
        tabb[,2] <- b.sd
		#se reemplaza la tercera columna  con el vector de coeficiente de variacion
        tabb[,3] <- b.t
		#se calcula la norma P y se coloca en la columna 4
        tabb[,4] <- pnorm(tabb[,3],lower.tail=FALSE)*2
    } else {
	#Si no hay bootstraping, tabb se le asigna tabc y a la solucion b se le añade qp (linea 174)
        tabb <- tabc
        b.cf <- qp 
    }
    #Se define la funcion predict
    predict.iSA <- function(x){
        #se busca el indice del elemento X en nS
        idx <- match(x, nS)
		#se aplica a entre los valores 1:lD (linea 107) la formula P[i,x]*b.cf[x]/pS[i], para obtener el valor mayor para cada uno de los indices idx obtenidos.
        sapply(idx, function(i) which.max(sapply(1:lD,function(x) P[i,x]*b.cf[x]/pS[i]))) -> aa
		#se guarda el resultado anterior
        nD[aa]
    }
    
    pred <- NULL
	#se obtiene las predicciones
    if(predict){
        pred <- predict.iSA(S)
    }
	#se obtiene el tiempo
    etime <- (proc.time()-ptm)[1]
    if(verbose)
    cat(sprintf("\nElapsed time: %.2f seconds\n",etime))
	#si ret.boot, se debuelve el bootstraping, de lo contrario se retorna null en ese argumento.
    if(ret.boot)
     return(list( est=qp, tab=tabc, best=b.cf, btab=tabb, boot=boot, pred=pred,time=etime ) )
    return(list( est=qp, tab=tabc, best=b.cf, btab=tabb, boot=NULL, pred=pred,time=etime ) )
}
