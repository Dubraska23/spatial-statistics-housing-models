
#Paquetes necesarios para el seguimiento del temario

#install.packages("raster")
#install.packages("dplyr")
#install.packages("sf")
#install.packages("ggplot2")
#install.packages("raster")
#install.packages("osmdata")
#install.packages("leaflet")
#install.packages("terra")  
#install.packages("spdep")  
#install.packages("knitr")
#install.packages("kableExtra")
#install.packages("spatstat")
#install.packages("geosphere")
#install.packages("data.table")
#install.packages("dbscan")
#install.packages("osrm")
#install.packages("SpatialEpi")
#install.packages("spatialreg")
#install.packages("earth")
#install.packages("spgwr")
#install.packages("pspatreg")
#install.packages("rpart")
#install.packages("SpatialML")
#install.packages("cluster")
#install.packages("mlr3")
#install.packages("mlr3learners")
#install.packages("mlr3spatial")


library(raster)
library(dplyr)
library(sf)
library(ggplot2)
library(raster)
library(osmdata)
library(leaflet)
library(terra)  
library(spdep)  
library(knitr)
library(kableExtra)
library(spatstat)
library(geosphere)
library(data.table)
library(dbscan)
library(osrm)
library(tidyr)
library(SpatialEpi)
library(spatialreg)
library(earth)
library(spgwr)
library(pspatreg)
library(rpart)
library(SpatialML) 
library(cluster) 
library(mlr3)
library(mlr3learners)
library(mlr3spatial)



eliminar_cercanos <- function(data, umbral_metros) {
  # Inicializar un vector para almacenar los índices de los puntos a conservar
  indices_a_conservar <- c(TRUE)  # Comenzamos con el primer punto
  
  # Iterar sobre cada fila del dataframe
  for (i in 1:nrow(data)) {
    # Si el índice ya está marcado como conservado, continuar
    if (!indices_a_conservar[i]) next
    
    # Calcular distancias desde el punto actual a los puntos restantes
    distancias <- distHaversine(c(data$lon[i], data$lat[i]), 
                                cbind(data$lon, data$lat))
    
    # Marcar los índices que están dentro del umbral (menos de 100 metros)
    indices_cercanos <- distancias < umbral_metros
    
    # Solo conservar el primer índice (el actual) y marcar los demás como FALSE
    indices_a_conservar <- indices_a_conservar | (1:nrow(data) %in% which(indices_cercanos) & (1:nrow(data) != i))
  }
  
  # Devolver el dataframe sin los puntos cercanos
  return(data[which(indices_a_conservar), ])
}

#Función Histogramas
Hist<-function(bbdd_fff,response,predicted,var,n,breaks=8){
  
  
  
  names<-colnames(bbdd_fff)
  
  bbdd_fff$predicted<-predicted
  
  bbdd_fff$response<-response
  
  bbdd_fff$VIVO<-1
  
  q1<- bbdd_fff %>%
    
    group_by(cut(var,breaks=breaks)) %>%
    
    summarise(Exposicion = sum(VIVO),
              
              Frecuencia = sum(response)/sum(VIVO),
              
              Predicted = sum(predicted)/sum(VIVO)
              
    )
  
  
  
  q1<-as.data.table(q1)
  
  
  
  c1<-q1[,1]
  
  q1$c1<-q1[,1]
  
  
  
  ff<-ggplot(q1,aes(x=c1)) +
    
    geom_bar(aes(y=(Exposicion/sum(Exposicion))*mean(Frecuencia)*2, fill="% de Exposicion"), stat = "identity")+
    
    geom_point(aes(y=Frecuencia, colour="Frecuencia"), group=1)+
    
    geom_line(aes(y=Frecuencia, colour="Frecuencia"), group=1)+
    
    geom_point(aes(y=Predicted, colour="Predicted"), group=1)+
    
    geom_line(aes(y=Predicted, colour="Predicted"), group=1)+    
    
    xlab("var") + ylab("") +
    
    #scale_y_continuous(limits = c(0, 0.45),breaks = c(0:45/100))+
    
    theme(legend.key=element_blank(),legend.title=element_blank(),
          
          legend.box="horizontal", legend.position = "bottom",
          
          axis.text.x = element_text(angle = 60)) +
    
    ggtitle(names[n])
  
  
  
  return(ff)
  
  
  
  
  
}

#Función Kulldorf
kulldorff_scan<-function(X,response,EXP){
  
  #INICIO
  
  X$response<-response
  X$response_exp<-EXP
  
  # Crear Polígonos Espaciales
  ###############################
  df_sf <- st_as_sf(X, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  bbox <- st_bbox(df_sf)
  bbox_sf <- st_as_sfc(bbox)
  cellsize <- 0.01 
  grid_hex <- st_make_grid(bbox_sf, cellsize = cellsize, square = FALSE) 
  grid_hex_sf <- st_sf(geometry = grid_hex)
  grid_hex_sp <- as(grid_hex_sf, "Spatial")
  data_hex <- data.frame(ID = 1:length(grid_hex_sp), row.names = as.character(1:length(grid_hex_sp)))
  Geo <- SpatialPolygonsDataFrame(grid_hex_sp, data = data_hex)
  
  # Verificar si Geo ya tiene un CRS
  if (is.na(proj4string(Geo))) {
    proj4string(Geo) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Solo si no tiene CRS
  } else {
    Geo <- spTransform(Geo, CRS("+proj=longlat +datum=WGS84 +no_defs"))  # Reproyectar correctamente
  }
  # Agregar ID a los datos
  Geo@data <- cbind(ID = 1:nrow(Geo@data), Geo@data)
  
  Geo_SC_Ag<-Geo
  
  # Puntos Espaciales
  ###############################
  bbdd2 <- X
  
  # Convertir a SpatialPoints
  coordinates(bbdd2) <- c("LONGITUDE", "LATITUDE")
  
  # Verificar si ya tiene un CRS asignado
  if (is.na(proj4string(bbdd2))) {
    proj4string(bbdd2) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Solo si no tiene CRS
  } else {
    bbdd2 <- spTransform(bbdd2, CRS("+proj=longlat +datum=WGS84 +no_defs"))  # Reproyectar correctamente
  }
  
  bbdd33<-bbdd2
  bbdd33$VIVO<-1
  
  
  #Agregamos Espacialmente Polígonos
  ###############################
  agrss<-over(Geo,bbdd33,fn = sum)
  
  agrss$VIVOS<-agrss$VIVO
  agrss$responses<-agrss$response/agrss$VIVO
  agrss$response_exp<-agrss$response_exp/agrss$VIVO
  agrss[is.na(agrss)] <- 0
  Geo_SC_Ag@data<-agrss
  
  Geo_SC_Ag<-Geo_SC_Ag[Geo_SC_Ag$VIVOS>0,]
  Geo_SC_Ag@data<-cbind(ID=1:nrow(Geo_SC_Ag@data),Geo_SC_Ag@data)
  row.names(Geo_SC_Ag@data) <- seq_len(nrow(Geo_SC_Ag@data))
  #Creamos Base de Datos
  ###############################
  LATITUDE<-centroid(Geo_SC_Ag)[,2]
  LONGITUDE<-centroid(Geo_SC_Ag)[,1]
  Geo_Y<-cbind(Geo_SC_Ag@data,LATITUDE,LONGITUDE)
  Geo_Y <- as.data.frame(sapply(Geo_Y, function(x) as.numeric(as.character(x))))
  
  ##Preparamos para Kulldorf
  ###############################
  geo2<-dplyr::select(Geo_Y,LONGITUDE,LATITUDE)
  geo2<-rename(geo2,x=LONGITUDE,y=LATITUDE)
  geo2 <- latlong2grid(geo2)
  
  ## Get aggregated counts of population and cases for each county
  population2 <- tapply(Geo_Y$VIVOS,Geo_Y$ID,sum)
  cases2 <- tapply(Geo_Y$responses,Geo_Y$ID,sum)
  expected.cases2 <- as.numeric(Geo_Y$response_exp)
  
  population2<-Geo_Y$VIVOS
  cases2<-Geo_Y$responses
  
  ## Set Parameters
  pop.upper.bound <- 0.6
  n.simulations <- 500
  alpha.level <- 0.1
  plot <- FALSE
  
  poisson <- SpatialEpi::kulldorff(geo = geo2, cases = cases2, population = population2,expected.cases =expected.cases2,pop.upper.bound= pop.upper.bound,n.simulations = n.simulations, alpha.level = alpha.level,plot = plot)
  
  
  cluster <- poisson$most.likely.cluster$location.IDs.included
  
  GR<-Geo_SC_Ag[cluster,]
  GR$clustering<-1
  agrss<-over(bbdd2,GR)$clustering
  
  
  lista<-list(Geo_SC_Ag,cluster,GR,agrss)
  
  return(lista)
  
}

.GWR_int_miguel <- function(fit.points, coords, gweight, y, x, weights, yhat,GWR_args) {
  
  if (GWR_args$predictions) {
    predx <- fit.points[, -c(1,2), drop=FALSE]
    fit.points <- fit.points[, c(1,2), drop=FALSE]
    # Maciej Kryza 130906 drop issue
  }
  
  n <- nrow(fit.points)
  m <- NCOL(x)
  x1 <- matrix(1, nrow=nrow(x), ncol=1)
  sum.w <- numeric(n)
  betas <- matrix(nrow=n, ncol=m)
  colnames(betas) <- colnames(x)
  if(!GWR_args$fp.given) {
    gwr.e <- numeric(n)
  } else {
    gwr.e <- NULL
  }
  if (GWR_args$se.fit) {
    betase <- matrix(nrow=n, ncol=m)
    colnames(betase) <- paste(colnames(x), "se", sep="_")
  } else {
    betase <- NULL
  }
  if (GWR_args$predictions) {
    pred <- numeric(n)
    if (GWR_args$se.fit) {
      pred.se <- numeric(n)
    } else {
      pred.se <- NULL
    }
  } else {
    pred <- NULL
    pred.se <- NULL
  }
  if (is.null(yhat)) {
    localR2 <- NULL
  } else {
    localR2 <- numeric(n)
  }
  
  if (!GWR_args$fp.given && GWR_args$hatmatrix) 
    lhat <- matrix(nrow=n, ncol=n)
  if (is.null(GWR_args$adapt)) {
    bw <- GWR_args$bandwidth
    bandwidth <- rep(GWR_args$bandwidth, n)
    
    bbbbbbbbbb <- GWR_args$distancias                                                                 # Ajusto el BandWith
    bandwidth<- GWR_args$distancias
  } else {
    bandwidth <- gw.adapt(dp=coords, fp=fit.points,quant=GWR_args$adapt, longlat=GWR_args$longlat)
    bw <- bandwidth
    
    bandwidth<-GWR_args$distancias
  }
  if (any(bandwidth < 0)) stop("Invalid bandwidth")
  for (i in 1:n) {
    dxs <- spDistsN1(coords, fit.points[i,], 
                     longlat=GWR_args$longlat)
    if (any(!is.finite(dxs)))
      dxs[which(!is.finite(dxs))] <- 0
    #		if (!is.finite(dxs[i])) dxs[i] <- 0
    w.i <- gweight(dxs^2, bandwidth[i])
    w.i <- w.i * weights
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- lm.wfit(x, y, w.i)
    sum.w[i] <- sum(w.i)
    betas[i,] <- coefficients(lm.i)
    ei <- residuals(lm.i)
    # prediction fitted values at fit point
    if (GWR_args$predictions) {
      pred[i] <- sum(predx[i,] * betas[i,])
    }
    # use of diag(w.i) dropped to avoid forming n by n matrix
    # bug report: Ilka Afonso Reis, July 2005
    
    # local R-squared as explained by Tomoki Nakaya, 090624 (in GWR4)
    # differs from local weighted regression R-squared
    
    if (!is.null(yhat)) {
      RSS <- sum(w.i * (y - yhat)^2)
      yss <- sum(w.i * (y - weighted.mean(y, w.i))^2)
      localR2[i] <- 1 - (RSS/yss)
    }
    
    if (GWR_args$se.fit) {
      p <- lm.i$rank
      if (p != m) {
        warning(paste("OLS fit not full rank at fit point", i))
      } else {
        p1 <- 1:p
        inv.Z <- chol2inv(lm.i$qr$qr[p1, p1, drop=FALSE])
        # p. 55 CC definition 
        if (GWR_args$se.fit.CCT) {
          C <- inv.Z %*% t(x) %*% diag(w.i)
          CC <- C %*% t(C)
          # only return coefficient covariance matrix diagonal raw values
          # for post-processing
          betase[i,] <- diag(CC)
        } else {
          betase[i,] <- diag(inv.Z)
        }
        # prediction "standard errors"
        # only return raw values for post-processing
        if (GWR_args$predictions && (p==m)) {
          if (GWR_args$se.fit.CCT) {
            pred.se[i] <- t(predx[i,]) %*% CC %*% predx[i,]
          } else {
            pred.se[i] <- t(predx[i,]) %*% inv.Z %*% predx[i,]
          }
        }
      }
    }
    # assigning residual bug Torleif Markussen Lunde 090529
    # cluster issue with fit.points 120505 Maximilian Spross
    if (!GWR_args$fp.given || GWR_args$fit_are_data) 
      gwr.e[i] <- ei[i]
    
    if (!GWR_args$fp.given && GWR_args$hatmatrix && (p==m)) 
      lhat[i,] <- t(x[i,]) %*% inv.Z %*% t(x) %*% diag(w.i)
  }
  df <- cbind(sum.w, betas, betase, gwr.e, pred, pred.se, localR2)
  if (!GWR_args$fp.given && GWR_args$hatmatrix) 
    return(list(df=df, lhat=lhat, bw=bw))
  else return(list(df=df, bw=bw))
} # GWR_int

gwr_miguel<-function (formula, data = list(), coords, bandwidth, gweight = gwr.Gauss,adapt = NULL, hatmatrix = FALSE, fit.points, longlat = NULL, se.fit = FALSE, weights, cl = NULL, predictions = FALSE, fittedGWRobject = NULL, se.fit.CCT = TRUE, distancias) {
  
  timings <- list()  #Registra el Tiempo
  .ptime_start <- proc.time()                                                     #Registra el Tiempo
  this.call <- match.call()                                                       #Registra el Tiempo
  p4s <- as.character(NA)                                                         #Vector en Blanco
  Polys <- NULL                                                                   #Poligonos Nulos
  
  if (is(data, "SpatialPolygonsDataFrame"))                                       #Creo Polígonos Espaciales
    Polys <- as(data, "SpatialPolygons")
  
  if (is(data, "Spatial")) {                                                      #Creo Mi Spatial Polygon Data
    if (!missing(coords)) 
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    p4s <- proj4string(data)
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
    data <- as(data, "data.frame")
  }
  
  if (is.null(longlat) || !is.logical(longlat))                                    #Control Log Lat
    longlat <- FALSE
  
  if (missing(coords))                                                             #Control Log Lat
    stop("Observation coordinates have to be given")
  
  if (is.null(colnames(coords)))                                                   #Archivo de Coordenadas
    colnames(coords) <- c("coord.x", "coord.y")
  
  mf <- match.call(expand.dots = FALSE)                                            #Objeto mf
  m <- match(c("formula", "data", "weights"), names(mf), 0)                        #
  mf <- mf[c(1, m)]                                                                #
  mf$drop.unused.levels <- TRUE                                                    #
  mf[[1]] <- as.name("model.frame")                                                #
  mf <- eval(mf, parent.frame())                                                   #
  mt <- attr(mf, "terms")                                                          # 
  dp.n <- length(model.extract(mf, "response"))                                    # 
  weights <- as.vector(model.extract(mf, "weights"))                               # 
  
  if (!is.null(weights) && !is.numeric(weights))                                   # Control sobre weights
    stop("'weights' must be a numeric vector")
  
  if (is.null(weights))                                                            # Si no hay Weights todos valen 1
    weights <- rep(as.numeric(1), dp.n)
  
  if (any(is.na(weights)))                                                         # Weights con NA
    stop("NAs in weights")
  
  if (any(weights < 0))                                                            # Weights con NEgativos
    stop("negative weights")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  lm <- lm.wfit(x, y, w = weights)
  lm$x <- x
  lm$y <- y
  gTSS <- c(cov.wt(matrix(y, ncol = 1), wt = weights, method = "ML")$cov *dp.n)
  
  if (hatmatrix) 
    se.fit <- TRUE
  if (hatmatrix) 
    predictions <- TRUE
  if (missing(fit.points)) {
    fp.given <- FALSE
    fittedGWRobject <- NULL
    predictions <- TRUE
    fit.points <- coords
    colnames(fit.points) <- colnames(coords)
    predx <- x
  }
  else fp.given <- TRUE
  
  griddedObj <- FALSE
  if (is(fit.points, "Spatial")) {
    if (predictions) {
      t1 <- try(slot(fit.points, "data"), silent = TRUE)
      if (inherits(t1, "try-error")) 
        stop("No data slot in fit.points")
      predx <- try(model.matrix(delete.response(mt), fit.points))
      if (inherits(predx, "try-error")) 
        stop("missing RHS variable in fit.points")
      if (ncol(predx) != ncol(x)) 
        stop("new data matrix columns mismatch")
    }
    Polys <- NULL
    if (is(fit.points, "SpatialPolygonsDataFrame")) {
      Polys <- as(fit.points, "SpatialPolygons")
      fit.points <- coordinates(fit.points)
    }
    else {
      griddedObj <- gridded(fit.points)
      fit.points <- coordinates(fit.points)
    }
  }
  else {
    if (predictions && fp.given) 
      stop("predictions not available for matrix fit points")
  }
  
  n <- NROW(fit.points)
  rownames(fit.points) <- NULL
  if (is.null(colnames(fit.points))) 
    colnames(fit.points) <- c("x", "y")
  if (predictions) {
    if (nrow(predx) != nrow(fit.points)) 
      stop("new data matrix rows mismatch")
    fit.points <- cbind(fit.points, predx)
  }
  fit_are_data <- isTRUE(all.equal(fit.points, coords, check.attributes = FALSE))
  input_predictions <- predictions
  if (fit_are_data && !predictions) {
    predictions <- TRUE
    input_predictions <- FALSE
    predx <- x
    fit.points <- cbind(fit.points, predx)
  }
  
  m <- NCOL(x)
  if (NROW(x) != NROW(coords)) 
    stop("Input data and coordinates have different dimensions")
  if (missing(bandwidth) && is.null(adapt)) 
    stop("Bandwidth must be given for non-adaptive weights")
  
  if (!is.null(adapt)) {
    stopifnot(is.numeric(adapt))
    stopifnot((adapt >= 0))
    stopifnot((adapt <= 1))
  }
  else {
    stopifnot(length(bandwidth) == 1)
  }
  
  if (missing(bandwidth)) 
    bandwidth <- NULL
  
  lhat <- NA
  yhat <- NULL
  
  if (!is.null(fittedGWRobject)) {
    yhat <- fittedGWRobject$SDF$pred
  }
  
  GWR_args <- list(fp.given = fp.given, hatmatrix = hatmatrix, 
                   longlat = longlat, bandwidth = bandwidth, adapt = adapt, 
                   se.fit = se.fit, predictions = predictions, se.fit.CCT = se.fit.CCT, 
                   fit_are_data = fit_are_data,distancias=distancias)
  
  timings[["set_up"]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  
  if (!is.null(cl) && length(cl) > 1 && fp.given && !hatmatrix) {                                       # Enmiezo con el GWR Optimization
    if (requireNamespace("parallel", quietly = TRUE)) {
      l_fp <- lapply(parallel::splitIndices(nrow(fit.points), 
                                            length(cl)), function(i) fit.points[i, , drop = FALSE])
      parallel::clusterEvalQ(cl, library(spgwr))
      varlist <- list("GWR_args", "coords", "gweight","y", "x", "weights", "yhat")
      env <- new.env()
      assign("GWR_args", GWR_args, envir = env)
      assign("coords", coords, envir = env)
      assign("gweight", gweight, envir = env)
      assign("y", y, envir = env)
      assign("x", x, envir = env)
      assign("weights", weights, envir = env)
      assign("yhat", yhat, envir = env)
      parallel::clusterExport(cl, varlist, env)
      res <- parallel::parLapply(cl, l_fp, function(fp) .GWR_int_miguel(fit.points = fp,coords = coords, gweight = gweight, y = y, x = x, 
                                                                        weights = weights, yhat = yhat, GWR_args = GWR_args))
      parallel::clusterEvalQ(cl, rm(varlist))
      rm(env)
      df <- list()
      df$df <- as.data.frame(do.call("rbind", lapply(res,function(x) x$df)))
      bw <- do.call("c", lapply(res, function(x) x$bw))
      results <- NULL
    }
    else {
      stop("parallel not available")
    }
  }
  else {
    df <- .GWR_int_miguel(fit.points = fit.points, coords = coords,gweight = gweight, y = y, x = x, weights = weights, 
                          yhat = yhat, GWR_args = GWR_args)
    if (!fp.given && hatmatrix) 
      lhat <- df$lhat
    bw <- df$bw
    results <- NULL
  }                                                                                                       # Fin de la GWR Optimization
  
  
  timings[["run_gwr"]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  if (predictions && !input_predictions) 
    predictions <- FALSE
  if (!fp.given && hatmatrix) {
    v1 <- sum(diag(lhat))
    B2 <- t(lhat) %*% lhat
    v2 <- sum(diag(B2))
    edf <- dp.n - 2 * v1 + v2
    B1 <- t(diag(dp.n) - lhat) %*% (diag(dp.n) - lhat)
    rss <- c(t(y) %*% B1 %*% y)
    delta1 <- sum(diag(B1))
    sigma2 <- rss/delta1
    odelta2 <- sum(diag(B1)^2)
    delta2 <- sum(diag(B1 %*% B1))
    nu1 <- sum(diag(B2))
    sigma2.b <- rss/dp.n
    AICb.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * pi) + (dp.n * ((n + v1)/(dp.n - 2 - v1)))
    AICh.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * pi) + dp.n + v1
    AICc.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * pi) + dp.n * ((delta1/delta2) * (dp.n + nu1))/((delta1^2/delta2) - 
                                                                                                               2)
    results <- list(v1 = v1, v2 = v2, delta1 = delta1, delta2 = delta2,sigma2 = sigma2, sigma2.b = sigma2.b, AICb = AICb.b, 
                    AICh = AICh.b, AICc = AICc.b, edf = edf, rss = rss, nu1 = nu1, odelta2 = odelta2, n = dp.n)
    timings[["postprocess_hatmatrix"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  
  if ((!fp.given || fit_are_data) && is.null(fittedGWRobject)) {
    localR2 <- numeric(n)
    if (is.null(adapt)) {
      
      bw <- bandwidth
      bandwidthR2 <- rep(bandwidth, n)
      
      bbbbbbbbbb <- distancias                                                                    # Ajusto el BandWith
      bandwidthR2<-distancias
    }
    else {
      
      bandwidthR2 <- gw_adapt_miguel(dp = coords, fp = fit.points[,1:2, drop = FALSE], quant = adapt, longlat = longlat)
      bw <- bandwidthR2
      
      bandwidthR2<-distancias
    }
    if (any(bandwidth < 0)) 
      stop("Invalid bandwidth")
    for (i in 1:n) {
      dxs <- spDistsN1(coords, fit.points[i, 1:2], longlat = GWR_args$longlat)
      if (any(!is.finite(dxs))) 
        dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
      w.i <- gweight(dxs^2, bandwidthR2[i])
      w.i <- w.i * weights
      if (any(w.i < 0 | is.na(w.i))) 
        stop(paste("Invalid weights for i:", i))
      RSS <- sum(w.i * (y - df$df[, "pred"])^2)
      yss <- sum(w.i * (y - weighted.mean(y, w.i))^2)
      localR2[i] <- 1 - (RSS/yss)
    }
    df$df <- cbind(df$df, localR2)
    timings[["postprocess_localR2"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  if (se.fit) {
    EDFS <- NULL
    normSigmaS <- NULL
    EDF <- NULL
    normSigma <- NULL
    if (fp.given && !is.null(fittedGWRobject)) {
      if (fittedGWRobject$hatmatrix) {
        EDF <- fittedGWRobject$results$edf
        normSigma <- sqrt(fittedGWRobject$results$rss/EDF)
        EDFS <- fittedGWRobject$results$n - fittedGWRobject$results$v1
        normSigmaS <- sqrt(fittedGWRobject$results$rss/EDFS)
      }
    }
    if (!fp.given && hatmatrix) {
      EDFS <- results$n - results$v1
      normSigmaS <- sqrt(results$rss/EDFS)
      EDF <- results$edf
      normSigma <- sqrt(results$rss/EDF)
    }
    ses <- grep("_se", colnames(df$df))
    senms <- colnames(df$df)[ses]
    betase <- df$df[, ses, drop = FALSE]
    df$df[, ses] <- NA
    if (predictions) {
      pred.se <- df$df[, "pred.se", drop = FALSE]
      df$df[, "pred.se"] <- NA
    }
    if (!is.null(EDF)) {
      betaseEDF <- normSigma * sqrt(betase)
      colnames(betaseEDF) <- paste(senms, "EDF", sep = "_")
      df$df[, ses] <- normSigmaS * sqrt(betase)
      df$df <- cbind(df$df, betaseEDF)
      if (predictions) {
        pred.se_EDF <- normSigma * sqrt(pred.se)
        df$df[, "pred.se"] <- normSigmaS * sqrt(pred.se)
        df$df <- cbind(df$df, pred.se_EDF)
      }
    }
    else {
      warning("standard errors set to NA, normalised RSS not available")
    }
    timings[["postprocess_SE"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  df <- as.data.frame(df$df)
  if (predictions) 
    fit.points <- fit.points[, 1:2, drop = FALSE]
  row.names(fit.points) <- row.names(df)
  SDF <- SpatialPointsDataFrame(coords = fit.points, data = df, 
                                proj4string = CRS(p4s))
  if (griddedObj) {
    gridded(SDF) <- TRUE
  }
  else {
    if (!is.null(Polys)) {
      df <- data.frame(SDF@data)
      rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i, 
                                                                       "ID"))
      SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
    }
  }
  timings[["final_postprocess"]] <- proc.time() - .ptime_start
  z <- list(SDF = SDF, lhat = lhat, lm = lm, results = results, 
            bandwidth = bw, adapt = adapt, hatmatrix = hatmatrix, 
            gweight = deparse(substitute(gweight)), gTSS = gTSS, 
            this.call = this.call, fp.given = fp.given, timings = do.call("rbind", 
                                                                          timings)[, c(1, 3)])
  class(z) <- "gwr"
  invisible(z)
}

get_5vecino <- function(row) {
  # Ordenar la fila de distancias y seleccionar la quinta más pequeña
  fifth_nearest_neighbor <- sort(row)[5]
  return(fifth_nearest_neighbor)
}



