# functions to be used in etas estimation/simulation

dmagn_mark <- function(x, data, params) {
  if (params[7] > 0) {
    lambda <- PtProcess::etas_gif	(data, x[, "time"], params = params[1:5])
    y <- dgamma(x[, "magnitude"], shape = 1 + sqrt(lambda) * params[7],
                rate = params[6], log = TRUE)
  } else y <- dexp(x[, "magnitude"], rate = params[6], log = TRUE)
  return(y)
}
rmagn_mark <- function(ti, data, params) {
  if (params[7]>0) {
    lambda <- PtProcess::etas_gif	(data, ti, params = params[1:5])
    y <- rgamma(1, shape = 1 + sqrt(lambda) * params[7], rate = params[6])
  } else y <- rexp(1, rate = params[6])
  return(list(magnitude = y))
}
expmap <- function(y, p) {
  y$params[1:5] <- exp(p)
  return(y)
}

# eq_catalog.R
# generate an earthquake catalog in the format that can be used for the estimation
# of Hawkes processes
##################################################################################

EQcatalog <- function(data, t_start = 1923, t_end = 2015, depth = 100, origin = "1970-01-01",
                      mag_min = 4, lat_range = c(15, 55), lng_range = c(118,161)) {
  # given raw sample, this function selects the sample with specified sample period,
  # magnitude threshold, depth, and coordinate ranges
  # choose only data with known hypocenter and are natural (Subsidiary information ="1")
  # and maximum intensity is known
  sample.data = data[(data$Precision!="N" | data$Date>="2015-01-01") &
                       data$Max_Intensity !=" " &
                       (data$Sub_info==1 | is.na(data$Sub_info)) &
                       data$Depth<=depth, ]
  # remove the all NA rows
  sample.data = sample.data[complete.cases(sample.data$Rec_Id), ]

  # choose sample period
  sample.start = as.Date(paste0(t_start, "-01-01"))
  sample.end = as.Date(paste0(t_end, "-12-31"))
  sample.data = sample.data[sample.data$Date >= sample.start &
                              sample.data$Date<=sample.end, ]

  # choose region
  sample.data = sample.data[sample.data$Latitude  <= lat_range[2] &
                              sample.data$Latitude >=lat_range[1] &
                              sample.data$Longitude <= lng_range[2] &
                              sample.data$Longitude >= lng_range[1], ]
  #  select magnitude>=thres, neglect all records where magnitude is unknown
  sample.data = sample.data[sample.data$Magnitude_1 >= mag_min & (!is.na(sample.data$Magnitude_1)),  ]
  n = nrow(sample.data)

  catalog = data.frame(date = character(n))
  catalog$date = as.Date(sample.data$Date)
  catalog$hour = sample.data$Hour
  catalog$minute = sample.data$Minute
  catalog$second = sample.data$Second
  catalog$latitude = sample.data$Latitude
  catalog$longitude = sample.data$Longitude
  catalog$magnitude = sample.data$Magnitude_1
  catalog$time = julian(catalog$date, origin = as.Date(origin)) +
    ((catalog$second/60 + catalog$minute)/60 + catalog$hour)/24
  catalog$depth = sample.data$Depth
  catalog$max_intensity = sample.data$Max_Intensity
  if(sum(duplicated(catalog$time))!=0){
    print("duplicated records")
    catalog = catalog[-which(duplicated(catalog$time)), ]
  }

  return(catalog)
}


# etas_estim.R
# this function estimates the parameters of the etas model
# input: initial parameters, data, convergence criterion
# output: mpp object x0 with fitted parameter
etas_estim <- function(catalog, city = NULL, city_range = NULL, magMin, params, t0, tN){
  # transform magnitude and time
  catalog$magnitude <- catalog$magnitude -  magMin
  catalog <- catalog[catalog$magnitude>=0, ]
  # time should be already in numeric form, with origin specified; origin is set to be start of the sample.
  # select estimation period and initial parameter for the average magnitude
  TT <- c(0, julian(as.Date(tN), origin = as.Date(t0)) )
  params[6] <- 1/mean(catalog$magnitude)
  # choose the range of coordinates of given city
  if(length(city)!=0){
    cityR <- city_range[city_range$City==city, ]
    # select only the events inside given city
    catalog <- catalog[catalog$latitude <= cityR$latMax &
                         catalog$latitude >= cityR$latMin &
                         catalog$longitude <= cityR$lngMax &
                         catalog$longitude >= cityR$lngMin, ]
  }
  # estimate parameters for given catalog
  x <- PtProcess::mpp(data = catalog, gif = PtProcess::etas_gif	,
           mark = list(dmagn_mark, rmagn_mark), params = params, TT = TT,
           gmap = expression(params[1:5]), mmap = expression(params))
  initial <- log(params[1:5])
  z <- optim(initial, neglogLik, object = x, pmap = expmap,
             control = list(trace = 0, maxit = 100))
  initial <- z$par
  z <- nlm(neglogLik, initial, object = x, pmap = expmap,
           print.level = 0, iterlim = 500, typsize = initial)
  x0 <- expmap(x, z$estimate)
  return(x0)
}

#etas_test.R
# this function performs KS-type tests on residuals of the ETAS model
# input: mpp object x with estimated parameters
# output: test stats and (if required) plots
etas_test <- function(object, plot = F){
  x0 <- object
  N <- nrow(x0$data)
  alp1 = 0.05
  alp2 = 0.01
  K_alp_1 = sqrt(-(log(alp1/2)/2))
  K_alp_2 = sqrt(-(log(alp2/2)/2))
   #kstat1 <- ks.test(residuals(x0)/N, 'punif', alternative = 'two.sided')
  yk = diff(residuals(x0))
  uk = 1-exp(-yk)
  kstat <- ks.test(uk, 'punif', alternative = 'two.sided')

  if(plot==T){
    par(mfrow=c(2,2))
    plot(residuals(x0),1:N, type='l',
         ylab = "Event Number", xlab = "Transformed Time",
         col='red')
    # title(paste("Year", tStart, "-", tEnd, ", Magnitude", magMin,
    #             "+, Latitude ", latMin,"-",latMax,", Longitude ", lngMin,"-",lngMax))
    lines(0:0.1:N, 0:0.1:N, type='l')
    lines(0:0.1:N, 0:0.1:N + K_alp_2*sqrt(N), type='l')
    lines(0:0.1:N, 0:0.1:N - K_alp_2*sqrt(N), type='l')
    lines(0:0.1:N, 0:0.1:N + K_alp_1*sqrt(N), type='l')
    lines(0:0.1:N, 0:0.1:N - K_alp_1*sqrt(N), type='l')
    #
    plot(ecdf(uk),do.points=F, xlim=c(0,1),ylim=c(0,1),
         ylab = "Cumulative Distribution", xlab = "U_(k)",
         col='red',main=NULL)
    lines(0:0.01:1, 0:0.01:1)
    lines(0:0.01:1, 0:0.01:1 + K_alp_2/sqrt(N-1))
    lines(0:0.01:1, 0:0.01:1 - K_alp_2/sqrt(N-1))
    lines(0:0.01:1, 0:0.01:1 + K_alp_1/sqrt(N-1))
    lines(0:0.01:1, 0:0.01:1 - K_alp_1/sqrt(N-1))

  }
  return(kstat)
}


# etas_prob.R
# This function calculates the eq probability forecasts
# within a given time period through simulation
# input: time period, number of simulation, parameters, data
# output: earthquake probability for threshold 1, threshold 2, ...
etas_prob <- function(object, n.sim = 10000,
                    time.length = 90, TT = NULL,
                    threshold = 0, mag.Min = 4.5,
                    aggregate = 1){
  # "object" is an mpp object (marked point process, see package
  # "PtProcess") containing event history and parameters
  # set the simulation threshold
  if(length(TT)==0){
    TT <- c(object$TT[2], object$TT[2] + time.length)
    object$TT <- TT
  } else if(length(TT)==2 & TT[2] > TT[1]){
    object$TT <- TT
  } else
    stop("undefined time interval TT")
  gif <- object$gif
  # simulate n times the point process in the interval TT
  n.rec <- data.frame(matrix(NA, nrow = n.sim, ncol = length(threshold)))
  colnames(n.rec) <- paste0("magnitude", mag.Min + threshold, "+")
  for(i in 1:n.sim){
    #set.seed(i)
    sim = simulate(object)$data
    # retrieve the simulated data
    sim_data = sim[which(sim$time<TT[2] & sim$time>=TT[1]), ]
    # compare simulated data with each threshold and store the number in n.rec
    compare = vapply(sim_data$magnitude, function(x) x>=threshold,
                     numeric(length(threshold)))
    if(length(threshold)>1){
      n.rec[i, ] = (rowSums(compare) !=0)
    } else {
      n.rec[i, ] = (sum(compare) !=0)
    }

  }
  # simulated probability is number of events exceeding threshold divided by n.sim
  prob.threshold <- data.frame(matrix(NA, nrow = 1, ncol = length(threshold)))
  colnames(prob.threshold) <- paste0("magnitude", mag.Min + threshold, "+")
  prob.threshold <- colSums(n.rec) / n.sim
  if(aggregate == 1){
    return(as.numeric(prob.threshold))
  } else{
    return(n.rec)
  }
}


###################################################################################
# gen_Xpsi_city.R
# parallelizable function wrapper to generate eq probabilities
gen_Xpsi_city <- function(Iter.val, list.Est, city_range, time, time.end,
                          date.start, threshold, magMin,
                          n.sim = 1e3, time.out=60){
  XpsiCity <- matrix(0, nrow=length(time), ncol=length(city_range$City))
  colnames(XpsiCity) <- city_range$City
  rownames(XpsiCity) <- time
  XpsiCity <- as.data.frame(XpsiCity)
  l.XpsiCity <- rep(list(XpsiCity), length(threshold))
  print(paste("# simulations:", n.sim))
  if(is.na(time.out)){
    for(i in 1:length(city_range$City)){
      city <- city_range$City[i]
      obj <- list.Est[[i]]
      for(j in 1:length(time)){
        TT1 <- time[j]
        TT2 <- ifelse(j<(length(time)), time[j+1], time.end)
        # simulation. start of simulation period is TT1(beginning of this quarter) and end is start of next quarter
        probs <- etas_prob(obj, n.sim = n.sim,
                         TT = c(julian(as.Date(TT1), origin = date.start),
                                julian(as.Date(TT2), origin = date.start)),
                         threshold = threshold - magMin[i],
                         mag.Min = magMin[i],
                         aggregate = 1)
        for(l in 1:length(threshold)){
          l.XpsiCity[[l]][j, i] <- probs[l]
        }
      }
    }
  }
  else {
    for(i in 1:length(city_range$City)){
      city <- city_range$City[i]
      obj <- list.Est[[i]]
      for(j in 1:length(time)){
        # print(time[j])
        TT1 <- time[j]
        TT2 <- ifelse(j<(length(time)), time[j+1], time.end)
        
        try1 <- tryCatch(withTimeout(etas_prob(obj, n.sim = n.sim,
                                               TT = c(julian(as.Date(TT1), origin = date.start),
                                                 julian(as.Date(TT2), origin = date.start)),
                                             threshold = threshold - magMin[i],
                                             mag.Min = magMin[i],
                                             aggregate = 1),
                                     timeout = time.out, elapsed = time.out, onTimeout="silent"),
                         error=function(e) {})
        while(length(try1) == 0){
          try1 <- tryCatch(withTimeout(etas_prob(obj, n.sim = n.sim,
                                                 TT = c(julian(as.Date(TT1), origin = date.start),
                                                        julian(as.Date(TT2), origin = date.start)),
                                                 threshold = threshold - magMin[i],
                                                 mag.Min = magMin[i],
                                                 aggregate = 1),
                                       timeout = time.out, elapsed = time.out, onTimeout="silent"),
                           error=function(e) {})
          
        }
          for(l in 1:length(threshold)){
            l.XpsiCity[[l]][j, i] <- try1[l]
          }
        
      }
    }
  }
  #print(paste("simulation", Iter.val, "finished"))
  names(l.XpsiCity) <- paste0("threshold_", threshold)
  return(l.XpsiCity)
}
