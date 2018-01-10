## ------------------------------------------------------------------------
# load packages

# data manipulation
library(dplyr)
library(lubridate)

# forecasting
library(forecast)


## ------------------------------------------------------------------------

# Assume data is of the following form and it's collected hourly, 
# where n_uniq_tip = # of unique target IPs a source_ip communicates with
#       date_hour     source_ip n_uniq_tip
#   date_hour      source_ip   n_uniq_tip
#1  2016-06-01T00  10.0.1.2         37
#2  2016-06-01T00  10.0.1.3         30
#3  2016-06-01T00  10.0.1.4         43
#4  2016-06-01T00  10.0.1.5        134


## ------------------------------------------------------------------------
# Burn-in

fit_model <- function(x){
  myts <- ts(x, frequency = 24)
  ets(myts, model = "ANA", lower = c(0.01, 0, 0.0001, 0))
}


hw_model <- function(dt_burnin){
  # This function estimates the initial states and parameters for a weekend and weekday holt-winters model for each IP
  # It's necessary to have at least two seasonal cycles (24 hours) for both the weekend and weekday model
  # dt_burnin: data.frame, date_hour, source_ip, n_uniq_tip
  # Returns two data.frames:
  # m_weekday: data.frame, updated holt-winters parameters for weekday model
  # m_weekend: data.frame, updated holt-winters parameters for weekend model
    
  dt_burnin$weekend <- as.numeric(wday(substring(dt_burnin$date_hour, 1, 10)) %in% c(1, 7))
    
  dt_weekend <- filter(dt_burnin, weekend == 1)
  dt_weekday <- filter(dt_burnin, weekend != 1)
    
  # fit weekday model for each ip
  dt_weekday <- dt_weekday[order(dt_weekday$date_hour), ]
  m_weekday <- with(dt_weekday, by(n_uniq_tip, factor(source_ip), fit_model))
  
  par <- lapply(m_weekday, "[[", "par")
  nm_par <- names(par[[1]])
  par <- data.frame(matrix(unlist(par), byrow = TRUE, ncol = length(par[[1]])))
  colnames(par) <- nm_par

  initstate <- lapply(m_weekday, "[[", "initstate")
  nm_initstate <- names(initstate[[1]])
  initstate <- data.frame(matrix(unlist(initstate), byrow = TRUE, ncol = length(initstate[[1]])))
  colnames(initstate) <- nm_initstate
    
  res <- sapply(m_weekday, "[[", "residuals")
  res_mean <- apply(res, 2, mean)
  res_sd <- apply(res, 2, sd)

  m_weekday <- data.frame(source_ip = names(m_weekday), initstate, alpha = par$alpha, gamma = par$gamma, res_mean = res_mean, res_sd = res_sd)
    
  # fit weekend model for each ip
  dt_weekend <- dt_weekend[order(dt_weekend$date_hour), ]
  m_weekend <- with(dt_weekend, by(n_uniq_tip, factor(source_ip), fit_model))

  par <- lapply(m_weekend, "[[", "par")
  nm_par <- names(par[[1]])
  par <- data.frame(matrix(unlist(par), byrow = TRUE, ncol = length(par[[1]])))
  colnames(par) <- nm_par

  initstate <- lapply(m_weekend, "[[", "initstate")
  nm_inistate <- names(initstate[[1]])
  initstate <- data.frame(matrix(unlist(initstate), byrow = TRUE, ncol = length(initstate[[1]])))
  colnames(initstate) <- nm_initstate
  
  res <- sapply(m_weekend, "[[", "residuals")
  res_mean <- apply(res, 2, mean)
  res_sd <- apply(res, 2, sd)

  m_weekend <- data.frame(source_ip = names(m_weekend), initstate, alpha = par$alpha, gamma = par$gamma, res_mean = res_mean, res_sd = res_sd)
    
  return(list(m_weekday = m_weekday, m_weekend = m_weekend))
}

## ------------------------------------------------------------------------
# Hourly Update

hw_update <- function(dt_hr, m_weekday, m_weekend){
  # This compares the observed value to the predicted value of the time series and updates 
  # the parameters of the holt winters model.
  # dt_hr: data.frame, aggregated hourly data
  # m_weekday: data.frame, current parameter values for each holt-winters model 
  # m_weekend: data.frame, current parameter values for each holt-winters model
  #
  # Returns three data.frames:
  # results: data.frame, date_hour, source_ip, observed value, predicted value
  # m_weekday: data.frame, updated holt-winters parameters for weekday model
  # m_weekend: data.frame, updated holt-winters parameters for weekend model
    
  weekend_ind <- as.numeric(wday(substring(dt_hr$date_hour[1], 1, 10)) %in% c(1, 7))  
  s_ind <- paste0("s", as.numeric(substring(dt_hr$date_hour[1], 12, 13)) + 1) 
  r_ind <- paste0("r", as.numeric(substring(dt_hr$date_hour[1], 12, 13)) + 1)

  if(weekend_ind != 1){
    pred <- m_weekday$l + m_weekday[, s_ind]

    m <- match(m_weekday$source_ip, dt_hr$source_ip)
    obs <- dt_hr$n_uniq_tip[m]
      
    std_res <- ((obs - pred) - m_weekday$res_mean) / m_weekday$res_sd

    a_new <- m_weekday$alpha * (obs - m_weekday[, s_ind]) + (1 - m_weekday$alpha) * m_weekday$l
    s_new <- m_weekday$gamma * (obs - a_new) + (1 - m_weekday$gamma) * m_weekday[, s_ind]

    m_weekday$l <- a_new
    m_weekday[, s_ind] <- s_new
    results <- data.frame(date_hour = dt_hr$date_hour, source_ip = m_weekday$source_ip, obs = obs, pred = pred, std_res = std_res)
  } else{
    pred <- m_weekend$l + m_weekend[, s_ind]

    m <- match(m_weekend$source_ip, dt_hr$source_ip)
    obs <- dt_hr$n_uniq_tip[m]
      
    std_res <- ((obs - pred) - m_weekend$res_mean) / m_weekend$res_sd

    a_new <- m_weekend$alpha * (obs - m_weekend[, s_ind]) + (1 - m_weekend$alpha) * m_weekend$l
    s_new <- m_weekend$gamma * (obs - a_new) + (1 - m_weekend$gamma) * m_weekend[, s_ind]

    m_weekend$l <- a_new
    m_weekend[, s_ind] <- s_new
    results <- data.frame(date_hour = dt_hr$date_hour, source_ip = m_weekend$source_ip, obs = obs, pred = pred, std_res = std_res)
  }
  
  return(list(results = results, m_weekday = m_weekday, m_weekend = m_weekend))
}
    

      
    



## ------------------------------------------------------------------------
###########################
# Run code on sample data
###########################

## ------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(fasttime)

## ------------------------------------------------------------------------
# load data and split into data used for burn-in and data used to run the model in real time
ts_dt <- read.table("training-data", FALSE, "|", stringsAsFactors = FALSE)
colnames(ts_dt) <- c("date_hour","source_ip","n_uniq_tip")

ts <- fastPOSIXct(ts_dt$date_hour, "GMT")

# use first 7 days for training/burn-in
sel <- ts < (min(ts) + 86400*7)
dt_burnin <- ts_dt[sel,]

# use remaining days for "real-time" prediction
dt_rt <- ts_dt[!sel,]

## ------------------------------------------------------------------------
dt_burnin <- dt_burnin[order(dt_burnin$date_hour), ]
rownames(dt_burnin) <- NULL
dt_burnin[1:10,]

dt_rt <- dt_rt[order(dt_rt$date_hour), ]
rownames(dt_rt) <- NULL
dt_rt[1:10, ]


## ------------------------------------------------------------------------
# Initialize the model (training)
current_model <- hw_model(dt_burnin)

## ------------------------------------------------------------------------

## =====================
## ==== Acumos part ====
## =====================

library(acumos)

## create the model component
invisible(compose(predict=function(...,
                                   inputs=lapply(dt_rt, class), ## auto-detect input types from the training data
                                   outputs=c(date_hour="character", source_ip="character", obs="numeric", pred="numeric", std_res="numeric")) {
    dt_hr = data.frame(..., stringsAsFactors=FALSE) ## turn inputs into a data frame
    str(dt_hr) ## print the inputs for debugging
    ## update the model (on-line trianing) and predict
    current_model <<- hw_update(dt_hr, m_weekday = current_model$m_weekday, m_weekend = current_model$m_weekend)
    res <- current_model$results
    for (i in 1:2) res[[i]] <- as.character(res[[i]]) ## convert factors to characters
    print(res) ## print the result (for debugging)
}, aux=list(hw_update=hw_update, current_model=current_model))) ## pass the update function and trained model state


##
## to run, use
## acumos:::run(runtime=list(input_port=8101))
## or similar - either set output_url= to the next component or run
## acumos:::run(runtime=list(input_port=8101, data_response=TRUE))
## for REST-API mode (as opposed to default push-push mode)
##
## to push the component to the server use:
## push(<url>)
## if authentication is needed use
## push(<url>, token=auth(<auth-url>, <username>, <password>))
## add file=... to push() if pushing other .amc files
##
