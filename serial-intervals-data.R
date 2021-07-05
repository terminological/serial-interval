here::i_am("serial-intervals-data.R")

library(tidyverse)
source(here::here("common-setup.R"))
cache = PassthroughFilesystemCache$new(here::here("cache"))

options("ukcovid.reproduce.at"=as.Date("2021-06-29"))

## Serial interval estimation ----

serialIntervals = readxl::read_excel(here::here("input/serial-interval-literature.xlsx"))

metaParams = cache$getSaved("SERIAL-INTERVAL-META", params=list(serialIntervals), orElse = function(...) {
  tmp = serialIntervals %>% mutate(yi = mean_si_estimate, sei = (mean_si_estimate_high_ci-mean_si_estimate_low_ci)/3.92) %>% filter(!is.na(sei)) %>% filter(assumed_distribution == "gamma" & estimate_type %>% stringr::str_starts("serial"))
  meanfit = suppressWarnings(metaplus::metaplus(tmp$yi, tmp$sei, slab=tmp$label, random="mixture"))
      
  tmp2 = serialIntervals %>% mutate(yi = std_si_estimate, sei = (std_si_estimate_high_ci-std_si_estimate_low_ci)/3.92) %>% filter(!is.na(sei)) %>% filter(assumed_distribution == "gamma" & estimate_type %>% stringr::str_starts("serial"))
  sdfit = suppressWarnings(metaplus::metaplus(tmp2$yi, tmp2$sei, slab=tmp2$label, random="mixture"))
  
  serialIntervalMetaAnalysis = tribble(
    ~param, ~mean, ~sd, ~lower, ~upper,
    "mean", meanfit$results["muhat",1], NA, meanfit$results["muhat",2], meanfit$results["muhat",3],
    "sd", sdfit$results["muhat",1], NA, sdfit$results["muhat",2], sdfit$results["muhat",3]
  ) %>% mutate(dist="gamma")
  
  # save the data:
  save(serialIntervalMetaAnalysis, file = here::here("output/serialIntervalMetaAnalysis.rda"), compress = "bzip2", version = 2)
  # usethis::use_data(serialIntervalMetaAnalysis, overwrite = TRUE)
  return(serialIntervalMetaAnalysis)
})

siMeta = SerialIntervalProvider$metaAnalysis(dpc, metaDf = metaParams)


#### from ff100 

ff100 = dpc$spim$getFF100()
ff100dfit = cache$getSaved("SERIAL-INTERVAL-FF100",params=list(ff100), orElse = function(...) {
  
  dists=c("gamma","norm","weibull")
  tmp = ff100 %>% inner_join(ff100, by=c("ContactOf_FF100_ID"="FF100_ID"),suffix=c(".infectee",".infector"))
  tmp = tmp %>% select(FF100_ID,ContactOf_FF100_ID,contains("date"))
  dfit = DistributionFit$new(distributions = dists)
  
  if ("gamma" %in% dists) {
    dfit$models$gamma$start$shape = 1.1
    dfit$models$gamma$lower$shape = 1
    dfit$models$gamma$support = c(.Machine$double.xmin,Inf)
  }
  
  if ("weibull" %in% dists) {
    dfit$models$weibull$start$shape = 1.1
    dfit$models$weibull$lower$shape = 1
    dfit$models$weibull$support = c(.Machine$double.xmin,Inf)
  }
  
  if ("nbinom" %in% dists) {
    dfit$models$nbinom$support = c(.Machine$double.xmin,Inf)
  }
  
  if ("lnorm" %in% dists) {
    dfit$models$lnorm$support = c(.Machine$double.xmin,Inf)
  }
  
  if ("norm" %in% dists) {
    dfit$models$norm$start = NULL
  }
  
  tmp = tmp %>% mutate(
    left = as.integer(date_onset.infectee - date_onset.infector)-0.5,
    right = as.integer(date_onset.infectee - date_onset.infector)+0.5
  )
  
  dfit$fromCensoredData(tmp,lowerValueExpr = left,upperValueExpr = right,truncate = TRUE)
  return(dfit)
  
})
siFF100 = FittedSerialIntervalProvider$new(providerController = dpc,offset = 0,dfit = ff100dfit)

resamples = cache$getSaved("SERIAL-INTERVAL-RESAMPLE-DATA", orElse = function() {

  shortestCredibleSI = -7
  longestCredibleSI = 28
  samples=100
  boot.samples = NULL
  set.seed(101)
      # bootIterations = 250
      
  dfit = DistributionFit$new()
      
  parameterisedSIs = serialIntervals %>% 
    filter(
      estimate_type %>% stringr::str_starts("serial") &
      !assumed_distribution %in% c("empirical","unknown"),
    ) %>% 
    group_by(i = row_number()) %>% 
    group_modify(function(d,g,...) {
      paramDf = tribble(
        ~param, ~mean, ~sd, ~lower, ~upper,
        "mean", d$mean_si_estimate, (d$mean_si_estimate_high_ci - d$mean_si_estimate_low_ci)/3.96, d$mean_si_estimate_low_ci, d$mean_si_estimate_high_ci,
        "sd", d$std_si_estimate, NA, d$std_si_estimate_low_ci, d$std_si_estimate_high_ci
      )
      dfit$withSingleDistribution(dist = d$assumed_distribution, paramDf = paramDf, bootstraps = samples, N=d$sample_size, source = g$i)
      tibble()
      # side effect in group modify. This deserves 100 Hail Mary's
  })
      
  dfit$generateSamples(sampleExpr = N, seed=101)
  sampleDf = dfit$samples %>% select(bootstrapNumber, value, N, source)
      
  # include raw data for Xu et al, who did not produce parameterised estimates - only empirical distribution.
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7276042/bin/73153-2020.03.02.20029868-1.xlsx
  xuXls = dpc$datasets$download(id = "XU_SERIAL_INTERVAL",url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7276042/bin/73153-2020.03.02.20029868-1.xlsx", type = "xlsx")
  xuData = readxl::read_excel(xuXls)
  xuSamples = suppressWarnings(xuData %>% mutate(dt = as.integer(`p_Date of onset`)-as.integer(`s_Date of onset`)) %>% filter(!is.na(dt)) %>% pull(dt))
  xuSrc = max(sampleDf$source)+1
  
  # bootstrapping 90% of total into 100 samples
  N = length(xuSamples)*0.9
  set.seed(101)
  for(i in 1:samples) {
    # raw data also available from Du et al
    # https://github.com/MeyersLabUTexas/COVID-19/blob/master/Table%20S5%20medrxiv.xlsx?raw=true
    # their analysis is included though as resulted in normal distribution
    # N.B. it is probably the same data as Xu
    sampleDf = sampleDf %>% bind_rows(tibble(
      bootstrapNumber = i, 
      value = sample(xuSamples, size=N),
      N = N,
      source = xuSrc
    ))
  }
      
  # We can now: 
  # use this to directly estimate a parameterised distribution for the serial interval with or without shift
  # Or create a set of discretised empirical distributions for epiestim to explore.
  out = sampleDf %>% ungroup() #group_by(N,source)
  serialIntervalResampling = out %>% filter(value > shortestCredibleSI & value <= longestCredibleSI)
  save(serialIntervalResampling, file = here::here("output/serialIntervalResampling.rda"), compress = "bzip2", version = 2)
  # usethis::use_data(serialIntervalResampling, overwrite = TRUE)
  return(out)
})

## Resampled SI

siResample = cache$getSaved("SERIAL-INTERVAL-RESAMPLED",params=resamples, orElse = function(...) {return(NonParametricSerialIntervalProvider$new(dpc, samples=resamples))})
siResample$getBasicConfig(quick = FALSE)$si_sample %>% write.table(here::here("output/resampled-truncated-empirical-si-sample.txt"))

siResampleDfit = cache$getSaved("SERIAL-INTERVAL-RESAMPLED-FIT", params=list(resamples), orElse = function(...) {
  
  dists = c("norm","gamma","weibull")
  samples = resamples
  
  dfit = DistributionFit$new(distributions = dists, shifted = 0) #offset)
  
  if ("gamma" %in% dists) {
    dfit$models$gamma$start$shape = 1.1
    dfit$models$gamma$lower$shape = 1
    #dfit$models$gamma$support = c(0.01,Inf)
  }
  
  if ("weibull" %in% dists) {
    dfit$models$weibull$start$shape = 1.1
    dfit$models$weibull$lower$shape = 1
    #dfit$models$weibull$support = c(0.01,Inf)
  }
  
  if ("nbinom" %in% dists) {
    #dfit$models$nbinom$support = c(.Machine$double.xmin,Inf)
  }
  
  if ("lnorm" %in% dists) {
    dfit$models$lnorm$support = c(0.01,Inf)
  }
  
  if ("norm" %in% dists) {
    dfit$models$norm$start = NULL
  }
  
  #dfit$fromBootstrappedData(samples %>% ungroup() %>% select(bootstrapNumber, value))
  dfit$fromBootstrappedCensoredData(samples %>% ungroup() %>% mutate(left = value-0.5,right=value+0.5) %>% select(bootstrapNumber, left,right))
  
  return(dfit)
})





# defaultSI = SerialIntervalProvider$default(dpc)

## Incubation period estimation ----

bopUrl = "https://github.com/beoutbreakprepared/nCoV2019/raw/d8a68075fd17707fc8f861e93f775b6ec6e0158b/latest_data/latestdata.tar.gz"
latestdata2 = cache$getSaved("BOP-DATAFILE", params = list(bopUrl), orElse = function(...) {
  
  latestDataFile = dpc$datasets$downloadAndUntar(id = "BE-OUTBREAK-PREPARED",url = bopUrl, pattern = "latestdata.csv")
  
  latestdata <- readr::read_csv(
    latestDataFile,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  latestdata2 = latestdata %>% filter(!is.na(date_onset_symptoms) & !is.na(travel_history_dates))
  latestdata2 = latestdata2 %>% select(date_onset_symptoms,travel_history_dates) %>% 
    mutate(
      date_onset_symptoms = as.Date(date_onset_symptoms, "%d.%m.%Y"),
    ) %>% separate(
      travel_history_dates, into=c("travel.from","travel.to"),sep="[^0-9\\.]+", fill="left"
    ) %>%
    mutate(
      travel.from = as.Date(travel.from, "%d.%m.%Y"),
      travel.to = as.Date(travel.to, "%d.%m.%Y")
    ) %>%
    mutate(
      travel.from = ifelse(is.na(travel.from),travel.to,travel.from),
      left = as.integer(date_onset_symptoms - travel.to)-0.5,
      right = as.integer(date_onset_symptoms - travel.from)+0.5
    ) %>% filter(travel.to<date_onset_symptoms)
  return(latestdata2)
})

#devtools::load_all("~/Git/uk-covid-datatools/")
bopFit = cache$getSaved("BOP-FIT", params = list(latestdata2), orElse = function(...) {
  bopFit = DistributionFit$new(c("gamma","lnorm","weibull"))
  bopFit$models$weibull$lower$shape = 1
  bopFit$models$weibull$start$shape = 1.1
  bopFit$models$gamma$lower$shape = 1
  bopFit$models$gamma$start$shape = 1.1
  # bopFit$models$gamma$support = c(0.01,Inf)
  # bopFit$models$weibull$support = c(0.01,Inf)
  # bopFit$models$lnorm$support = c(0.01,Inf)
  bopFit$fromCensoredData(latestdata2,lowerValueExpr = left, upperValueExpr = right,truncate = TRUE,bootstraps = 200)
  return(bopFit)
})


# and from FF100

ff100 = dpc$spim$getFF100()

incubFF100Fit = cache$getSaved("INCUB-FIT", params=list(ff100), orElse = function(...) {
  
  censIncub = ff100 %>% filter(!is.na(date_exposure_first)) %>% select(date_exposure_first,date_exposure_last,date_onset) %>%
    mutate(
      right = as.numeric(date_onset-date_exposure_first)+0.5, 
      left = as.numeric(date_onset-date_exposure_last)-0.5) %>%
    mutate( 
      left = ifelse(left<=0,NA_integer_,left)
    ) %>% filter( right > 0)
  
  incubFF100Fit = DistributionFit$new(distributions = c("gamma","lnorm","weibull"))
  incubFF100Fit$models$weibull$lower$shape = 1
  incubFF100Fit$models$weibull$start$shape = 1.1
  incubFF100Fit$models$gamma$lower$shape = 1
  incubFF100Fit$models$gamma$start$shape = 1.1
  incubFF100Fit$fromCensoredData(censIncub,lowerValueExpr = left,upperValueExpr = right,truncate = TRUE, bootstraps = 200)
  return(incubFF100Fit)
  #incubFF100Fit$plot(xlim=c(0,7))
})




## Generation interval estimation ----

incubDist = cache$getSaved("INCUBATION-PERIOD-ERROR-SIMULATION", params=bopFit, orElse=function(...) {
  
  incubFit = bopFit$filterModels(aic == min(aic))
  incubFit$bootstraps = incubFit$bootstraps %>% filter(bootstrapNumber <= 100)
  
  # generate a set of samples from best fitting incubation period distribution
  # get samples from incubation and split into 2 groups - one for index patient and one for affected patient
  incubFit$generateSamples(sampleExpr = 2000,seed = 101)
  incubSamples = incubFit$samples %>%
    mutate(sampleCat = (sampleNumber-1) %/% 1000 + 1, sampleNumber = ((sampleNumber-1) %% 1000)+1) %>%
    pivot_wider(names_from = sampleCat, values_from = value, names_prefix = "incub")
  incubDist = incubSamples %>% mutate(delayOffset = incub2, incubError = incub1-incub2, transition = "infection to onset", from="infection", to = "onset")
  return(incubDist)
  
})

expt = cache$getSaved("GENERATION-INTERVAL-DATASET", params=list(incubDist,siResample), orElse = function(...) {
  
  actualSI = siResample$bootstrapSamples %>% select(bootstrapNumber,value)
  errorGItoSI = incubDist %>% filter(incub1 < 14 & incub2 < 14)  %>% select(bootstrapNumber, error=incubError) 
  expt = actualSI %>% group_by(bootstrapNumber) %>% nest(actual=value) %>% inner_join(
    errorGItoSI %>% group_by(bootstrapNumber) %>% nest(error=error),
    by="bootstrapNumber")
  return(expt)
  
})


## setup grid search and error miminisation

simSi = function(shape, rate, errors) {
  N = length(errors)
  genSim = rgamma(N*100,shape,rate = rate)+rep(errors,100)
  return(genSim)
}

errorFunction = function(simulated, actuals) {
  return(
    #abs(mean(actuals)-mean(simulated))+
    abs(IQR(actuals)-IQR(simulated))
  )
}

minimise = function(shape, rate, errors, actuals) {
  simulated = simSi(shape,rate,errors)
  out = errorFunction(simulated, actuals)
  return(out)
}

estimateParams = function(errors,actuals) {

  mu = mean(actuals)
  search = function(sdLim=c(mu/5,mu), grid = NULL) {
    #browser()
    sdWidth=(sdLim[2]-sdLim[1])/10
    if(sdWidth<0.00001) return(grid)
    for(sd in seq(sdLim[1]+sdWidth,sdLim[2]-sdWidth,length.out = 5)) {
      shape = mu^2/sd^2
      rate = mu/sd^2
      grid = grid %>% bind_rows(
        tibble(
          mean = mu,
          sd = sd,
          shape=shape,
          rate=rate,
          error=minimise(shape,rate,errors,actuals)
        )
      )
    }
    gridmin = grid %>% filter(error==min(error))
    grid = gridmin %>% group_modify(function(d,g,..) {
      return(search(
        d$sd+c(-1,1)*sdWidth,
        grid
      ))
    })
    return(grid %>% distinct())
  }

  return(search() %>% filter(error == min(error)))
}

## execute grid search

genIntTmp = dpc$getSaved("GENERATION-INTERVAL-OPTIM", params = expt, orElse = function(...) {
  genInt = expt %>% group_by(bootstrapNumber) %>% group_modify(function(d,g,...) {
    actuals=d$actual[[1]]$value
    errors = d$error[[1]]$error
    #initMean = mean(actuals)
    #initSd = sd(actuals)
    message(".",appendLF = FALSE)
    #out = optim(par = c(initMean,initSd), fn = minimise, lower = c(0,0), method = "L-BFGS-B", errors=errors, actuals=actuals)
    out = estimateParams(errors=errors, actuals=actuals)
    #return(tibble(mean = out$par[1],sd = out$par[2], err = out$value))
    return(out)
  })
  #genInt = genInt %>% mutate(shape = mean^2/sd^2, rate = mean/sd^2)
  generationIntervalSimulation = genInt %>% select(bootstrapNumber,shape,rate) %>% pivot_longer(cols = c(shape,rate), names_to = "param", values_to = "value") %>% mutate(dist="gamma")
  save(generationIntervalSimulation, file = here::here("output/generationIntervalSimulation.rda"), compress = "bzip2", version = 2)
  #usethis::use_data(generationIntervalSimulation, overwrite = TRUE)
  return(generationIntervalSimulation)
})

genIntFit = cache$getSaved("GENERATION-INTERVAL-FIT", params = genIntTmp, orElse = function(...) {
  genIntFit = DistributionFit$new()
  genIntFit$fromBootstrappedDistributions(fittedDistributions = genIntTmp)
  return(genIntFit)
})

siGeneration = cache$getSaved("GENERATION-INTERVAL", params = genIntFit, orElse = function(...) return(FittedSerialIntervalProvider$new(dpc, dfit = genIntFit)))

# export gen int sample probability distributions
tmp = siGeneration$getCustomConfigs(quick=FALSE)[[1]][[1]]$si_sample
tmp %>% write.table(here::here("output/generation-interval-fitted-si-sample.txt"))

# export distribution parameters
siGeneration$dfit$bootstraps %>% pivot_wider(names_from = "param",values_from = "value") %>% 
  mutate(mean = shape/rate, sd = sqrt(shape/(rate^2))) %>% write_csv(here::here("output/generation-interval-fitted-parameters.csv"))


## Impact of using the serial interval on estimation of *R~t~* ----



ukts = cache$getSaved("PHE-API", orElse = function(...) {
  ukts = dpc$datasets$getPHEApi(areaName = "england") %>% filter(date >= as.Date("2020-03-01") & date<=as.Date("2020-06-30") & type=="incidence") 
  
  ukts = ukts %>% bind_rows(
    ukts %>% group_by(statistic) %>% filter(date == min(date)) %>% select(-date) %>% mutate(value=0) %>% left_join(
      tibble(date = as.Date(as.Date("2020-02-16"):as.Date("2020-02-29"),"1970-01-01")), by=character()
    )
  )
  
  ukts = ukts %>% 
    tsp$imputeAndWeeklyAverage() %>% tsp$logIncidenceStats(smoothingWindow = 7)
  
  return(ukts)
})



siDates = tibble(
  `Start date` = as.Date(c("2020-03-19","2020-04-12","2020-05-12","2020-06-23")),
  `End date` = NA,
  `Label` = c("ascending","peak","early descending","late descending")
)

meanlim=c(2,8)
sdlim=c(1,6)

siImpact = cache$getSaved("SERIAL-INTERVAL-IMPACT", params = list(ukts,siDates,meanlim,sdlim), orElse = function(...) {
  # perform the Rt estimate for the various mean sd combinations
  
  ukcases = ukts %>% filter(statistic=="case")
  
  I = ukcases %>% pull(Est.value)
  dates = ukcases %>% pull(date)
  
  out = NULL
  
  for(endDate in siDates$`Start date`) {
  
    event = siDates$Label[siDates$`Start date` == endDate]
    start = match(endDate,dates)-3
    end = match(endDate,dates)+3
    
    for (siMean in seq(meanlim[1],meanlim[2],length.out = 40)) {
      #siMean = 4
      for (siSd in seq(sdlim[1],sdlim[2],length.out = 40)) {
      #siSd = 4  
        
        cfg_tmp = EpiEstim::make_config(mean_si = siMean, std_si = siSd, t_end = end, t_start = start, method="parametric_si")
        rEst = EpiEstim::estimate_R(I, method="parametric_si", config = cfg_tmp)
        out = bind_rows(
          out,
          tibble(
            startDate = as.Date(endDate-7,"1970-01-01"),
            endDate = as.Date(endDate,"1970-01-01"),
            label = event,
            siMean = siMean,
            siSd = siSd,
            medianR = rEst$R$`Median(R)`
          )                    
        )
      }
    }
    
  }
  return(out)
})



#serialIntervals = readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRdVV2wm6CcqqLAGymOLGrb8JXSe5muEOotE7Emq9GHUXJ1Fu2Euku9d2LhIIK5ZvrnGsinH11ejnUt/pub?gid=0&single=true&output=csv")

siFits = bind_rows(
  siResample$getSummary() %>% mutate(name = "resample"),
  siFF100$getSummary() %>% mutate(name = "ff100"),
  siMeta$getSummary()%>% mutate(name = "random effects"),
  siGeneration$getSummary() %>% mutate(name = "generation")
)

tmp = siImpact %>% left_join(siFits, by=character()) %>% group_by(name) %>% filter(abs(meanOfMean-siMean) == min(abs(meanOfMean-siMean)) & abs(meanOfSd-siSd) == min(abs(meanOfSd-siSd)))
tmp2 = tmp %>% ungroup() %>% group_by(startDate, label) %>% summarise(min = min(medianR), max= max(medianR)) %>% mutate(delta = max-min, percent = delta*2/(max+min)*100) %>%
  mutate(desc = sprintf("%s - %1.0f%% variation (*R~t~*: %1.2f to %1.2f)",label,percent,min,max))
summaryImpact = paste0(tmp2$desc,collapse = "; ")



## Time delays from infection to case identification, admission, and death ----

CHESS = dpc$spim$getCHESS()
CHESSClean = CHESS %>% dpc$chessProcessor()$chessAdmissionSubset(updatedWithin = 21)
#plots=list()
#tables=list()

totalCHESS = nrow(CHESSClean)
percentageCHESSWithDates = sprintf("%1.1f%%",length(na.omit(CHESSClean$estimateddateonset))/nrow(CHESSClean)*100)

rawOnset = cache$getSaved("CHESS-DATA", params = CHESSClean, orElse = function(...) {
  onsetToTest = CHESSClean %>% 
      filter(age>10 & !is.na(estimateddateonset) & !is.na(infectionswabdate)) %>% 
      mutate(
        transition = "onset to test",
        time = as.integer(infectionswabdate - estimateddateonset)+0.01
      ) %>% select(caseid,transition,time) %>% filter(time < 28 & time > -14) %>% group_by(transition)
  
  onsetToAdmission = CHESSClean %>% 
    # Onset ot admission
      filter(age>10 & !is.na(estimateddateonset) & !is.na(hospitaladmissiondate)) %>% 
      mutate(
        transition = "onset to admission",
        time = as.integer(hospitaladmissiondate - estimateddateonset)+0.01
      ) %>% select(caseid,transition,time) %>% filter(time < 100 & time > 0) %>% group_by(transition)
  
  onsetToDeath = CHESSClean %>% 
    filter(age>10 & !is.na(estimateddateonset) & !is.na(finaloutcomedate) & finaloutcome=="Death") %>% 
      mutate(
        transition = "onset to death",
        time = as.integer(finaloutcomedate - estimateddateonset)
      ) %>% select(caseid,transition,time) %>% filter(time < 100 & time > 0) %>% group_by(transition)
  
  rawOnset = bind_rows(
    onsetToTest, onsetToAdmission, onsetToDeath
  ) %>% group_by(transition)
  
  return(rawOnset)
})



# generate a observation-observation delay error function
# this is resampling raw CHESS delays to have a fixed number of bootstraps of given length & generating a delay distribution by looking at diffence between random resamples

infectionToObservationDelay = cache$getSaved("INFECTION-OBSERVATION-DELAY-DATA", params=list(rawOnset,incubDist), orElse=function(...) {
  rawDelay = rawOnset %>% mutate(transition = stringr::str_replace(transition,"onset to ","infection to ")) %>% group_by(transition) %>% group_modify(function(d,g,...) {
    boot = 1:100
    tmp = lapply(boot, function(b) tibble(
      bootstrapNumber = b,
      sampleNumber = 1:1000,
      delay1 = sample(d$time,size=1000),
      delay2 = sample(d$time,size=1000),
      delayOffset = delay2-delay1
    ))
    return(bind_rows(tmp))
  })
  
  # combine incubation period + observation delay of infectee period to get time correction
  combinedDelay = rawDelay %>% 
    left_join(incubDist, by=c("bootstrapNumber","sampleNumber"), suffix=c("",".incub")) 
  
  combinedDelay %>% 
    mutate(timeDelay = incub2+delay2) %>% ungroup() %>%
    select(transition,timeDelay,bootstrapNumber,sampleNumber) %>% return()
})

infectionToObservationFit = cache$getSaved("INFECTION-OBSERVATION-DELAY-FIT", params=infectionToObservationDelay, orElse = function(...) {
  fit = DistributionFit$new(distributions = c("lnorm","gamma","weibull"))
  fit$fromBootstrappedData(infectionToObservationDelay %>% select(-sampleNumber) %>% group_by(transition), valueExpr = timeDelay)
  return(fit)
})





timeCorrection2 = infectionToObservationDelay %>%
  mutate(to = stringr::str_remove(transition,"infection to ")) %>%
  group_by(to) %>%
  summarise(
    mean = mean(timeDelay),
    sd = sd(timeDelay),
    lower = quantile(timeDelay,0.025),
    upper = quantile(timeDelay,0.975),
  ) %>% 
  bind_rows(
    incubDist %>% group_by(to) %>% summarise(
      mean = mean(delayOffset),
      sd = sd(delayOffset),
      lower = quantile(delayOffset,0.025),
      upper = quantile(delayOffset,0.975),
    )
  ) 

timeCorrection2 %>% write.csv(file=here::here("output/TIME-CORRECTION.csv"))

ukCovidObservationDelays = timeCorrection2 %>%
  mutate(
    `Mean delay (days)` = sprintf("%1.2f",mean),
    `SD (days)` = sprintf("%1.2f",sd),
    `95% quantiles (days)` = sprintf("%1.2f; %1.2f",lower, upper)
  )
save(ukCovidObservationDelays, file = here::here("output/ukCovidObservationDelays.rda"), compress = "bzip2", version = 2)



## Estimation of *R~t~*  using formal and pragmatic methods ----


formalRt = cache$getSaved("FORMAL-RT", params = list(siGeneration,infectionToObservationFit,ukts), orElse = function(...) {

  convolutionsFit = infectionToObservationFit$clone()
  convolutionsFit = convolutionsFit$filterModels(aic == min(aic))
  pmfs = convolutionsFit$discreteProbabilities(q = 0:30)
  pmfs = pmfs %>% mutate(statistic = case_when(
    transition == "infection to test" ~ "case",
    transition == "infection to admission" ~ "hospital admission",
    transition == "infection to death" ~ "death",
    TRUE ~ NA_character_
  ))
  
  ukts2 = ukts %>% filter(type=="incidence") %>% group_by(statistic) %>% group_modify(function(d,g,...) {
  
    caseSts = surveillance::sts(
      observed = d %>% pull(Imputed.value) ,
      epoch = d %>% pull(date)
    )
    
    casePmf = pmfs %>% filter(statistic == g$statistic) %>% pull(Mean.discreteProbability)
    
    # Back-propagate to get an estimate of the original, assuming our observations
    # are a Poisson process
    
    bpnp.control <- list(k=0,eps=rep(0.005,2),iter.max=rep(250,2),B=-1) #,verbose=TRUE)
    bpnp.control2 <- modifyList(bpnp.control, list(hookFun=NULL,k=2,B=100,eq3a.method="C"))
    
    #Fast C version (use argument: eq3a.method="C")! 
    sts.bp <- surveillance::backprojNP(caseSts, incu.pmf=casePmf, control=bpnp.control2)
    
    d = d %>% mutate(
      Est.value = sapply(1:dim(sts.bp@lambda)[1], FUN=function(x) {mean(sts.bp@lambda[x,,])}),
      Est.Quantile.0.975.value = sapply(1:dim(sts.bp@lambda)[1], FUN=function(x) {quantile(sts.bp@lambda[x,,],0.975)}),
      Est.Quantile.0.025.value = sapply(1:dim(sts.bp@lambda)[1], FUN=function(x) {quantile(sts.bp@lambda[x,,],0.025)}),
      Est.Quantile.0.75.value = sapply(1:dim(sts.bp@lambda)[1], FUN=function(x) {quantile(sts.bp@lambda[x,,],0.75)}),
      Est.Quantile.0.25.value = sapply(1:dim(sts.bp@lambda)[1], FUN=function(x) {quantile(sts.bp@lambda[x,,],0.25)})
    )
    
    return(d)
  
  })
  
  # incidence = (ukts2 %>% tsp$plotIncidenceQuantiles(colour = statistic, dates = c("2020-03-01","2020-07-01"))) + 
  #   geom_line(aes(y=Deconv.value,colour=statistic),linetype="21") + 
  #   geom_ribbon(aes(ymin=Deconv.0.025.value,ymax=Deconv.0.975.value,group=statistic),fill="black",colour=NA,alpha=0.2) + 
  #   guides(fill="none",colour="none")+
  #   facet_wrap(vars(statistic))
  
  
  ukts3 = ukts2 %>%
    mutate(value=Est.value) %>%
    tsp$estimateRt(valueVar = value, serialIntervalProvider = siGeneration, window = 14)
  #ggplot(result,aes(x=date,fill=statistic,colour=statistic))+geom_point(aes(y=observed))+geom_line(aes(y=smoothed))+geom_line(aes(y=Deconv.value),linetype="dashed")+geom_ribbon(aes(ymin=lower,ymax=upper),fill="blue",colour=NA,alpha=0.2)
  return (ukts3)
  
})


pragmaticRt = cache$getSaved("PRAGMATIC-RT", params = list(timeCorrection2, siResample, ukts), orElse = function(...) {
  
  offsetAssumptions = timeCorrection2 %>% inner_join(
    tribble( ~to, ~name,
             "admission", "hospital admission",
             "admission", "icu admission",
             "death", "death",
             "test", "case",
             "onset", "symptom",
             "onset", "triage"
             ), by="to") %>% 
    select(name,mean) %>% deframe()
  
  ukts4 = ukts %>%
    tsp$logIncidenceStats(valueVar = value, growthRateWindow = 14, smoothingWindow = 7) %>%
    tsp$estimateRt(valueVar = value, serialIntervalProvider = siResample, window = 14) %>%
    tsp$adjustRtDates(window = 0, offsetAssumptions = offsetAssumptions)
  return(ukts4)
})

uktsComparison = bind_rows(
  formalRt %>% mutate(subgroup = "formal"),
  pragmaticRt %>% mutate(subgroup = "pragmatic")
)



