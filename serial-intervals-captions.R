
.def = list(
  fig = captioner::captioner(prefix="Figure"),
  tab = captioner::captioner(prefix="Table"),
  sfig = captioner::captioner(prefix="Supplemental figure"),
  stab = captioner::captioner(prefix="Supplemental table")
)

.index = list(
  fig = function(name,caption) .def$fig(name,caption,display=FALSE),
  tab = function(name,caption) .def$tab(name,caption,display=FALSE),
  sfig = function(name,caption) .def$sfig(name,caption,display=FALSE),
  stab = function(name,caption) .def$stab(name,caption,display=FALSE)
)

cap = list(
  fig = function(name) .def$fig(name,display="full"),
  tab = function(name) .def$tab(name,display="full"),
  sfig = function(name) .def$sfig(name,display="full"),
  stab = function(name) .def$stab(name,display="full")
)

ref = list(
  fig = function(name) .def$fig(name,display="cite"),
  tab = function(name) .def$tab(name,display="cite"),
  sfig = function(name) .def$sfig(name,display="cite"),
  stab = function(name) .def$stab(name,display="cite")
)

.index$fig("infection-timeline","A timeline of events associated with a single infector-infectee pair in a transmission chain")
.index$tab("rt-methods","Comparison of 2 approaches for estimating the reproduction number and parameters needed to support each approach")
.index$tab("lit-review","Sources of serial interval estimates from a literature search")
.index$fig("si-estimates","Panel A - Days between infected infectee disease onset based on a resampling of published estimates from the literature and Panel B - Estimates of serial interval from FF100 data. The histogram in panel A shows the combined density of all sets of samples within the original research.")
.index$fig("incub-period","Incubation period distributions reconstructed from Open COVID-19 Data Working Group and from FF100 data. Histogram data is approximate due to interval censoring.")
.index$tab("incub-period-gof","Goodness of fit statistics for incubation period distributions reconstructed from Open COVID-19 Data Working Group and from FF100 data")
.index$fig("generation-interval", "Estimated generation interval distributions, from resampled serial intervals as predictor, and estimated serial intervals from incubation period combined with samples from a generation interval assumed as a gamma distributed quantity.")
.index$fig("epi-curve","Epidemic curve for cases, deaths and hospital admissions are used for analysis in this paper. Dashed vertical lines show dates at which we conduct our analysis, chosen to represent the ascending, peak, early and late descending phases of cases during the first wave in the UK.")
.index$fig("si-impact","Time varying reproduction number given various assumptions on the serial interval mean and standard deviation. The blue points show the central estimate of serial intervals from the literature, whereas the coloured error bars show the mean and standard deviation of the 2 serial interval (green, violet) and 1 generation interval (orange) estimates presented in this paper. Contours show the *R~t~* estimate for that combination of mean and standard deviation serial interval. The four panels represent the 4 different time points investigated.")
.index$fig("chess-delay","Panel A: Time delay distributions from symptom onset to test (diagnosis or case identification), admission or death, estimated from CHESS data set, plus in Panel B estimated delays from infection to observation, and can be negative in certain cases. based on incubation period and observation delay. These  can be used for deconvolution.")
.index$tab("rt-offset","Estimated time delays between infection and various observations over the course of an infection, based on the combination of incubation period and symptom onset to observation delay")
.index$fig("deconv-compare","Panel A: The epidemic incidence curves in England for different observations (orange - formal) and inferred estimates of infection rates (green - pragmatic) based on a deconvolution of the time delay distributions. Panel B, the resulting *R~t~* values calculated either using infection rate estimates and generation interval (formal subgroup) or unadjusted incidence of observation, and serial interval (pragmatic). The *R~t~* estimated direct from observed incidence curves (pragmatic) have their dates adjusted by the mean delay estimate")

.index$stab("resampling-algo","An algorithm for resampling representative serial intervals based on paramterised distributions and raw data from published studies")
.index$stab("chess-trusts","The NHS trusts in the CHESS dataset who report all admissions")
.index$sfig("meta-analysis","Forest plot for serial interval studies for the mean of the serial interval, using the normal mixture random effect model, and from studies identified in the literature which give confidence intervals")
.index$stab("dfit-serial-interval-resample","Parametererised serial interval distributions from resampling the literature. Gamma and weibull estimates are from data truncated at zero. AIC estimates are not comparable to those for Normal distribution which is fitting all data, including negative serial intervals, and hence has a lower mean.")
.index$stab("dfit-serial-interval-ff100","Parametererised serial interval distributions from FF100. Gamma and weibull estimates are from data truncated at zero. AIC estimates are not comparable to those for Normal distribution which is fitting all data, including negative serial intervals, and hence has a lower mean.")
.index$stab("incub-period-detail","Distribution details for estimated incubation period distributions reconstructed from Open COVID-19 Data Working Group and from FF100 data")
.index$sfig("gof-incub-fit","Graphical goodness of fit for parameterised incubation period distributions fitted to the Be Outbreak Prepared dataset")
.index$sfig("parameter-distributions","The distribution of the parameters of the fitted generation interval estimates")
.index$stab("chess-delay-params","Time delay distributions estimated from CHESS data set, for both transitions from disease onset to case, admission or death, and presumed infection and case, admission or death")

