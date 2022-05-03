library(gfoRmula)
library(optmatch)
library(survey)
library(svyVGAM)
library(dplyr)
library(data.table)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
####################                               ######################
####################   survey panel data           ######################
####################                               ######################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

######### example data analysis


{
  outcome_type <- 'binary_eof'
  id <- 'id_num'
  time_name <- 'time'
  covnames <- c('cov1', 'cov2', 'treat')
  outcome_name <- 'outcome'
  histories <- c(lagged, cumavg)
  histvars <- list(c('treat', 'cov1', 'cov2'), c('cov1', 'cov2')) # time varying covariate or confounder
  covtypes <- c('binary', 'zero-inflated normal', 'normal')
  covparams <- list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 +
                                    cov3 + time,
                                  cov2 ~ lag1_treat + cov1 + lag1_cov1 +
                                    lag1_cov2 + cov3 + time,
                                  treat ~ lag1_treat + cumavg_cov1 +
                                    cumavg_cov2 + cov3 + time))
  ymodel <- outcome ~ treat + cov1 + cov2 + lag1_cov1 + lag1_cov2 + cov3
  intvars <- list('treat', 'treat')
  interventions <- list(list(c(static, rep(0, 7))),
                        list(c(static, rep(1, 7))))
  # list(c(threshold, 1, Inf)))
  int_descript <- c('Never treat', 'Always treat')
  nsimul <- 10000
  ncores <- 2
}

gform_bin_eof <- gformula(obs_data = binary_eofdata,
                          outcome_type = outcome_type, id = id,
                          time_name = time_name, covnames = covnames,
                          outcome_name = outcome_name, covtypes = covtypes,
                          covparams = covparams, ymodel = ymodel,
                          intvars = intvars, interventions = interventions,
                          int_descript = int_descript, histories = histories,
                          histvars = histvars, basecovs = c("cov3"),
                          seed = 1234, parallel = TRUE, nsamples = 5,
                          nsimul = nsimul, ncores = ncores)
gform_bin_eof



#####################################
#####################################
#####################################
#####################################
#####################################
############ dat_anl_2 ##############
#####################################
#####################################
#####################################
#####################################
#####################################


# chronic~binge+age_cat+sex+income_cat+edu+stress : total

## data preprocessing 
{
  dat_anl_2 <- dat_anl_2 %>% filter(year>2013)
  dd <- dat_anl_2 %>% group_by(PIDWON) %>% summarise(a=max(diag_year)>2013) %>% 
    filter(a==TRUE)
  dd <- dd$PIDWON
  temp <- (dat_anl_2 %>% group_by(PIDWON) %>% summarise(nn=if(all(chronic==2)) 0 else 1))
  nope <- temp[temp$nn==0,]$PIDWON
  dat_anl_2 <- dat_anl_2 %>% filter(PIDWON %in% c(dd,nope))
  dat_anl_2$chronic_p <- ifelse(dat_anl_2$chronic %in% c(2,4),0,1)
  
  # dat_anl_2 <- dat_anl_2 %>% filter(!is.na(I_WGL_08))
  # dd <- (dat_anl_2 %>% filter(!(I_WGL_08 %>% is.na)))$PIDWON %>% unique
  # dat_anl_2 <- dat_anl_2 %>% filter(PIDWON %in% dd)
  
  ## time 1~4
  df <- dat_anl_2 %>% group_by(PIDWON) %>% summarise(year,time=year-min(year))
  dat_anl_2 <- dat_anl_2 %>% left_join(df,by=c("PIDWON","year"))
  dat_anl_2 <- dat_anl_2 %>% arrange(PIDWON,year)
  min_time <- min((dat_anl_2 %>% as.data.table)[["time"]])
  correct_time_indicator<- tapply((dat_anl_2 %>% as.data.table)[["time"]], (dat_anl_2 %>% as.data.table)[["PIDWON"]],
                                  FUN = function(x){
                                    all(x == min(min_time, 0):(length(x)+min(min_time, 0)-1))
                                  })
  dat <- correct_time_indicator %>% as.data.frame
  dat[,2] <- rownames(dat)
  dat <- data.frame(PIDWON=dat[,2],b=dat[,1])
  dd <- dat[dat$b,1]
  dat_anl_2 <- dat_anl_2 %>% filter(PIDWON %in% dd)
  lastid <- (dat_anl_2%>% group_by(PIDWON) %>% summarise(n=n(),t=n>3) %>% filter(t==T))$PIDWON
  dat_anl_2 <- dat_anl_2 %>% filter(PIDWON %in% lastid)
  
  ## stress na impute........
  for(i in lastid){
    dat_anl_2[dat_anl_2$PIDWON==i,"stress"][1:2] <- 
      dat_anl_2[dat_anl_2$PIDWON==i,"stress"][3:4] %>% mean %>% round
  }
  dat_anl_2$stress <- dat_anl_2$stress %>% as.factor

}

# compare the proportion of each covariate between case and control group
((dat_anl_2 %>% group_by(PIDWON) %>% summarise(smoke_new,binge)) %>% unique)[,2:3] %>% 
  table(useNA = "ifany")
(dat_anl_2 %>% group_by(PIDWON) %>% summarise(smoke_new,binge) %>% unique)[,2:3] %>% table

## total population
{
  outcome_type <- 'binary_eof'
  id <- 'PIDWON'
  time_name <- 'time'
  covnames <- c('chronic_p','binge','stress','who_walk',"income","edu","bmi","smoke_new")
  covtypes <- c('binary','binary',"categorical","binary","normal","normal","normal","binary")
  outcome_name <- 'chronic_p'
  histories <- c(lagged,cumavg)
  histvars <- list(c('chronic_p','binge','stress','who_walk','bmi','smoke_new'),c("bmi","income","edu")) # time-varying confounder
  covparams <- list(covmodels = c(chronic_p~1,
                                  binge ~ lag1_chronic_p+lag1_binge+smoke_new+stress+lag1_stress+edu+income+sex+max_age+time,
                                  stress~lag1_chronic_p+lag1_stress+who_walk+income+cumavg_bmi+sex+max_age+time,
                                  who_walk~lag1_chronic_p+lag1_who_walk+max_age+cumavg_bmi+income+edu+sex+time,
                                  income~cumavg_income+cumavg_edu+sex+max_age+time,
                                  edu~cumavg_edu+sex+max_age,
                                  bmi~cumavg_bmi+stress+binge+who_walk+lag1_who_walk+income+sex+max_age+time,
                                  smoke_new~lag1_smoke_new+binge+stress+income+edu+sex+max_age))#,
  ymodel <- chronic_p ~ lag1_chronic_p+binge+lag1_binge+smoke_new+lag1_smoke_new+stress+lag1_stress+who_walk+lag1_who_walk+income+edu+bmi+sex+max_age
  intvars <- list('binge', 'binge')
  interventions <- list(list(c(static, rep(0, 4))),
                        list(c(static, rep(1, 4))))
  int_descript <- c('not binge',"so binge")
  nsimul <- 20000
  ncores <- 3
}

## binge drinking vs no drinking(binge 4 vs 0)
{
  outcome_type <- 'binary_eof'
  id <- 'PIDWON'
  time_name <- 'time'
  covnames <- c('chronic_p','binge','stress','who_walk',"income","edu","bmi","smoke_new")
  covtypes <- c('binary','binary',"categorical","binary","normal","normal","normal","binary")
  outcome_name <- 'chronic_p'
  histories <- c(lagged,cumavg)
  histvars <- list(c('chronic_p','binge','stress','who_walk','bmi',"smoke_new"),c("bmi","income","edu")) # time-varying confounder
  covparams <- list(covmodels = c(chronic_p~1,
                                  binge ~ lag1_chronic_p+smoke_new+stress+lag1_stress+edu+income+sex+max_age+time,
                                  stress~lag1_chronic_p+lag1_stress+who_walk+income+cumavg_bmi+sex+max_age+time,
                                  who_walk~lag1_chronic_p+lag1_who_walk+max_age+cumavg_bmi+income+edu+sex+time,
                                  income~cumavg_income+cumavg_edu+sex+max_age+time,
                                  edu~cumavg_edu+sex+max_age,
                                  bmi~cumavg_bmi+stress+binge+who_walk+lag1_who_walk+income+sex+max_age+time,
                                  smoke_new~lag1_smoke_new+binge+stress+income+edu+sex+max_age))
  ymodel <- chronic_p ~ lag1_chronic_p+binge+stress+lag1_stress+who_walk+lag1_who_walk+income+edu+bmi+sex+max_age
  intvars <- list('binge', 'binge')
  interventions <- list(list(c(static, rep(0, 4))),
                        list(c(static, rep(1, 4))))
  int_descript <- c('not binge',"so binge")
  nsimul <- 20000
  ncores <- 3
}

gform_bin_eof_temp_2 <- gformula(obs_data = dat_anl_2 %>% as.data.table,
                                   outcome_type = outcome_type, id = id,
                                   time_name = time_name, covnames = covnames,
                                   outcome_name = outcome_name, covtypes = covtypes,
                                   covparams = covparams, ymodel = ymodel,
                                   intvars = intvars, interventions = interventions,
                                   int_descript = int_descript, histories = histories,
                                   histvars = histvars, basecovs = c("sex","max_age"),
                                   seed = 1234, parallel = TRUE, nsamples = 100,
                                   nsimul = nsimul, ncores = ncores)

gform_bin_eof_temp_2$result

### causal estimand and plot

xx <- gform_bin_eof_temp_2$result[2:3,4:5]
del_me <- function(x) { return(1/( x*(1-x) )) }
odds_se <- rep(0,2)
odds_rt <- c(lodds(xx[2,1]),lodds(xx[1,1])) %>% as.numeric
odds_se <- c(xx[2,2]*del_me(xx[2,1]), xx[1,2]*del_me(xx[2,1])) %>% as.numeric

plot(x=seq(-3,3,by=0.1),y=seq(1:3))

a <- c(
  xx[1,1]-qnorm(0.975)*xx[1,2]/10,
  xx[1,1]+qnorm(0.975)*xx[1,2]/10,
  xx[2,1]-qnorm(0.975)*xx[2,2]/10,
  xx[2,1]+qnorm(0.975)*xx[2,2]/10,
  xx[2,1]-xx[1,1] - qnorm(0.975)*sqrt(xx[1,2]^2+xx[2,2]^2)/10,
  xx[2,1]-xx[1,1] + qnorm(0.975)*sqrt(xx[1,2]^2+xx[2,2]^2)/10,
odds_rt[1] - qnorm(0.975)*odds_se[1]/10,
odds_rt[1] + qnorm(0.975)*odds_se[1]/10,
odds_rt[2] - qnorm(0.975)*odds_se[1]/10,
odds_rt[2] + qnorm(0.975)*odds_se[1]/10,
(odds_rt[1]-odds_rt[2]) - qnorm(0.975)*sqrt(odds_se^2 %>% sum),
(odds_rt[1]-odds_rt[2]) + qnorm(0.975)*sqrt(odds_se^2 %>% sum)
)

plot(x=c(a[1:4]),y=rep(1,each=4), type="p",
     main="Mean 95% CI", xlab="Mean", ylab="", xlim=c(0,0.1))
lines(x=a[1:2],y=rep(1,2), lwd=2)
lines(x=a[3:4],y=rep(1,2),lwd=2)

plot(x=c(a[5:6]),y=rep(1,each=2),main="Mean Difference 95% CI",
     type="l",lwd=2, xlim=c(0,0.03), xlab="Diff")
abline(v=0,col="red",lty=3,lwd=2)

plot(x=c(a[7:10]),y=rep(1,each=4),xlim=c(-3.05,-2.2), type="p",
     main="Log odds 95% CI", xlab="Log Odds", ylab="")
lines(x=a[7:8],y=rep(1,2), lwd=2)
lines(x=a[9:10],y=rep(1,2),lwd=2)

plot(x=c(a[11:12]),y=rep(1,each=2),xlim=c(-0.0005,0.5),main="Log OR 95% CI",
     type="l",lwd=2, xlab="Log OR")
abline(v=0,col="red",lty=3,lwd=2)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
####################### only hypertension patients who participated on survey from 2014 to 2017
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
{
  dd <- dat_all %>% filter(year==2014) %>% filter(CD2 %in% c(1,3))
  dd <- dd$PIDWON
  patient <- dat_all %>% filter(PIDWON %in% dd)
  patient <- patient[,-3]
  colnames(patient)[1:23] <- 
    c("PIDWON","year","chronic","doctor","diag_year","medicine","smoke","walk","walk_time",
      "stress1","stress2","stress3","sleep_week","sleep_weekend","height","kg","avg_alc","bool_alc",
      "sex","age","marry","edu","income")
  patient <- patient %>% filter(year>2013)
}

## patient categorization
{
  dd <- patient %>% group_by(PIDWON) %>% summarise(chronic_new=if(all(chronic==2)) 0 else 1)
  patient <- patient %>% left_join(dd,by=c("PIDWON"))
  patient$medicine_new <- ifelse(patient$medicine==1,2,ifelse(patient$medicine==2,1,0)) 
  patient$smoke_new <- ifelse(patient$smoke==4,0,1) # no smoking history:0, have smoking history:1
  patient$who_walk <- ifelse(patient$walk==3,ifelse(patient$walk_time %in% 5:6,1,0),
                             ifelse(patient$walk==4,ifelse(patient$walk_time %in% 4:6,1,0),
                                    ifelse(patient$walk==5,ifelse(patient$walk_time %in% 3:6,1,0),
                                           ifelse(patient$walk==6,ifelse(patient$walk_time %in% 2:6,1,0),0) 
                                    )
                             )
  ) # walking more than 150 minutes per week : 1
  dd <- patient %>% group_by(PIDWON) %>% summarise(max_age=max(year)-unique(age))
  patient <- patient %>% left_join(dd,by="PIDWON")
  # patient$max_age %>% quantile()
  patient$age_cat <- ifelse(patient$max_age<30,0,ifelse(patient$max_age<65,1,2)) # AGE < 30: 0, 30~64 : 1, 65 >= : 2
  
  patient$marry_new <- ifelse(patient$marry==1,1,0) # marry and live together : 1,  o.w. :0
  
  # stress
  dd <- patient %>% group_by(PIDWON,year) %>% summarise(stress=if(all(is.na(stress1),is.na(stress2),is.na(stress3))) NA 
                                                        else{
                                                          if(!any(is.na(stress1),is.na(stress2),is.na(stress3))){
                                                            if(sum(c(stress1,stress2,stress3) %in% 1)==0) 0
                                                            else if(sum(c(stress1,stress2,stress3) %in% 1)==1) 1
                                                            else if(sum(c(stress1,stress2,stress3) %in% 1)==2) 2
                                                            else if(sum(c(stress1,stress2,stress3) %in% 1)==3) 3
                                                          }
                                                        })
  patient <- patient %>% left_join(dd, by=c("PIDWON","year"))
  patient <- patient %>% arrange(PIDWON,year)
  patient$edu_cat <- ifelse(patient$edu<24,0,ifelse(patient$edu<46,1,2)) 
  patient$bmi_cat <- ifelse(patient$bmi<25,0,ifelse(patient$bmi<30,1,2))
  patient$sex <- ifelse(patient$sex==1,1,0)
  patient$chronic_p <- ifelse(patient$chronic %in% c(2,4),0,1)
  lastid <- (patient%>% group_by(PIDWON) %>% summarise(n=n(),t=n>3) %>% filter(t==T))$PIDWON
  patient <- patient %>% filter(PIDWON %in% lastid)

  ## stress na impute........
  for(i in lastid){
    patient[patient$PIDWON==i,"stress"][1:2] <- 
      patient[patient$PIDWON==i,"stress"][3:4] %>% mean %>% round
  }
  
  df <- patient %>% group_by(PIDWON) %>% summarise(year,time=year-min(year))
  patient <- patient %>% left_join(df,by=c("PIDWON","year"))
  patient <- patient %>% arrange(PIDWON,year)
  patient$stress <-  as.factor(patient$stress)
  patient$medicine_new <- as.factor(patient$medicine_new)
}

### g-formula modeling for hypertension patients
{
  outcome_type <- 'binary_eof'
  id <- 'PIDWON'
  time_name <- 'time'
  covnames <- c('binge','stress','medicine_new',"income","bmi","smoke_new")
  covtypes <- c('binary',"categorical",'categorical',"normal","normal","binary")
  outcome_name <- 'chronic_p'
  histories <- c(lagged,cumavg)
  histvars <- list(c('binge','stress','medicine_new','smoke_new'),c("bmi","income")) # 시간에 따라 변하는? 요인?
  covparams <- list(covmodels = c(binge ~ lag1_binge+smoke_new+lag1_medicine_new+stress+income+sex+max_age+time,
                                  stress~lag1_stress+income+cumavg_bmi+max_age+time,
                                  medicine_new~lag1_medicine_new+lag1_binge+time,
                                  # income~cumavg_income+cumavg_edu+time,
                                  # edu~cumavg_edu,
                                  # bmi~cumavg_bmi+stress+income+max_age+time,
                                  income~cumavg_income+sex+max_age+time,
                                  # edu~cumavg_edu+sex+max_age,
                                  bmi~cumavg_bmi+stress+binge+income+sex+max_age+time,
                                  smoke_new~lag1_smoke_new+binge+stress+medicine_new+sex+max_age+time))
  ymodel <- chronic_p ~ binge + lag1_binge+smoke_new+lag1_smoke_new+lag1_medicine_new+stress+income+bmi+sex+max_age
  intvars <- list('binge', 'binge')
  interventions <- list(list(c(static, rep(0, 4))),
                        list(c(static, rep(1, 4))))
  int_descript <- c('not binge',"so binge")
  nsimul <- 2000
  ncores <- 3
}

gform_bin_eof_temp_2 <- gformula(obs_data = dat_anl_2 %>% as.data.table,
                                   outcome_type = outcome_type, id = id,
                                   time_name = time_name, covnames = covnames,
                                   outcome_name = outcome_name, covtypes = covtypes,
                                   covparams = covparams, ymodel = ymodel,
                                   intvars = intvars, interventions = interventions,
                                   int_descript = int_descript, histories = histories,
                                   histvars = histvars, basecovs = c("sex","max_age"),
                                   seed = 1234, parallel = TRUE, nsamples = 100,
                                   nsimul = nsimul, ncores = ncores)
