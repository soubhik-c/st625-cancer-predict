#install.packages("librarian", repos = "http://cran.us.r-project.org")

# libs --------------------------------------------------------------------
librarian::shelf(rstudioapi,
                 renv,
                 mltools,
                 GGally,
                 ggplot2,
                 ggcorrplot,
                 MASS,
                 dplyr,
                 data.table,
                 caret,
                 stringr,
                 olsrr,
                 leaps,
                 doMC,
                 gbm,
                 CombMSC,
                 cran_repo = 'https://cran.r-project.org')

#renv::clean()

# init --------------------------------------------------------------------
set.seed(12673)
setwd(dirname(getSourceEditorContext()$path))
cancer_reg <- read.csv("cancer_reg.csv")

timeit <- function(body) {
  print(paste("Running", sep=""))
  start_time <- Sys.time() 
  res <- body()
  end_time <- Sys.time()
  print(paste("time taken is", format(end_time - start_time)))
  return(res)
}

# data cleanups -----------------------------------------------------------

cancer_reg$binnedInc<-as.factor(cancer_reg$binnedInc)
cancer_reg[c('county', 'state')] <- str_split_fixed(cancer_reg$Geography, ',', 2)
paste("States ", length(unique(cancer_reg$state)), 
      " County", length(unique(cancer_reg$county)))
# cancer_reg$state<-as.factor(cancer_reg$state)
# cancer_reg$county<-as.factor(cancer_reg$county)

# md.pattern(cancer_reg)
VIM::aggr(cancer_reg)
# Loose NA columns
imputed_df<-subset(cancer_reg, select = -c(PctEmployed16_Over,PctPrivateCoverageAlone,PctSomeCol18_24) )

regionmap <- list()
regionmap[['North East']] = c(
  'Connecticut', 
  'Maine',
  'Massachusetts', 
  'New Hampshire', 
  'Rhode Island',
  'Vermont',
  'New Jersey',
  'New York',
  'Pennsylvania'
)
regionmap[['Mid West']] = c(
  'Illinois', 
  'Indiana',
  'Michigan', 
  'Ohio',
  'Wisconsin',
  'Iowa', 
  'Kansas',
  'Minnesota',
  'Missouri',
  'Nebraska', 
  'North Dakota', 
  'South Dakota'
)

regionmap[['South']] = c(
  'Delaware',
  'Florida',
  'Georgia',
  'Maryland',
  'North Carolina',
  'South Carolina',
  'Virginia',
  'District of Columbia',
  'West Virginia',
  'Alabama',
  'Kentucky', 
  'Mississippi',
  'Tennessee',
  'Arkansas',
  'Louisiana',
  'Oklahoma',
  'Texas'
)

regionmap[['West']] = c(
  'Arizona',
  'Colorado',
  'Idaho',
  'Montana',
  'Nevada',
  'New Mexico',
  'Utah',
  'Wyoming',
  'Alaska',
  'California',
  'Hawaii',
  'Oregon',
  'Washington'
)

regionmap=reshape2::melt(regionmap, value.name = c("state"))
names(regionmap) = c("state", "region")
regionmap$state=tolower(regionmap$state)

map_states <- function(s) {
  s=tolower(str_trim(s))
  regionmap[which(regionmap$state == s),2]
}

imputed_df['region']=sapply(imputed_df$state, map_states)

# Check for normality ----------------
plot(Hmisc::describe(imputed_df)) 
names(imputed_df)
drp_cols=-c(9,13,32:34)
data.table(Column=names(imputed_df[,drp_cols]),
           P=apply(imputed_df[,drp_cols], 2,
           function(x) shapiro.test(x)$p.value)
           )

# correlation matrix ---------
fsel=imputed_df[,drp_cols]

fsel<-cbind(fsel,
            one_hot(
              as.data.table(
                lapply(imputed_df[,c("region", "binnedInc")], as.factor)
                )
              ))

fsel
corMtrx <- cor(fsel)
print(corMtrx)
highcorr=findCorrelation(corMtrx, cutoff=0.5)
print(highcorr)
# interim[,highcorr]

corr = round(cor(fsel), 2)
ggcorrplot(corr,title="Correlation heatmap among cancer variables")


# # Split data -------------
# selected<-createDataPartition(interim$TARGET_deathRate,p=.90,list=F)
# fsel<-interim[selected]
# ftest<-interim[-selected]

# Attempt3-------
subset.model=regsubsets(TARGET_deathRate~.,data=fsel,nbest=20,method=c("exhaustive"))
summary(subset.model)
# 
plot(subset.model, scale="r2")
plot(subset.model, scale="adjr2")
plot(subset.model, scale="bic")
plot(subset.model, scale="Cp")

#stepAIC
step.model <- stepAIC(lm(TARGET_deathRate ~ ., data=fsel), direction = "both", 
                      trace = F)
summary(step.model)

#olsrr
# reg=lm(TARGET_deathRate ~ ., data=fsel)
# summary(reg)
# bestsubset=ols_step_best_subset(reg)
# bestsubset
# allpos=ols_step_all_possible(reg)
# allpos
# plot(allpos)

#RandomForests
# registerDoMC(cores=8)
trmod <- train(TARGET_deathRate ~ ., 
                  data=fsel,
                  method="gbm",
                  trControl=trainControl(
                    method="repeatedcv",
                    number=20,
                    repeats = 3,
                    allowParallel=T
                  ),
                )
fselImp <- varImp(trmod)
print(fselImp)

features=c(
  "incidenceRate",
  "povertyPercent",
  "PctHS18_24",
  "PctBachDeg25_Over",
  "PctPrivateCoverage",
  "PctEmpPrivCoverage",
  "PctOtherRace",
  "PctMarriedHouseholds",
  "region_Mid West",
  "region_South",
  "TARGET_deathRate")

imp_feat<-c("PctPublicCoverageAlone",
            "povertyPercent",
            "AvgHouseholdSize",
            "PctUnemployed16_Over",
            "PctBlack",
            "PctHS18_24",
            "TARGET_deathRate")

# length(intersect(names(fsel), features))
# fsel[,features]

eda_feat<-function(apply_feat) {
  print(GGally::ggpairs(fsel[,apply_feat]))
  summary(fsel[,apply_feat])
}

do_model<-function(apply_feat) {
  ln=length(apply_feat)
  seq_along(apply_feat[-ln])
  
  models<-list()
  for(i in seq_along(apply_feat[-ln])){
    models[[i]]<- lm(TARGET_deathRate~.,data=fsel[,apply_feat[c(1:i,ln)]])
  }
  return(
    list(
      models,
      lapply(models,summary),
      do.call(anova,models),
      lapply(models, PRESS.lm)
    )
  )
}


# eda_feat(features)
leapft<-do_model(features)

gbm_feat<-c(
  "incidenceRate",
  "PctBachDeg25_Over",
  # "avgDeathsPerYear",
  # "medIncome",
  # "PctHS25_Over",
  # "popEst2015",
  # "povertyPercent",
  "PctPublicCoverageAlone",
  # "avgAnnCount",
  # "PctPrivateCoverage",
  "region_West",
  "region_South",
  "PctUnemployed16_Over",
  "PctOtherRace",
  "PctHS18_24",
  # "MedianAgeFemale",
  # "AvgHouseholdSize",
  # "PctMarriedHouseholds",
  "PercentMarried",
  "PctBlack",
  "TARGET_deathRate"
)
do_model(gbm_feat)


gbm_feat2<-c(
  "incidenceRate",
  "PctBachDeg25_Over",
  # "avgDeathsPerYear",
  # "medIncome",
  # "PctHS25_Over",
  # "popEst2015",
  "povertyPercent",
  "PctPublicCoverageAlone",
  "avgAnnCount",
  # "PctPrivateCoverage",
  "region_West",
  "region_South",
  "PctUnemployed16_Over",
  # "PctOtherRace",
  # "PctHS18_24",
  # "MedianAgeFemale",
  # "AvgHouseholdSize",
  # "PctMarriedHouseholds",
  # "PercentMarried",
  # "PctBlack",
  "TARGET_deathRate"
)

extractmodel<-function(mft, modnum) {
  return(list(mft[[1]][[modnum]],
              mft[[4]][[modnum]]
  ))
}

calcr2jack<-function(pressval) {
  SSyy=(nrow(fsel)-1)*var(fsel$TARGET_deathRate)
  return (1-(pressval/SSyy))
}

# eda_feat(gbm_feat)
glmft2<-do_model(gbm_feat2)
glmft2

par(mfrow=c(2,2))

x<-extractmodel(leapft, 10)
y<-extractmodel(glmft2, 8)

x.selmod=x[[1]]
x.pressval=x[[2]]

y.selmod=y[[1]]
y.pressval=y[[2]]

plot(x.selmod)
plot(y.selmod)

x.pressval
y.pressval

calcr2jack(x.pressval)
calcr2jack(y.pressval)

par(mfrow=c(1,1))

results<-predict(y.selmod,fsel)
data.frame(results,actual=fsel$TARGET_deathRate) %>%
  ggplot(aes(x=results,y=actual)) + 
  geom_point()+
  stat_smooth(method="lm",show.legend = T)


hist(fsel$TARGET_deathRate)

