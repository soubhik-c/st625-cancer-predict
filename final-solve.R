#install.packages("librarian", repos = "http://cran.us.r-project.org")

# libs --------------------------------------------------------------------
librarian::shelf(rstudioapi,
                 renv,
                 mltools,
                 ggplot2,
                 ggcorrplot,
                 MASS,
                 dplyr,
                 data.table,
                 caret,
                 stringr,
                 olsrr,
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

regionmap[which(regionmap$state == tolower("Washington")),2]

map_states <- function(s) {
  s=tolower(str_trim(s))
  regionmap[which(regionmap$state == s),2]
}

imputed_df['region']=sapply(imputed_df$state, map_states)

# Check for normality ----------------
plot(Hmisc::describe(imputed_df)) 
drp_cols=-c(9,13,32:34)
data.table(Column=names(imputed_df[,drp_cols]),
           P=apply(imputed_df[,drp_cols], 2,
           function(x) shapiro.test(x)$p.value)
           )

# correlation matrix ---------
fsel=imputed_df[,drp_cols]
fsel<-cbind(fsel, one_hot(as.data.table(as.factor(imputed_df$region))))

fsel
corMtrx <- cor(fsel)
print(corMtrx)
highcorr=findCorrelation(corMtrx, cutoff=0.5)
print(highcorr)
#fsel[,highcorr]

corr = round(cor(fsel), 2)
ggcorrplot(corr,title="Correlation heatmap among cancer variables")

# Attempt3-------
subset.model=regsubsets(TARGET_deathRate~.,data=fsel, nbest=20,method=c("exhaustive"))
summary(subset.model)

plot(subset.model, scale="r2")
plot(subset.model, scale="adjr2")
plot(subset.model, scale="bic")
plot(subset.model, scale="Cp")

features=c("incidenceRate","povertyPercent","PctHS18_24","PctBachDeg25_Over",
  "PctPrivateCoverage","PctEmpPrivCoverage","PctOtherRace","PctMarriedHouseholds",
  "V1_Mid West", "V1_South",
  "TARGET_deathRate")

ln=length(features)
seq_along(features[-ln])

models<-list()
for(i in seq_along(features[-ln])){
  models[[i]]<- lm(TARGET_deathRate~.,data=fsel[,features[c(1:i,ln)]])
}
lapply(models,summary)
lapply(models,anova)

