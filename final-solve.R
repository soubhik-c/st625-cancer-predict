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
                 cowplot,
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

hist(imputed_df$TARGET_deathRate)

# correlation matrix ---------
cor.df=imputed_df[,drp_cols]

corMtrx <- cor(cor.df)

colnames(corMtrx)
highcorr=findCorrelation(corMtrx, names=T, exact=F, cutoff=0.7)
names(cor.df[,highcorr])

dims=length(colnames(corMtrx))
corrDF <- expand.grid(row = 1:dims, col = 1:dims)
corrDF$correlation <- as.vector(corMtrx)
levelplot(correlation ~ row + col, corrDF)

corr = round(cor(cor.df), 2)
ggcorrplot(corr,title="Correlation heatmap among cancer variables")


# Feature Selection-------
fsel.df <- as.data.frame(unclass(imputed_df), stringsAsFactors = TRUE)
names(fsel.df)
fsel<-fsel.df[,-c(13,32,33)]
names(fsel)

subset.model=regsubsets(TARGET_deathRate~., really.big=T, data=fsel,nbest=20,method=c("exhaustive"))
summary(subset.model)

plot(subset.model, scale="r2")
plot(subset.model, scale="adjr2")
plot(subset.model, scale="bic")
plot(subset.model, scale="Cp")

features=c(
  "incidenceRate",
  "PctHS18_24",
  "PctBachDeg25_Over",
  "PctPublicCoverageAlone",
  "PctOtherRace",
  "PctMarriedHouseholds",
  "region",
  "TARGET_deathRate")

# EDA--------------
eda.df<-fsel[,features]

# rank distribution
ggplot(data=eda.df,aes(x=TARGET_deathRate,fill=PctOtherRace))+
  geom_bar(position='dodge')+
  facet_grid(. ~ region)+
  ggtitle("Bar-plot showing rank distribution, separated by gender and discipline")

#histogram
ggplot(data=eda.df,aes(x=TARGET_deathRate,fill=PctOtherRace))+
  geom_density(alpha=0.4)+
  facet_grid(. ~ region)+
  ggtitle("Density curves showing salary distribution, separated by gender and discipline")

# boxplots
ggplot(data=eda.df,aes(x=PctOtherRace,y=TARGET_deathRate,fill=PctOtherRace))+
  geom_boxplot(notch=TRUE)+
  facet_grid(. ~ region)+
  ggtitle("Notched boxplots showing deathRate distribution, separated by region and discipline")

# violinplots
ggplot(data=eda.df,aes(x=PctOtherRace,y=TARGET_deathRate,fill=PctOtherRace))+
  geom_violin()+
  geom_boxplot(fill='darkred',width=0.1,notch=TRUE)+
  geom_point(position='jitter',size=1)+
  facet_grid(. ~ region)+
  ggtitle("Notched violinplots showing salary distribution, separated by gender and discipline")

plotpervar<-function(featrs) {
  ln=length(featrs)
  plots.elems<-list()
  for(i in seq_along(featrs[-ln])) {
    plots.elems[[i]]<-
      ggplot(eda.df,aes_string(featrs[c(ln)],featrs[c(i)],colour="region"))+
      geom_point()+geom_rug()
  }

  return (plots.elems)
}

plot_grid(plotlist=plotpervar(intersect(names(eda.df), features)))

ggpairs(eda.df,mapping=ggplot2::aes(colour = region))


# Modelling -----------------

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

leapft<-do_model(features)
leapft

contrasts(fsel$region)

extractmodel<-function(mft, modnum) {
  return(list(mft[[1]][[modnum]],
              mft[[4]][[modnum]]
  ))
}

calcr2jack<-function(pressval) {
  SSyy=(nrow(fsel)-1)*var(fsel$TARGET_deathRate)
  return (1-(pressval/SSyy))
}

par(mfrow=c(2,2))
x<-extractmodel(leapft, 7)

x.selmod=x[[1]]
x.pressval=x[[2]]

plot(x.selmod)

x.pressval

calcr2jack(x.pressval)

par(mfrow=c(1,1))
results<-predict(x.selmod,fsel)
data.frame(results,actual=fsel$TARGET_deathRate) %>%
  ggplot(aes(x=results,y=actual)) + 
  geom_point()+
  stat_smooth(method="lm",show.legend = T)



