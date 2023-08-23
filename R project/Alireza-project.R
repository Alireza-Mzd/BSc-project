#Statistics and Probabilities Project - Dr. Lorvand
#Alireza Mohammadzadeh
#Mehdi Jafarian

data_petrol=read.csv("D:\\Alireza\\E-Classes\\401-1\\statistics and Probabibilities\\R class 1401_2\\AlirezaProject\\second data\\APdata.csv",header=TRUE)
str(data_petrol)

head(data_petrol,5)
tail(data_petrol,5)
#part1.1

#it kept giving errors that (as ‘lib’ is unspecified) therefore I wrote all of them
library(foreign)
library(survival)
library(MASS)
library(nnet)
install.packages("epiDisplay")
library(epiDisplay)
tab1(data_petrol$Grid, sort.group = "increasing", cum.percent = TRUE)
tab1(data_petrol$FluidType, sort.group = "increasing", cum.percent = TRUE)
#part 1.2

#summary statistics 
summary(data_petrol)

#Mode
install.packages("modeest")
library(modeest)
mfv(data_petrol$Porosity)
mfv(data_petrol$Viscosity)


#median
median(data_petrol$Porosity)
median(data_petrol$Viscosity)

#mean
mean(data_petrol$Porosity)
mean(data_petrol$Viscosity)

#Range 
range(data_petrol$Porosity)
range(data_petrol$Viscosity)

#variance
var(data_petrol$Porosity)
var(data_petrol$Viscosity)

#standard deviation
sd(data_petrol$Porosity)
sd(data_petrol$Viscosity)

#part 1.3
#graphing the data
Grid=table(data_petrol$Grid)
barplot(Grid)
hist(data_petrol$Viscosity)
hist(data_petrol$Porosity)
boxplot(data_petrol$Porosity)
boxplot(data_petrol$Viscosity)
grid=table(data_petrol$Grid)
pie(grid)
type=table(data_petrol$FluidType)
pie(type)

#part 1.4
#covariance
cov(data_petrol$Porosity, data_petrol$Viscosity )

#correlation coefficient 
cor(data_petrol$Porosity, data_petrol$Viscosity )


#part 2.1 
sample.mean=mean(data_petrol$Porosity)
print(sample.mean)

sample.n=length(data_petrol$Porosity)
sample.sd=sd(data_petrol$Porosity)
sample.se=sample.sd/sqrt(sample.n)
print(sample.se)

alpha= 0.05
degrees.freedom=sample.n-1
t.score=qt(p=1-(alpha/2), df=degrees.freedom, lower.tail = TRUE)
margin.error=t.score*sample.se
L= sample.mean - margin.error
U= sample.mean+margin.error
print((c(L, U)))

#part 2.2
sample.var = var(data_petrol$Porosity)
print(sample.var)

sample.n1 = length(data_petrol$Porosity)
alpha = 0.05
degrees.freedom = sample.n1 - 1
L = qchisq(1-(alpha/2),df=degrees.freedom,lower.tail = TRUE)
U = qchisq(alpha/2,df=degrees.freedom,lower.tail = TRUE)
l = (degrees.freedom - sample.var)/L
u = (degrees.freedom - sample.var)/U
print(c(l,u))

#part 2.3
alpha=0.05
sample.mean=mean(data_petrol$Porosity)
sample.n=length(data_petrol$Porosity)
sample.sd=sd(data_petrol$Porosity)
sample.se=sample.sd/sqrt(sample.n)
print(sample.se)

t.stat=(sample.mean-10)/sample.se
degrees.freedom=sample.n-1
t.score=qt(p=1-alpha, df=degrees.freedom, lower.tail = TRUE)
t.stat
t.score

#saving graphs and plots
getwd()
setwd("D:/Alireza/E-Classes/401-1/statistics and Probabibilities/R class 1401_2/AlirezaProject")
#however, I used 'export tab' in the plots section to save them.

