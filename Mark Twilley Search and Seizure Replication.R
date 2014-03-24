

setwd("~/R/search and seizure")

#  Remove the stored objects and load my libraries

rm(list = ls())     # remove all stored objects   

library(rpart) #library for running trees
library(Rcpp) # mark says this is a required package for "arm" to run
library(MASS) # mark says theese are required packages
library(Matrix) # mark says this is a required package for "arm" to run
library(lme4) # mark says this is a required package for "arm" to run
library(arm) #for using "attach.all" function
library(foreign)    # allow for import of STATA datset





#added a space between the arrow and read.dta
search.court <- read.dta("search_data_replication.dta", convert.underscore = T)
search.court$except <- ifelse(search.court$except >= 1, 1,0) #turn except count variable into binary 0,1

search.court$sctdec <- ifelse(search.court$sctdec == 1, 0, 1)#switch DV so that 1 is a liberal decision

#consistent w/ Benesh data
search.court$dec <- factor(search.court$sctdec, levels=0:1, labels=c("R", "U")) #use this factor for DV 

# 0 = reasonable
#1 = unreasonable




#The code didnt originally provide an attach.all function so I found this online and it seems to work
attach.all <- function (x, overwrite = NA, name = "attach.all")  {
  rem <- names(x) %in% ls(.GlobalEnv)
  if (!any(rem)) overwrite <- FALSE
  rem <- names(x)[rem]
  if (is.na(overwrite)) {
    question <- paste("The following objects in .GlobalEnv will mask\nobjects in the attached database:\n", paste(rem, collapse = ", "), "\nRemove these objects from .GlobalEnv?", sep = "")
    if (interactive()) {
      if (.Platform$OS.type == "windows")  overwrite <- "YES" == winDialog(type = "yesno",  question)
      else overwrite <- 1 == menu(c("YES", "NO"), graphics = FALSE, title = question)
    }
    else overwrite <- FALSE
  }
  if (overwrite) remove(list = rem, envir = .GlobalEnv)
  attach(x, name = name)
}


#Apply my fancy new attach.all function
attach.all(search.court)



#First break in the provided code; The next section makes, prunes, then cross validates the first tree



keep <- term <= 72
serc1<-rpart(dec~lctinc + lctaft + warrant + house + person + business + car + full + 
               except, cp=.01, subset = keep, minsplit= 15)
printcp(serc1, digits = 2) #print cp table
post(serc1, file = "") #plot tree in graphics window

#Prune Tree - use .p to indicate pruned trees
serc1.p <- prune(serc1, cp = .02)
printcp(serc1.p, digits = 2)

#PLOT TREE 
post(serc1.p, file = "")#note: numbers under nodes need to be reversed to present

#modal category 1st 
# the default is level 1/level 2, i.e. R/U
#Resubstitution Predictions
re.serc1 <- ifelse((predict(serc1.p, type="vector")) == 2, 1, 0)#Resubsitution prediction

#Recode prediction vector to 0,1 not 1,2
table(re.serc1,dec[keep])

#Proportional reduction in error
#100 x [(% correctly predicted - % in modal category )/(100%-% in modal category.]
modal.cat.serc1 <- ifelse((length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]) > .5, #since modal cat is 0
                          (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]),                            #need ifelse command
                          1 - (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]))                       
cor.re.serc1 <- ifelse(re.serc1==sctdec[keep],1,0)
pred.cor.re.serc1  <- sum(cor.re.serc1)/length(sctdec[keep])
pred.cor.re.serc1 #Resubstitution prediction rate
pre.re.serc1 <- 100*((100*pred.cor.re.serc1 - 100*modal.cat.serc1)/(100-100*modal.cat.serc1))
pre.re.serc1 #Resubstitution PRE

#Cross-Validated Prediction
#simulate 1000 times and take mean
n.sims <- 100
pre.xv.serc1 <- rep(NA, length(n.sims))
pred.correct <- rep(NA, length(n.sims))
for (i in 1:n.sims){
  xv.serc1 <- ifelse((xpred.rpart(serc1.p)) == 2, 1,0) #cross-validated prediction; 
  #recode prediction vector to 0,1 not 1,2
  correct.xv.serc1 <- ifelse(xv.serc1==sctdec[keep],1,0)
  pred.cor.xv.serc1  <- sum(correct.xv.serc1)/length(sctdec[keep])
  pred.correct[i] <- pred.cor.xv.serc1
  pre.xv.serc1[i] <- 100*((100*pred.cor.xv.serc1 - 100*modal.cat.serc1)/(100-100*modal.cat.serc1))
}
print(summary(pred.correct)) # Cross-validated prediction rate 
print(summary(pre.xv.serc1)) # Cross-validated PRE

############################################
############################################
#Break in the code for Figure 6
############################################
############################################

#create a property variable and run tree just on full search and property

property.dummy <- ifelse(house == 1 | car==1| business == 1| person ==1, 1,0)
serc1.prop<-rpart(dec~property.dummy + full, cp=.0001, subset = keep)
printcp(serc1.prop, digits = 2)
post(serc1.prop, file = "")
re.serc.property <- ifelse((predict(serc1.prop, type="vector")) == 2, 1, 0)#resubstitution rate

#plot partitition with actual cases statistics
#get counts for each "region" -- starting w/ bottom left, and going clockwise
keep <- term<=72
region.1 <- ifelse(full==0 & property.dummy==0 & keep, 1,0)#identify cases in region 1 
sum.actual.1 <- sum(region.1) #get total number of cases in region 1 in 1st period, using [keep]
correct.region.1 <- ifelse(re.serc.property[region.1==1] == sctdec[region.1==1 & keep], 1,0)
sum.correct.region.1 <- sum(correct.region.1)

region.2 <- ifelse(full==1 & property.dummy==0 & keep, 1,0)#identify cases in region 2
sum.actual.2 <- sum(region.2) 
correct.region.2 <- ifelse(re.serc.property[region.2==1] == sctdec[region.2==1 & keep], 1,0)
sum.correct.region.2 <- sum(correct.region.2)

region.3 <- ifelse(full==1 & property.dummy==1 & keep, 1,0)#identify cases in region 3 
sum.actual.3 <- sum(region.3) 
correct.region.3 <- ifelse(re.serc.property[region.3==1] == sctdec[region.3==1 & keep], 1,0)
sum.correct.region.3 <- sum(correct.region.3)

region.4 <- ifelse(full==0 & property.dummy==1 & keep, 1,0)#identify cases in region 4 
sum.actual.4 <- sum(region.4) 
correct.region.4 <- ifelse(re.serc.property[region.4==1] == sctdec[region.4==1 & keep], 1,0)
sum.correct.region.4 <- sum(correct.region.4)

print(sum.actual.1)
print(sum.correct.region.1)
print(sum.actual.2)
print(sum.correct.region.2)
print(sum.actual.3)
print(sum.correct.region.3)
print(sum.actual.4)
print(sum.correct.region.4)

pdf("Figure_6.pdf", height = 4, width = 4)

par(mfrow = c(1,1), mar=c(5,4,1,1), pty = "s")#make plot square
plot(0, 0, xlab="",ylab="", type="n", axes = F, 
     xaxs = "i", yaxs = "i", xlim = c(0,10), ylim = c(0,10)) # call empty plot
box()
axis(2, at = c(2.5,7.5), tick = T, label = c("No","Yes"), las = 1)
mtext("Full Search?", 2, line =2.5, cex = 1.2, font = 2)
axis(1, at = c(2.5,7.5), tick = T, label = c("No","Yes"), las = 1)
mtext("Property Interest?", 1, line = 2.5, cex = 1.2, font = 2)
#indicate exclude area with shading
polygon(x=c(5, 10, 10,5),
        y = c(5,5,10,10), #left-hand block of "not liable" region
        col= "gray90",  border=T)
abline(v=5, lty = 2)
abline(h=5, lty =2 )
#add counts for each region: first number of cases that tree correctly predicts,
#then number of cases that fall into each region
text(2.5,2.5, "0 / 0")
text(2.5,7.5, "5 / 7")
text(7.5,7.5, "33 / 49")
text(7.5,2.5, "8 / 9")

dev.off()

################################################################
################################################################
#Figure 7: Period 2, 1962-1976 terms
################################################################
################################################################

detach()
attach.all(search.court)
keep <- term<=76


#.0000001
serc2<-rpart(dec~lctinc + lctaft + warrant + house + person + business + car + full + except, 
             cp=.0000001, subset = keep, minsplit = 15)
printcp(serc2, digits = 2)
post(serc2, file = "")

#Prune Tree - use .p to indicate pruned trees   .04
serc2.p <- prune(serc2, cp = .01)
printcp(serc2.p, digits =2)

#PLOT TREE
post(serc2.p, file = "") #Note: tree is flipped in paper for consistency with Figure 5
#i.e. no and yes split in the same direction from 1st node

#Resubstitution Predictions
re.serc2 <- ifelse((predict(serc2.p, type="vector")) == 2, 1, 0)#Resubsitution prediction

#Recode prediction vector to 0,1 not 1,2
table(re.serc2,dec[keep])

#Proportional reduction in error
#100 x [(% correctly predicted - % in modal category )/(100%-% in modal category.]
modal.cat.serc2 <- ifelse((length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]) > .5, #since modal cat is 0
                          (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]),                            #need ifelse command
                          1 - (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]))                       
cor.re.serc2 <- ifelse(re.serc2==sctdec[keep],1,0)
pred.cor.re.serc2  <- sum(cor.re.serc2)/length(sctdec[keep])
pred.cor.re.serc2 #Resubstitution prediction rate
pre.re.serc2 <- 100*((100*pred.cor.re.serc2 - 100*modal.cat.serc2)/(100-100*modal.cat.serc2))
pre.re.serc2  #Resubstitution PRE

#Cross-Validated Prediction
#simulate 1000 times and take mean
n.sims <- 100
pre.xv.serc2 <- rep(NA, length(n.sims))
pred.correct <- rep(NA, length(n.sims))
for (i in 1:n.sims){
  xv.serc2 <- ifelse((xpred.rpart(serc2.p)) == 2, 1, 0) #cross-validated prediction; 
  #recode prediction vector to 0,1 not 1,2
  #table(xv.serc2,dec[keep])
  #xv.serc2 Ex Ad
  #       0 24 18
  #       1 12 11
  #Proportional reduction in error
  correct.xv.serc2 <- ifelse(xv.serc2==sctdec[keep],1,0)
  pred.cor.xv.serc2  <- sum(correct.xv.serc2)/length(sctdec[keep])
  pred.correct[i] <- pred.cor.xv.serc2
  pre.xv.serc2[i] <- 100*((100*pred.cor.xv.serc2 - 100*modal.cat.serc2)/(100-100*modal.cat.serc2))
}
print(summary(pred.correct)) #Cross-validated prediction rate 
print(summary(pre.xv.serc2)) #Cross-validated PRE

################################################################
################################################################
#Figure 8
################################################################
################################################################

property.dummy <- ifelse(house == 1 | car==1| business == 1| person ==1, 1,0)
serc2.prop<-rpart(dec~property.dummy + except, cp=.01, subset = keep)
printcp(serc2.prop, digits = 2)
post(serc2.prop, file = "")
re.serc.property <- ifelse((predict(serc2.prop, type="vector")) == 2, 1, 0)#resubstitution rate

#plot partitition with actual cases statistics
#get counts for each "region" -- starting w/ bottom left, and going clockwise
keep <- term<=76
region.1 <- ifelse(except==0&property.dummy==0 & keep, 1,0)#identify cases in region 1 
sum.actual.1 <- sum(region.1) #get total number of cases in region 1 in 1st period, using [keep]
correct.region.1 <- ifelse(re.serc.property[region.1==1] == sctdec[region.1==1 & keep], 1,0)
sum.correct.region.1 <- sum(correct.region.1)

region.2 <- ifelse(except==1&property.dummy==0 & keep, 1,0)#identify cases in region 1 
sum.actual.2 <- sum(region.2) 
correct.region.2 <- ifelse(re.serc.property[region.2==1] == sctdec[region.2==1 & keep], 1,0)
sum.correct.region.2 <- sum(correct.region.2)

region.3 <- ifelse(except==1&property.dummy==1 & keep, 1,0)#identify cases in region 1 
sum.actual.3 <- sum(region.3) 
correct.region.3 <- ifelse(re.serc.property[region.3==1] == sctdec[region.3==1 & keep], 1,0)
sum.correct.region.3 <- sum(correct.region.3)

region.4 <- ifelse(except==0&property.dummy==1 & keep, 1,0)#identify cases in region 1 
sum.actual.4 <- sum(region.4) 
correct.region.4 <- ifelse(re.serc.property[region.4==1] == sctdec[region.4==1 & keep], 1,0)
sum.correct.region.4 <- sum(correct.region.4)

print(sum.actual.1)
print(sum.correct.region.1)
print(sum.actual.2)
print(sum.correct.region.2)
print(sum.actual.3)
print(sum.correct.region.3)
print(sum.actual.4)
print(sum.correct.region.4)

pdf("Figure_8.pdf", height = 4, width = 4)

par(mfrow = c(1,1), mar=c(5,4,1,1), pty = "s")#make plot square
plot(0, 0, xlab="",ylab="", type="n", axes = F, 
     xaxs = "i", yaxs = "i", xlim = c(0,10), ylim = c(0,10)) # call empty plot
box()
axis(2, at = c(2.5,7.5), tick = T, label = c("No","Yes"), las = 1)
mtext("Exception?", 2, line =2.5, cex = 1.2, font =2)
axis(1, at = c(2.5,7.5), tick = T, label = c("No","Yes"), las = 1)
mtext("Property Interest?", 1, line = 2.5, cex = 1.2, font =2)
#indicate exclude area with shading
polygon(x=c(5, 10, 10,5),
        y = c(5,5,0,0), #left-hand block of "not liable" region
        col= "gray90",  border=T)
abline(v=5, lty = 2)
abline(h=5, lty =2 )
#add counts for each region: first number of cases that tree correctly predicts,
#then number of cases that fall into each region
text(2.5,2.5, "5 / 7")
text(2.5,7.5, "1 / 1")
text(7.5,7.5, "20 / 26")
text(7.5,2.5, "33 / 55")

dev.off()

###############################################################################################
#Figure 9, Period 3: 1962-1985 terms
#####################################################################################


detach()
attach.all(search.court)
keep <- term<=85

serc3<-rpart(dec~lctinc + lctaft + warrant + house + person + business + car + full + except, 
             cp=.01, subset = keep, minsplit = 15)
printcp(serc3, digits = 2)
post(serc3, file = "")

#Prune Tree - use .p to indicate pruned trees
serc3.p <- prune(serc3, cp = .02)
printcp(serc3.p, digits = 2)
#Plot Tree
post(serc3.p, file = "") ##Note: tree is flipped in paper for consistency with Figure 5
#i.e. no and yes split in the same direction from 1st node
#Resubstitution Predictions
re.serc3 <- ifelse((predict(serc3.p, type="vector")) == 2, 1, 0)#Resubsitution prediction
#Recode prediction vector to 0,1 not 1,2
table(re.serc3,dec[keep])# = 70% correct

#Proportional reduction in error
#100 x [(% correctly predicted - % in modal category )/(100%-% in modal category.]
modal.cat.serc3 <- ifelse((length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]) > .5, #since modal cat is 0
                          (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]),                            #need ifelse command
                          1 - (length(sctdec[keep])-sum(sctdec[keep]))/length(sctdec[keep]))                       
cor.re.serc3 <- ifelse(re.serc3==sctdec[keep],1,0)
pred.cor.re.serc3  <- sum(cor.re.serc3)/length(sctdec[keep])
pred.cor.re.serc3 #Resubstitution prediction rate
pre.re.serc3 <- 100*((100*pred.cor.re.serc3 - 100*modal.cat.serc3)/(100-100*modal.cat.serc3))
pre.re.serc3 # Resubstitution PRE

#Cross-Validated Prediction
#simulate 1000 times and take mean
n.sims <- 1000
pre.xv.serc3 <- rep(NA, length(n.sims))
pred.correct <- rep(NA, length(n.sims))


for (i in 1:n.sims){
  xv.serc3 <- ifelse((xpred.rpart(serc3.p)) == 2, 1, 0) 
  #Proportional reduction in error
  correct.xv.serc3 <- ifelse(xv.serc3==sctdec[keep],1,0)
  pred.cor.xv.serc3  <- sum(correct.xv.serc3)/length(sctdec[keep])
  pred.correct[i] <- pred.cor.xv.serc3
  pre.xv.serc3[i] <- 100*((100*pred.cor.xv.serc3 - 100*modal.cat.serc3)/(100-100*modal.cat.serc3))
}



print(summary(pred.correct)) #Cross-validation prediction rate
print(summary(pre.xv.serc3)) #Cross-validation PRE






