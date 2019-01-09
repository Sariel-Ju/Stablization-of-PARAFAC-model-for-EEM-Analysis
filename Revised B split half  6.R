rm(list=ls())
library(R.matlab)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Amelia)
library(parallel)
library(quadprog)
library(multiway)
library(gridExtra)
library(reshape2)
path <- ("D:\\Exchange program")
pathname <- file.path(path, "fluordata.mat")
testdata <- readMat(pathname)

#Set split of Rayleigh Scattering
split <- 10


#testdata   
#       $fluor (name, type, date, data, $label(8), $axisscale(8), title, 
#                 class, include, descriotion)
#       $ceonc (name, type, date, data, label, axisscale, title, class)
#fluo$data: 19(ex) * 136(em) * 5 * 405
#conc$data: 405 * 6}


#Sampling of specific data
sample <- which(testdata$conc$data[, 4] == 0 &
                testdata$conc$data[, 5] == 0)
test.temp <- testdata$fluor$data[sample,,,]
testmean <- (test.temp[,,,1] + test.temp[,,,2] + test.temp[,,,3] + 
         test.temp[,,,4] + test.temp[,,,5])/5
samp <- sample(1 : nrow(testmean), size = 30, replace = F)
testmean <- testmean[samp, , ]
em <- testdata$fluor$axisscale[2]
ex <- testdata$fluor$axisscale[3]
ex <- as.data.frame(ex)
em <- as.data.frame(em)


#Setting missing value and zero value
samex <- NULL
samem <- NULL
expos <- NULL
empos <- NULL
i <- 1
for(i in 1 : 19)
   {j <- 1
    for(j in 1 : 136)
       {if((em[1,j] - split <= ex[1, i] ) & ( em[1, j] + split >= ex[1, i]))
           {
              samex <- cbind(samex, ex[1, i])
              expos <- cbind(expos, i)
              samem <- cbind(samem, em[1, j])
              empos <- cbind(empos, j)
            }
        j <- j + 1
       }
    i <- i + 1
   }
i <- NULL
j <- NULL

raylwav <- rbind(samex, samem)
raylpos <- rbind(expos, empos)


i <- 1
for(i in 1 : (length(raylpos) / 2))
   {
      testmean[, (1: raylpos[2, i]), (raylpos[1, i] : 19)] <- 0
      i <- i + 1
   }
i <- 1
for(i in 1 : (length(raylpos) / 2))
   {
      testmean[, raylpos[2, i], raylpos[1, i]] <- NA
      i <- i + 1
   }
i <- NULL



#Imputing missing value
imputed <- array(, c(dim(testmean)[1],dim(testmean)[2],dim(testmean)[3]))
i <- 1
for(i in 1 : dim(testmean)[1]) 
   {  limit <- matrix(nrow = dim(testmean)[3], ncol = 1, c(1: dim(testmean)[3]))
      limit <- cbind(limit, 0, Inf)
      imp <- amelia(testmean[i,,], m = 1, bounds = limit)
      imputed[i,,] <- imp$imputations$imp1
      i <- i + 1
   }



#parafac on original data
nfac <- 4
parafac_nonn <- parafac(imputed, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
attach(parafac_nonn)
parafac_nonn <- rescale(parafac_nonn, mode = "B", 
                        newscale = 1/(nrow(B)^0.5), absorb = "A")
parafac_nonn <- rescale(parafac_nonn, mode = "C", 
                        newscale = 1/(nrow(C)^0.5), absorb = "A")
i <- 1
calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
for(i in 1: nfac)
   { calcu <- calcu + A[,i]%o%B[,i]%o%C[,i]
     i <- i + 1
   } 
i <- NULL
e <- imputed - calcu
detach(parafac_nonn)

nonn_em <- as.data.frame(parafac_nonn$B)
nonn_ex <- as.data.frame(parafac_nonn$C)
nonn_em <- cbind(t(em),nonn_em )
nonn_ex <- cbind(t(ex),nonn_ex )
names(nonn_em) <- c("wavel", 1 : nfac)
names(nonn_ex) <- c("wavel", 1 : nfac)
melt1 <- melt(nonn_em, id.vars = "wavel")
nonn_emimg <- ggplot(data = melt1, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Nonnegative emission spectra for original data", x = "Wavalength") 
melt2 <- melt(nonn_ex, id.vars = "wavel")
nonn_eximg <- ggplot(data = melt2, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Nonnegative excitation spectra for original data", x = "Wavalength") 


#Generation of new residuals
enew <- array(, c(dim(e)[1], dim(e)[2], dim(e)[3]))
i <- 1
for(i in 1 : nrow(e))
   {enew[i,,] <- matrix(
                        rnorm((dim(e)[2] * dim(e)[3]), 
                              mean = mean(e[i,,]), sd = sd(e[i,,])), 
                       nrow = dim(e)[2])
   i <- i + 1
   }

corA <- cor(parafac_nonn$A)
corB <- cor(parafac_nonn$B)
corC <- cor(parafac_nonn$C)

# 1. Revised components in C with original matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Anew <- parafac_nonn$A[, Acomp]
Bnew <- parafac_nonn$B[, Bcomp]
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Ccomp <- order
              Cnew <- parafac_nonn$C[, Ccomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("C", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\C", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\C\\C.csv")



# 2. Revised components in C with imputed B matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Anew <- parafac_nonn$A[, Acomp]
Bnew <- parafac_nonn$B[, Bcomp]
Bcomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$B)), 1, 0.5)
Bcomp4 <- Bcomp4/sqrt((sum(Bcomp4 ^ 2)))
Bnew <- cbind(parafac_nonn$B[, -4], Bcomp4)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
corA <- cor(Anew)
corB <- cor(Bnew)
corC <- cor(parafac_nonn$C)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Ccomp <- order
              Cnew <- parafac_nonn$C[, Ccomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("C", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\Bnew-C", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\Bnew-C\\Bnew-C.csv")





# 3. Revised components in C with imputed ABC matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Acomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$A)), 1, 0.5)
Acomp4 <- Acomp4/sqrt((sum(Acomp4 ^ 2)))*sqrt((sum(parafac_nonn$A[,4] ^ 2)))
Ag <- cbind(parafac_nonn$A[, -4], Acomp4)
Anew <- Ag
Bnew <- parafac_nonn$B[, Bcomp]
Bcomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$B)), 1, 0.5)
Bcomp4 <- Bcomp4/sqrt((sum(Bcomp4 ^ 2)))
Bnew <- cbind(parafac_nonn$B[, -4], Bcomp4)
Ccomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$C)), 1, 0.5)
Ccomp4 <- Ccomp4/sqrt((sum(Ccomp4 ^ 2)))
Cg <- cbind(parafac_nonn$C[, -4], Ccomp4)
corA <- cor(Anew)
corB <- cor(Bnew)
corC <- cor(Cg)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Ccomp <- order
              Cnew <- Cg[, Ccomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("C", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\ABCnew-C", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\ABCnew-C\\ABCnew-C.csv")





# 4. Revised components in C with imputed C matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Anew <- parafac_nonn$A[, Acomp]
Bnew <- parafac_nonn$B[, Bcomp]
Ccomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$C)), 1, 0.5)
Ccomp4 <- Ccomp4/sqrt((sum(Ccomp4 ^ 2)))
Cg <- cbind(parafac_nonn$C[, -4], Ccomp4)
corA <- cor(Anew)
corB <- cor(Bnew)
corC <- cor(Cg)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Ccomp <- order
              Cnew <- Cg[, Ccomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("C", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\Cnew-C", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\Cnew-C\\Cnew-C.csv")





# 5. Revised components in B with imputed C matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Anew <- parafac_nonn$A[, Acomp]
Ccomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$C)), 1, 0.5)
Ccomp4 <- Ccomp4/sqrt((sum(Ccomp4 ^ 2)))
Cnew <- cbind(parafac_nonn$C[, -4], Ccomp4)
corA <- cor(Anew)
corB <- cor(parafac_nonn$B[, Bcomp])
corC <- cor(Cnew)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Bcomp <- order
              Bnew <- parafac_nonn$B[, Bcomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("B", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\Cnew-B", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Cnew[, p], calcu_nonn$C[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\Cnew-B\\Cnew-B.csv")







# 6. Revised components in B with imputed B matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Ccomp <- c(1, 2, 3, 4)
Bnew <- parafac_nonn$B[, Bcomp]
Bcomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$B)), 1, 0.5)
Bcomp4 <- Bcomp4/sqrt((sum(Bcomp4 ^ 2)))
Bg <- cbind(parafac_nonn$B[, -4], Bcomp4)
Anew <- parafac_nonn$A[, Acomp]
Cnew <- parafac_nonn$C[, Ccomp]
corA <- cor(Anew)
corB <- cor(Bg)
corC <- cor(Cnew)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Bcomp <- order
              Bnew <- Bg[, Bcomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("B", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\Bnew-B", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Cnew[, p], calcu_nonn$C[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\Bnew-B\\Bnew-B.csv")





# 7. Revised components in B with imputed ABC matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Acomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$A)), 1, 0.5)
Acomp4 <- Acomp4/sqrt((sum(Acomp4 ^ 2)))*sqrt((sum(parafac_nonn$A[,4] ^ 2)))
Ag <- cbind(parafac_nonn$A[, -4], Acomp4)
Anew <- Ag
Bnew <- parafac_nonn$B[, Bcomp]
Bcomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$B)), 1, 0.5)
Bcomp4 <- Bcomp4/sqrt((sum(Bcomp4 ^ 2)))
Bg <- cbind(parafac_nonn$B[, -4], Bcomp4)
Ccomp4 <- dnorm(seq(0,5,length.out=nrow(parafac_nonn$C)), 1, 0.5)
Ccomp4 <- Ccomp4/sqrt((sum(Ccomp4 ^ 2)))
Cnew <- cbind(parafac_nonn$C[, -4], Ccomp4)
corA <- cor(Anew)
corB <- cor(Bg)
corC <- cor(Cnew)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Bcomp <- order
              Bnew <- Bg[, Bcomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("B", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\ABCnew-B", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Cnew[, p], calcu_nonn$C[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\ABCnew-B\\ABCnew-B.csv")





# 8. Revised components in B with original matrices.


standard <- 0.95
table <- matrix(, nrow = 10)
Acomp <- c(1, 2, 3, 4)
Bcomp <- c(1, 2, 3, 4)
Anew <- parafac_nonn$A[, Acomp]
Bg <- parafac_nonn$B[, Bcomp]
Cnew <- parafac_nonn$C[, Ccomp]
corA <- cor(Anew)
corB <- cor(Bg)
corC <- cor(Cnew)
sample <- sample(1 : nrow(calcu), size = nrow(calcu) / 2)
i <- 1
for(i in 1 : 4)
{  j <- 1
   for(j in 1 : 4)
      {order <- c(1, 2, 3, 4)
       if(i != j)          
             {print(i)
              print(j)
              order[j] <- order[i]
              Bcomp <- order
              Bnew <- Bg[, Bcomp]
              n <- 1
              calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Entire Revised Data", x = "Wavelength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Entire Revised Data", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Emission Spectra for Part 1 of Revised Data", x = "Wavelength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "Excitation Spectra for Part 1 of Revised Data", x = "Wavelength")
              name <- paste("B", order[1], order[2], order[3], order[4], sep = "")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), path = "D:\\Exchange program\\Final\\B", p)
              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Cnew[, p], calcu_nonn$C[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 5)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[i, j]
              temp[6] <- corB[i, j]
              temp[7] <- corC[i, j]
              temp[8] <- as.numeric(temp[5]) * as.numeric(temp[6]) * as.numeric(temp[7])
              temp[9] <- abs(as.numeric(temp[5])) + abs(as.numeric(temp[6])) + abs(as.numeric(temp[7]))
              temp[10] <- status
              table <- cbind(table, temp)
           
             }
              j <- j + 1
         } 
    i <- i + 1
}
write.csv(table, "D:\\Exchange program\\Final\\B\\B.csv")



#9. Peak 1 of C's position changing with stabled exchanged B of 1231
Anew <- parafac_nonn$A
Border <- c(1, 2, 3, 1)
Bnew <- parafac_nonn$B[, Border]
corA <- cor(Anew)
corB <- cor(Bnew)
Ccall <- array(NA, 28)
Ccall[6: 24] <- parafac_nonn$C[, 1]
table <- matrix(, nrow = 10)
i <- 1
for(i in 1 : 10)
{   print(i)
    Ccomp <- Ccall[i: (i+ 18)]
    Cnew <- parafac_nonn$C
    Ccomp <- na.spline(Ccomp)
    Ccomp[Ccomp < 0] <- 0
    Cnew[,1] <- Ccomp
    corC <- cor(Cnew)
    n <- 1
    calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
    
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative emission spectra", x = "Wavalength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative excitation spectra", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 emission spectra", x = "Wavalength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 onnegative excitation spectra", x = "Wavalength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 emission spectra", x = "Wavalength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 excitation spectra", x = "Wavalength")
              name <- paste("B at position", i, sep = " ")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), 
              path = "D:\\Exchange program\\Final\\Exchanged B of 1231", p)

              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 10)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[1, 4]
              temp[6] <- corB[1, 4]
              temp[7] <- corC[1, 4]
              temp[8] <- corC[2, 1]
              temp[9] <- corC[3, 1]
              temp[10] <- status
              table <- cbind(table, temp)
              
              i <- i + 1
}
    write.csv(table, "D:\\Exchange program\\Final\\Exchanged B of 1231\\Exchanged B of 1231.csv")





#10. Peak 2 of C's position changing with stabled exchanged B of 1134
Anew <- parafac_nonn$A
Border <- c(1, 1, 3, 4)
Bnew <- parafac_nonn$B[, Border]
corA <- cor(Anew)
corB <- cor(Bnew)
Ccall <- array(NA, 28)
Ccall[6: 24] <- parafac_nonn$C[, 2]
table <- matrix(, nrow = 10)
i <- 1
for(i in 1 : 10)
{   print(i)
    Ccomp <- Ccall[i: (i+ 18)]
    Cnew <- parafac_nonn$C
    Ccomp <- na.spline(Ccomp)
    Ccomp[Ccomp < 0] <- 0
    Cnew[,2] <- Ccomp
    corC <- cor(Cnew)
    n <- 1
    calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
    
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative emission spectra", x = "Wavalength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative excitation spectra", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 emission spectra", x = "Wavalength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 onnegative excitation spectra", x = "Wavalength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 emission spectra", x = "Wavalength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 excitation spectra", x = "Wavalength")
              name <- paste("B at position", i, sep = " ")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), 
              path = "D:\\Exchange program\\Final\\Exchanged B of 1134", p)

              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 10)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[1, 4]
              temp[6] <- corB[1, 4]
              temp[7] <- corC[1, 4]
              temp[8] <- corC[2, 1]
              temp[9] <- corC[3, 1]
              temp[10] <- status
              table <- cbind(table, temp)
              
              i <- i + 1
}
    write.csv(table, "D:\\Exchange program\\Final\\Exchanged B of 1134\\Exchanged B of 1134.csv")





#11. Peak 2 of C's position changing with stabled exchanged B of 3234
Anew <- parafac_nonn$A
Border <- c(3, 2, 3, 4)
Bnew <- parafac_nonn$B[, Border]
corA <- cor(Anew)
corB <- cor(Bnew)
Ccall <- array(NA, 28)
Ccall[6: 24] <- parafac_nonn$C[, 3]
table <- matrix(, nrow = 10)
i <- 1
for(i in 1 : 10)
{   print(i)
    Ccomp <- Ccall[i: (i+ 18)]
    Cnew <- parafac_nonn$C
    Ccomp <- na.spline(Ccomp)
    Ccomp[Ccomp < 0] <- 0
    Cnew[,3] <- Ccomp
    corC <- cor(Cnew)
    n <- 1
    calcu <- array(0, c(nrow(testmean), length(em), length(ex)))
              for(n in 1: nfac)
                 {calcu <- calcu + Anew[,n]%o%Bnew[,n]%o%Cnew[,n]
                  n <- n + 1
                 } 
              n <- NULL
              calcu <- calcu + enew
    
              calcu_nonn <- parafac(calcu, nfac = nfac, nstart = 30, 
                      const = c(0, 2, 2))
              calcu_nonn <- rescale(calcu_nonn, mode = "B", 
                            newscale = 1/(nrow(Bnew)^0.5), absorb = "A")
              calcu_nonn <- rescale(calcu_nonn, mode  = "C", 
                            newscale = 1/(nrow(Cnew)^0.5), absorb = "A")
              calcu_em <- as.data.frame(calcu_nonn$B)
              calcu_ex <- as.data.frame(calcu_nonn$C)
              calcu_em <- cbind(t(em),calcu_em )
              calcu_ex <- cbind(t(ex),calcu_ex )
              names(calcu_em) <- c("wavel", 1 : nfac)
              names(calcu_ex) <- c("wavel", 1 : nfac)

              melt3 <- melt(calcu_em, id.vars = "wavel")

              calcu_emimg <- ggplot(data = melt3, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative emission spectra", x = "Wavalength") 
              melt4 <- melt(calcu_ex, id.vars = "wavel")

              calcu_eximg <- ggplot(data = melt4, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "calcu nonnegative excitation spectra", x = "Wavalength") 

              #Split-half of enew
              calcu_sp1 <- calcu[sample,,]
              calcu_sp2 <- calcu[-sample,,]
              calcu1 <- parafac(calcu_sp1, nfac = nfac, nstart = 20, 
                        const = c(0, 2, 2))
              calcu1 <- rescale(calcu1, mode = "B", 
                        newscale = 1/(nrow(calcu1$B)^0.5), absorb = "A")
              calcu1 <- rescale(calcu1, mode = "C", 
                        newscale = 1/(nrow(calcu1$C)^0.5), absorb = "A")
              calcu2 <- parafac(calcu_sp2, nfac = nfac, nstart = 20, 
                          const = c(0, 2, 2))
              calcu2 <- rescale(calcu2, mode = "B", 
                        newscale = 1/(nrow(calcu2$B)^0.5), absorb = "A")
              calcu2 <- rescale(calcu2, mode = "C", 
                        newscale = 1/(nrow(calcu2$C)^0.5), absorb = "A")
              sp1_em <- as.data.frame(calcu1$B)
              sp1_ex <- as.data.frame(calcu1$C)
              sp1_em <- cbind(t(em),sp1_em )
              sp1_ex <- cbind(t(ex),sp1_ex )
              names(sp1_em) <- c("wavel", 1 : nfac)
              names(sp1_ex) <- c("wavel", 1 : nfac)

              sp2_em <- as.data.frame(calcu2$B)
              sp2_ex <- as.data.frame(calcu2$C)
              sp2_em <- cbind(t(em),sp2_em )
              sp2_ex <- cbind(t(ex),sp2_ex )
              names(sp2_em) <- c("wavel", 1 : nfac)
              names(sp2_ex) <- c("wavel", 1 : nfac)

              spmelt1_em <- melt(sp1_em, id.vars = "wavel")
              sp1_emimg <- ggplot(data = spmelt1_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 emission spectra", x = "Wavalength") 
              spmelt1_ex <- melt(sp1_ex, id.vars = "wavel")
              sp1_eximg <- ggplot(data = spmelt1_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp1 onnegative excitation spectra", x = "Wavalength") 
              spmelt2_em <- melt(sp2_em, id.vars = "wavel")
              sp2_emimg <- ggplot(data = spmelt2_em, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 emission spectra", x = "Wavalength") 
              spmelt2_ex <- melt(sp2_ex, id.vars = "wavel")
              sp2_eximg <- ggplot(data = spmelt2_ex, aes(x = wavel, y = value, color = variable)) + 
                   geom_line(size = 1) +
                   labs(title = "sp2 excitation spectra", x = "Wavalength")
              name <- paste("B at position", i, sep = " ")
              p <- grid.arrange(nonn_emimg, nonn_eximg,  calcu_emimg,  calcu_eximg, sp1_emimg, sp1_eximg, 
                                sp2_emimg, sp2_eximg, ncol = 2)
              ggsave(file = paste(name, ".pdf", sep = "", collapse = ""), 
              path = "D:\\Exchange program\\Final\\Exchanged B of 3234", p)

              corr <- matrix(, ncol = 4, nrow = 4)
              p <- 1
              for(p in 1: 4){
                  q <- 1
                  for(q in 1: 4){
                     corr[p, q] <- cor(Bnew[, p], calcu_nonn$B[, q])
                     q <- q + 1      }
                     p <- p + 1      }
              p <- 1
              temp <- array(, 4)
              for(p in 1: 4){
                  temp1 <- corr[p, corr[p,] > standard]
                  temp2 <- corr[corr[, p] > standard, p]
                  if((length(temp1) == 1)&(length(temp2) == 1))
                  temp[p] <- TRUE
                  else temp[p] <- FALSE
                  p <- p + 1
                             }
              status <- all(temp)
              temp <- matrix(, nrow = 10)
              temp[1] <- name
              temp[2] <- corcondia(imputed, calcu_nonn, divisor = "core")/100
              temp[3] <- corcondia(calcu_sp1, calcu1, divisor = "core")/100
              temp[4] <- corcondia(calcu_sp2, calcu2, divisor = "core")/100
              temp[5] <- corA[1, 4]
              temp[6] <- corB[1, 4]
              temp[7] <- corC[1, 4]
              temp[8] <- corC[2, 1]
              temp[9] <- corC[3, 1]
              temp[10] <- status
              table <- cbind(table, temp)
              
              i <- i + 1
}
    write.csv(table, "D:\\Exchange program\\Final\\Exchanged B of 3234\\Exchanged B of 3234.csv")



