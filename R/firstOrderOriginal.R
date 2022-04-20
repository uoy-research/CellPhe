## codes for first order variables

FOmean<-vector(length = length(mini_image))
FOsd<-vector(length = length(mini_image))
FOskew<-vector(length = length(mini_image))

for(i in c(1:length(mini_image)))
{
  justintensities<-subset(mini_image[[i]][c(4:length(mini_image[[i]])+3)], mini_image[[i]][c(4:length(mini_image[[i]])+3)] != -1)
  FOmean[i] = mean(justintensities)
  FOsd[i] = (var(justintensities))^(1/2)
  FOskew[i] = e1071::skewness(justintensities)
}



