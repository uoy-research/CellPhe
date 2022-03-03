data <- read.csv("/Users/laurawiggins/Desktop/05062019_B3/05062019_B3_3_Phase-FullFeatureTable.csv", header=FALSE)
cellID = 175
data<-data[-1,]
names(data)<-as.matrix(data[1,])
data<-data[-1,]
row.names(data)<-NULL

featuretable<-subset(data, data[,3]==cellID)
featuretable[,1]<-as.numeric(as.character(featuretable[,1]))
numframes<-(max(featuretable[,1])-min(featuretable[,1])+2)

##feature table output for cell cellID is featuretable
##numframes is the total number of frames for cell cellID, missing frames accounted for

roi<-RImageJROI::read.ijzip("/Users/laurawiggins/Desktop/05062019_B3/05062019_B3_3_Phase.zip")
roiData<-list()
for(j in c(1:length(roi)))
{
  roiData[[j]]<-list()
  roiData[[j]][[1]]<-strsplit(roi[[j]]$name, "-")[[1]][1]
  roiData[[j]][[2]]<-strsplit(roi[[j]]$name, "-")[[1]][2]
  roiData[[j]][[3]]<-roi[[j]]$coords
  print(j)
}

justcell<-subset(roiData,lapply(roiData,'[[', 2)==cellID)

maxboundarylength=0
for(i in c(1:length(justcell)))
{
  if(length(justcell[[i]][[3]][,1]) > maxboundarylength) 
    maxboundarylength = length(justcell[[i]][[3]][,1])
}

##maxboundarylength is the maximum boundary length for cell cellID

boundarycoords<-matrix(ncol = (maxboundarylength*2)+3, nrow = length(justcell))

for(i in c(1:length(justcell)))
{
  print(i)
  boundarycoords[i,1] = justcell[[i]][[1]]
  boundarycoords[i,2] = justcell[[i]][[2]]
  boundarycoords[i,3] = dim(justcell[[i]][[3]])[1]
  index = 4
  for(j in c(1:length(justcell[[i]][[3]][,1])))
  {
    print(j-1)
    print(j)
    boundarycoords[i,index] = justcell[[i]][[3]][[j,1]]
    boundarycoords[i,index+1] = justcell[[i]][[3]][j,2]
    index = index + 2
  }
}

##boundarycoords is the table to be used as the boundary file for cell cellID

## function to produce miniImage
storeInteriorIntensities<-function(whichcell, images, frame)
{
  intensities<-matrix(ncol = 3)
  orderedintensities<-vector()
  k = 0
  ##only keep interior coordinates
  for(i in c(min(whichcell[[3]][,1]):max(whichcell[[3]][,1])))
  {
    for(j in c(min(whichcell[[3]][,2]):max(whichcell[[3]][,2])))
    {
      if(ptinpoly::pip2d(whichcell[[3]],cbind(i,j)) == 1 || ptinpoly::pip2d(whichcell[[3]],cbind(i,j)) == 0)
      {
        toadd <- cbind(i,j,images[[frame]][j+1,i+1])
        intensities<-rbind(intensities,toadd)
        k = k+1
      }
      if(ptinpoly::pip2d(whichcell[[3]],cbind(i,j)) == -1)
      {
        toadd <- cbind(i,j,-1)
        intensities<-rbind(intensities,toadd)
        k = k+1
      }
    }
  }
  intensities<-na.omit(intensities)
  
  minx = min(intensities[,1])
  miny = min(intensities[,2])
  width = (max(intensities[,1])-min(intensities[,1]))+1
  height = (max(intensities[,2])-min(intensities[,2]))+1
  
  for(i in c(1:dim(intensities)[1]))
  {
    intensities[i,1] = intensities[i,1]-minx
    intensities[i,2] = intensities[i,2]-miny
  }
  
  indices<-vector()
  
  for(i in c(1:dim(intensities)[1]))
  {
    indices[i] = intensities[i,1] + (width * intensities[i,2])
  }
  
  intensities<-cbind(intensities,indices)
  intensities<-intensities[order(indices),]
  
  orderedintensities<-c(frame,width, height)
  
  for(i in c(1:dim(intensities)[1]))
  {
    orderedintensities[3+i] = intensities[i,3]
  }
  return(orderedintensities)
}

writeMiniImage<-function(cellID)
{
  print(cellID)
  justcell<-subset(roiData,lapply(roiData,'[[', 2)==cellID)
  
  miniImage<-list()
  k = 1
  
  for(j in c(1:length(justcell)))
  {
    cell = justcell[[j]][[2]]
    frame = justcell[[j]][[1]]
    image<-tiff
    miniImage[[k]] = storeInteriorIntensities(justcell[[j]], image, as.numeric(frame))
    k = k+1
    
  }
  return(miniImage)
}

multiReadTiff<-function(foldername)
{
  tiffFiles<-list.files(foldername, pattern="*tif$", full.name=TRUE)
  tiffList <- lapply(tiffFiles, tiff::readTIFF)
  for(i in c(1:length(tiffList)))
  {
    min = 0
    max = 1
    scale = 255*(max-min)
    tiffList[[i]]=(tiffList[[i]] - min)*scale
  }
  return(tiffList)
}

tiff<-multiReadTiff("/Users/laurawiggins/Desktop/05062019_B3/05062019_B3_3_imagedata")

miniImage<-writeMiniImage(cellID)

##miniImage stores the mini image file for cell cellID, rows in the same order as the boundary file
  
