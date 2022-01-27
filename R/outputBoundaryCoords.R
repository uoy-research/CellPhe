#' Extract a table of cell boundary coordinates from a zip file of ROIs
#'
#' @param zipFileOfROIs Name of zipped folder containing each ROI for boundary coordinate extraction, note that filenames must be entered in quotation marks and should exclude the .zip extension
#'
#' @return This function outputs a data table of extracted cell boundary coordinates and saves this as a global variable called \code{boundarycoords} where:
#' @return The value in \code{boundarycoords[1,1]} is the maximum boundary length of all cells within the data set
#' @return \code{boundarycoords[,1]} is a column of frame numbers 
#' @return \code{boundarycoords[,2]} is a column of cell IDs
#' @return \code{boundarycoords[,3]} is a column of boundary lengths
#' @return Remaining columns list x and y coordinates alternately until all boundary coordinates have been listed, in which case 0s are used to fill the remainder of the row
#'
#' @examples outputBoundaryCoords("myZipFolderOfROIs")
#' @export
outputBoundaryCoords<-function(zipFileOfROIs)
{
  zipFileOfROIs<-paste(zipFileOfROIs,".zip",sep='')
  roi<-RImageJROI::read.ijzip(zipFileOfROIs)
  maxboundarylength = 0
  
  for(i in c(1:length(roi)))
  {
    if(roi[[i]]$n > maxboundarylength) 
      maxboundarylength = roi[[i]]$n
    
  }
  
  boundarycoords<-matrix(ncol = (maxboundarylength*2)+3, nrow = length(roi)+1)
  boundarycoords[1,1] = maxboundarylength
  
  
  for(i in c(1:length(roi)))
  {
    print(i)
    boundarycoords[i+1,1] = strsplit(roi[[i]]$name, "-")[[1]][1]
    boundarycoords[i+1,2] = strsplit(roi[[i]]$name, "-")[[1]][2]
    boundarycoords[i+1,3] = (roi[[i]]$n)
    
    for(j in c(1:length(roi[[i]]$coords)))
    {
      
      boundarycoords[i+1,4+(j-1)] = as.list(t(roi[[i]]$coords))[[j]]
    }
  }
  boundarycoords[is.na(boundarycoords)]=0
  zipFileOfROIs = strsplit(zipFileOfROIs, "[.]")[[1]][1]
  zipFileOfROIs = paste(zipFileOfROIs,"-boundaryfile",sep='')
  assign(zipFileOfROIs, boundarycoords, envir = .GlobalEnv)
}