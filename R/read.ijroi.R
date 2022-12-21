##    https://gist.github.com/tractatus/9baea36c3c53b937dcc7f6aa06a1f082

##    readroi.R - Read ImageJ files in to R
##    Copyright (C) 2011 David C. Sterratt <david.c.sterratt@ed.ac.uk>

##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

read.ijroi <- function(file, verbose=FALSE) {
  ## Define internal helper functions
  getByte <- function(con) {
    pos <- seek(con) 
    n <- readBin(con, raw(0), 1, size=1)
    if (verbose)
      message(paste("Pos ", pos , ": Byte ", n, sep=""))
    return(as.integer(n))
  }

  getShort <- function(con) {
    pos <- seek(con) 
    n <- readBin(con, integer(0), 1, size=2, signed=TRUE, endian="big")
    if (n < -5000) {
      seek(con, -2, origin="current")
      n <- readBin(con, integer(0), 1, size=2, signed=FALSE, endian="big")
    }
    if (verbose)
      message(paste("Pos ", pos , ": Short ", n, sep=""))
    return(n)
  }
  
  getInt <- function(con)  {
    pos <- seek(con) 
    n <- readBin(con, integer(0), 1, size=4, signed=TRUE, endian="little")
    if (verbose)
      message(paste("Pos ", pos , ": Integer ", n, sep=""))
    return (n);
  }

  getFloat <- function(con)  {
    pos <- seek(con) 
    n <- readBin(con, double(0), 1, size=4, signed=TRUE, endian="big")
    if (verbose)
      message(paste("Pos ", pos , ": Float ", n, sep=""))
    return (n);
  }
  
  ## subtypes
  subtypes <- list(TEXT    = 1,
                   ARROW   = 2,
                   ELLIPSE = 3)

  ## options
  opts <- list(SPLINE_FIT    = as.raw(1),
               DOUBLE_HEADED = as.raw(2),
               OUTLINE       = as.raw(4))
  
  ## types
  types <- list(polygon  = 0,
                rect     = 1,
                oval     = 2,
                line     = 3,
                freeline = 4,
                polyline = 5,
                noRoi    = 6,
                freehand = 7,
                traced   = 8,
                angle    = 9,
                point    = 10)

  ## Main code
  name <- NULL
  if (!is.null(file)) {
    size <- file.info(file)$size
    if (!grepl(".roi$", file) && size>5242880)
      stop("This is not an ROI or file size>5MB)")
    name <- basename(file)
  }
  
  ## Open the connection
  con <- file(file, "rb")

  ## Test that it's the right kind of file
  if (getByte(con) != 73 || getByte(con) != 111) {  ## "Iout"
    stop("This is not an ImageJ ROI");
  }

  if (verbose)
    message("Reading format data")
  
  ## Create place to store data
  r <- list()

  ## Get the data. This all has to be in the order corresponding to the
  ## positions mentioned at the top of the file
  getShort(con)                         # Unused
  r$version <-      getShort(con)
  r$type <-         getByte(con)
  getByte(con)                          # Unused
  r$top <-          getShort(con)       # TOP
  r$left <-         getShort(con)       # LEFT
  r$bottom <-       getShort(con)       # Bottom
  r$right <-        getShort(con)       # RIGHT
  r$width <-    with(r, right-left)
  r$height <-   with(r, bottom-top)
  r$n <-            getShort(con)       # N_COORDINATES
  r$x1 <-           getFloat(con)     
  r$y1 <-           getFloat(con)     
  r$x2 <-           getFloat(con)     
  r$y2 <-           getFloat(con)
  r$strokeWidth <-  getShort(con)       # STROKE_WIDTH
  r$shapeRoiSize <- getInt(con)         # SHAPE_ROI_SIZE
  r$strokeColor <-  getInt(con)
  r$fillColor <-    getInt(con)
  r$subtype <-      getShort(con)       # SUBTYPE
  r$options < as.raw(getShort(con))     # OPTIONS
  if (r$type == "oval") {
    r$aspectRatio <- getFloat(con)      # ELLIPSE_ASPECT_RATIO
  } else {
    r$style <-      getByte(con)        # ARROW_STYLE
    r$headSize <-   getByte(con)        # ARROW_HEAD_SIZE    
    r$arcSize <-    getShort(con)       # ROUNDED_RECT_ARC_SIZE
  }
  r$position <-     getInt(con)         # POSITION
  getShort(con)                         # Unused
  getShort(con)                         # Unused

  if (verbose)
    message("Reading coordinate data")
  
  if (!is.null(name) && (grepl(".roi$", name)))
    r$name <- substring(name, 1, nchar(name) - 4)
    
  isComposite <- (r$shapeRoiSize >0);
  if (isComposite) {
    stop("Composite ROIs not supported")
    ## roi = getShapeRoi();
    ## if (version>=218) getStrokeWidthAndColor(roi);
    ##          roi.setPosition(position);
    ## return roi;
  }
  ## if (r$type %in% types[c("rect")]) {

  ## if (r$type %in% types[c("oval")]) {

  if (r$type %in% types["line"]) {
    if (r$subtype %in% types["ARROW"]) {
      r$doubleHeaded <- (r$options & opts$DOUBLE_HEADED) !=0
      r$outline <- (r$options & opts$OUTLINE) !=0
      ##                     if (style>=Arrow.FILLED && style<=Arrow.OPEN)
      ##                         ((Arrow)roi).setStyle(style);

      ##                     if (headSize>=0 && style<=30)
      ##                         ((Arrow)roi).setHeadSize(headSize);
      ##                 } else
      ##                     roi = new Line(x1, y1, x2, y2);
    }
  }

  ## Read in coordinates
  if (r$type %in% types[c("polygon", "freehand", "traced", "polyline", "freeline", "angle", "point")]) {
    r$coords <- matrix(NA, r$n, 2)
    if (r$n > 0) {
      for (i in 1:r$n) {
        r$coords[i, 1] <- getShort(con)
      }
      for (i in 1:r$n) {
        r$coords[i, 2] <- getShort(con)
      }
      r$coords[r$coords<0] <- 0
      r$coords[,1] <- r$coords[,1] + r$left
      r$coords[,2] <- r$coords[,2] + r$top
    }
  } 

  r$strType <- names(types)[which(types==r$type)]
  r$types <- types
  close(con)
  class(r) <- "IJROI"
  return(r)
} 
