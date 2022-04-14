/*******************************************************************************
 * FILE:        singlecell.cpp
 * AUTHOR:      Julie Wilson <julie.wilson@york.ac.uk>
 * AUTHOR:      Killian Murphy <killian.murphy@york.ac.uk>
 * DESCRIPTION: TODO
 ******************************************************************************/

/*******************************************************************************
 * HEADERS
 ******************************************************************************/
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Rcpp.h>
/*******************************************************************************
 * MACROS
 ******************************************************************************/
#define PI 3.14159265
/*******************************************************************************
 * TYPEDEFS
 ******************************************************************************/
typedef struct varnames {
  char var[100];
} NAMES;

typedef struct inputvars {
  double *frame;
  double *stats;
  double *vars;
} INVARS;

typedef struct bpix {
  int *xpix;
  int *ypix;
  int blength;
  int xlim1;
  int xlim2;
  int ylim1;
  int ylim2;
} BOUND;

typedef struct apix {
  int *xpix;
  int *ypix;
  int *intensity;
  int npix;
  int width;
  int height;
  int *image;
  int *mask;
  int lev0num;
  int lev1num;
  int lev2num;
  double *lev0Pix;
  double *lev1Pix;
  double *lev2Pix;
  double *cooc01;
  double *cooc12;
  double *cooc02;
} AREA;
/*******************************************************************************
 * FUNCTION DECLARATIONS
 ******************************************************************************/
int readLine(FILE *fp);
void cooccur(AREA *object, int cooccurrence_levels, int nframes,
             int *missingframe);
void minBox(BOUND *boundary, INVARS *input, int inputnum, int nframes,
            int *missingframe);
void varFromCentre(BOUND *boundary, INVARS *input, int maximum_boundary_length,
                   int inputnum, int nframes, int *missingframe);
void curvature(BOUND *boundary, INVARS *input, int inputnum, int nframes,
               int *missingframe);
void atob(BOUND *boundary, INVARS *input, int inputnum, int nframes,
          int *missingframe);
void polyClass(BOUND *boundary, INVARS *input, int inputnum, int nframes,
               int *missingframe);
void textureVariables(AREA *object, BOUND *boundary, INVARS *input,
                      int inputnum, int nframes, int *missingframe);
void firstOrderOriginal(AREA *object, INVARS *input, int maximum_cell_area,
                        int inputnum, int nframes, int *missingframe);
void cooccurVariables(AREA *object, int cooccurrence_levels, INVARS *input,
                      int inputnum, int nframes, int *missingframe);
void interpolate(INVARS *input, int numinput, int nframes, int *missingframe);
double timeSeriesVars(INVARS *input, int max_number_of_frames,
                      int number_of_wavelet_levels, int numinput, int nframes,
                      int *missingframe);
void getCoocMatrix(double *cooc, int cooccurrence_levels, double *bigimage,
                   int *bigmask, double *smallimage, int *smallmask, int w,
                   int h, int numlevs);
void waveTran(double *inputArray, int number_of_wavelet_levels, int arraynum,
              double *outputArray, int *detlength);
void haralick(double *cooc, int cooccurrence_levels, int num, double *hf);
void writedata(INVARS *input, int number_of_wavelet_levels, int numinput,
               int nstats, double areaoftraject, NAMES *vname,
               char *classlabel);
std::vector<std::string> create_column_labels(NAMES *vname,
                                              int number_of_variables);
//' Extract Features
//'
//' @param feature_table TODO
//' @param boundary_coordinates TODO
//' @param mini_image TODO
//' @param class_label TODO
//' @param max_number_of_frames TODO
//' @param maximum_boundary_length TODO
//' @param maximum_cell_area TODO
//' @param cooccurrence_levels TODO
//' @param number_of_wavelet_levels TODO
//'
//' @return TODO
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame
extract(Rcpp::NumericMatrix feature_table, Rcpp::List boundary_coordinates,
        Rcpp::List mini_image, const std::string class_label,
        const int max_number_of_frames, const int maximum_boundary_length,
        const int maximum_cell_area, const int cooccurrence_levels,
        const int number_of_wavelet_levels);

/*******************************************************************************
 * FUNCTION DEFINITIONS
 ******************************************************************************/
Rcpp::DataFrame
extract(Rcpp::NumericMatrix feature_table, Rcpp::List boundary_coordinates,
        Rcpp::List mini_image, const std::string class_label,
        const int max_number_of_frames, const int maximum_boundary_length,
        const int maximum_cell_area, const int cooccurrence_levels,
        const int number_of_wavelet_levels) {
  int j, k, ind, ix, iy, xpix, ypix, numvars, nframes, framenum;
  int tmp, dtmp, dtmp1, dtmp2, dtmp3, dtmp4;
  float X, Y, Volume, Thickness, Radius, Area, Sphericity, Vel1, Vel2;
  float Length, Width, Orientation, Mass, Displacement, Velocity, TrackLength;

  int startframe;
  double areaoftraject;

  int numinput = 49;
  int nstats = 3;

  FILE *fp = NULL;
  char ftfile[100];
  char bfile[100];
  char imfile[100];

  char ch;
  char classlabel[100];

  strcpy(classlabel, class_label.c_str());

  int *missingframe = NULL;
  missingframe = (int *)malloc(max_number_of_frames * sizeof(int));
  for (j = 0; j < max_number_of_frames; j++) {
    missingframe[j] = 1;
  }

  INVARS *input = NULL;
  input = (INVARS *)malloc(numinput * sizeof(INVARS));
  for (j = 0; j < numinput; j++) {
    input[j].frame = (double *)malloc(max_number_of_frames * sizeof(double));
    input[j].stats = (double *)malloc(nstats * sizeof(double));
    input[j].vars =
        (double *)malloc((number_of_wavelet_levels + 1) * 3 * sizeof(double));
  }

  NAMES *vname = NULL;
  vname = (NAMES *)malloc(numinput * sizeof(NAMES));

  BOUND *boundary = NULL;
  boundary = (BOUND *)malloc(max_number_of_frames * sizeof(BOUND));
  for (j = 0; j < max_number_of_frames; j++) {
    boundary[j].xpix = (int *)malloc(maximum_boundary_length * sizeof(int));
    boundary[j].ypix = (int *)malloc(maximum_boundary_length * sizeof(int));
    for (k = 0; k < maximum_boundary_length; k++) {
      boundary[j].xpix[k] = 0;
      boundary[j].ypix[k] = 0;
    }
  }

  AREA *object = NULL;
  object = (AREA *)malloc(max_number_of_frames * sizeof(AREA));
  for (j = 0; j < max_number_of_frames; j++) {
    object[j].xpix = (int *)malloc(maximum_cell_area * sizeof(int));
    object[j].ypix = (int *)malloc(maximum_cell_area * sizeof(int));
    object[j].intensity = (int *)malloc(maximum_cell_area * sizeof(int));
    object[j].image = (int *)malloc(maximum_cell_area * sizeof(int));
    object[j].mask = (int *)malloc(maximum_cell_area * sizeof(int));
    object[j].lev0Pix = (double *)malloc(maximum_cell_area * sizeof(double));
    object[j].lev1Pix = (double *)malloc(maximum_cell_area * sizeof(double));
    object[j].lev2Pix = (double *)malloc(maximum_cell_area * sizeof(double));
    object[j].cooc01 = (double *)malloc(cooccurrence_levels *
                                        cooccurrence_levels * sizeof(double));
    object[j].cooc12 = (double *)malloc(cooccurrence_levels *
                                        cooccurrence_levels * sizeof(double));
    object[j].cooc02 = (double *)malloc(cooccurrence_levels *
                                        cooccurrence_levels * sizeof(double));
  }

  startframe = -1;

  // TODO: iterate through input feature matrix here and copy variables
  for (int row_index = 0; row_index < feature_table.nrow(); ++row_index) {
    framenum = (feature_table(row_index, 0)) - 1;
    if (startframe == -1) {
      startframe = framenum;
    }
    framenum -= startframe;

    X = feature_table(row_index, 4);
    Y = feature_table(row_index, 5);
    xpix = feature_table(row_index, 6);
    ypix = feature_table(row_index, 7);
    Volume = feature_table(row_index, 8);
    Thickness = feature_table(row_index, 9);
    Radius = feature_table(row_index, 10);
    Area = feature_table(row_index, 11);
    Sphericity = feature_table(row_index, 12);
    Length = feature_table(row_index, 13);
    Width = feature_table(row_index, 14);
    Orientation = feature_table(row_index, 15);
    Mass = feature_table(row_index, 16);
    Displacement = feature_table(row_index, 17);
    Velocity = feature_table(row_index, 18);
    Vel1 = feature_table(row_index, 19);
    Vel2 = feature_table(row_index, 20);
    TrackLength = feature_table(row_index, 21);

    input[0].frame[framenum] = X;
    input[1].frame[framenum] = Y;
    input[2].frame[framenum] = Volume;
    strcpy(vname[2].var, "Vol\0");
    input[3].frame[framenum] = Radius;
    strcpy(vname[3].var, "Rad\0");
    input[4].frame[framenum] = Sphericity;
    strcpy(vname[4].var, "Sph\0");
    input[5].frame[framenum] = Length;
    strcpy(vname[5].var, "Len\0");
    input[6].frame[framenum] = Width;
    strcpy(vname[6].var, "Wid\0");
    input[7].frame[framenum] = Velocity;
    strcpy(vname[7].var, "Velocity\0");
    input[8].frame[framenum] = Displacement;
    strcpy(vname[8].var, "Dis\0");
    input[9].frame[framenum] = TrackLength;
    strcpy(vname[9].var, "Trac\0");
    input[10].frame[framenum] = 0.0;
    if (TrackLength > 0.0)
      input[10].frame[framenum] = Displacement / TrackLength;
    strcpy(vname[10].var, "D2T\0");
    missingframe[framenum] = 0;
  }

  numvars = 11;
  nframes = framenum + 1;

  /* read in boundary data */
  for (const Rcpp::IntegerVector coordinates : boundary_coordinates) {
    framenum = coordinates[0] - 1 - startframe;
    boundary[framenum].blength = coordinates[2];

    int stride = 3;
    for (k = 0; k < boundary[framenum].blength; ++k) {
      boundary[framenum].xpix[k] = coordinates[stride];
      boundary[framenum].ypix[k] = coordinates[stride + 1];
      stride += 2;
    }
  }

  /* read in mini image */
  for (const Rcpp::IntegerVector frame : mini_image) {
    ind = 0;
    framenum = frame[0] - 1 - startframe;
    object[framenum].npix = 0;
    object[framenum].width = frame[1];
    object[framenum].height = frame[2];

    //    TODO: check why this doesn't work
    //    for (k = 0; k < (frame.end() - (frame.begin() + 3)); ++k) {
    for (k = 0; k < (frame[1] * frame[2]); ++k) {
      object[framenum].mask[k] = 1;
      object[framenum].image[k] = frame[k + 3];

      if (object[framenum].image[k] == -1) {
        object[framenum].image[k] = 0;
        object[framenum].mask[k] = 0;
      } else {
        object[framenum].intensity[object[framenum].npix] = frame[k + 3];

        ix = ind % object[framenum].width;
        iy = (ind - ix) / object[framenum].width;

        object[framenum].xpix[object[framenum].npix] = ix;
        object[framenum].ypix[object[framenum].npix] = iy;

        object[framenum].npix++;
      }

      ind++;
    }
  }

  cooccur(object, cooccurrence_levels, nframes, missingframe);
  varFromCentre(boundary, input, maximum_boundary_length, numvars, nframes,
                missingframe);
  strcpy(vname[11].var, "VfC\0");
  numvars = numvars + 1;
  curvature(boundary, input, numvars, nframes, missingframe);
  strcpy(vname[12].var, "Cur\0");
  numvars = numvars + 1;
  for (k = 0; k < nframes; k++) {
    input[13].frame[k] = object[k].npix;
  }
  strcpy(vname[13].var, "Area\0");
  numvars = numvars + 1;
  atob(boundary, input, numvars, nframes, missingframe);
  strcpy(vname[14].var, "A2B\0");
  numvars = numvars + 1;
  minBox(boundary, input, numvars, nframes, missingframe);
  strcpy(vname[15].var, "Box\0");
  strcpy(vname[16].var, "Rect\0");
  numvars = numvars + 2;
  polyClass(boundary, input, numvars, nframes, missingframe);
  strcpy(vname[17].var, "poly1\0");
  strcpy(vname[18].var, "poly2\0");
  strcpy(vname[19].var, "poly3\0");
  strcpy(vname[20].var, "poly4\0");
  numvars = numvars + 4;
  textureVariables(object, boundary, input, numvars, nframes, missingframe);
  strcpy(vname[21].var, "IQ1\0");
  strcpy(vname[22].var, "IQ2\0");
  strcpy(vname[23].var, "IQ3\0");
  strcpy(vname[24].var, "IQ4\0");
  strcpy(vname[25].var, "IQ5\0");
  strcpy(vname[26].var, "IQ6\0");
  strcpy(vname[27].var, "IQ7\0");
  strcpy(vname[28].var, "IQ8\0");
  strcpy(vname[29].var, "IQ9\0");
  numvars = numvars + 9;
  firstOrderOriginal(object, input, maximum_cell_area, numvars, nframes,
                     missingframe);
  strcpy(vname[30].var, "FOmean\0");
  strcpy(vname[31].var, "FOsd\0");
  strcpy(vname[32].var, "FOskew\0");
  numvars = numvars + 3;
  cooccurVariables(object, cooccurrence_levels, input, numvars, nframes,
                   missingframe);
  strcpy(vname[33].var, "Hf1_01\0");
  strcpy(vname[34].var, "Hf2_01\0");
  strcpy(vname[35].var, "Hf3_01\0");
  strcpy(vname[36].var, "Hf4_01\0");
  strcpy(vname[37].var, "Hf5_01\0");
  strcpy(vname[38].var, "Hf1_12\0");
  strcpy(vname[39].var, "Hf2_12\0");
  strcpy(vname[40].var, "Hf3_12\0");
  strcpy(vname[41].var, "Hf4_12\0");
  strcpy(vname[42].var, "Hf5_12\0");
  strcpy(vname[43].var, "Hf1_02\0");
  strcpy(vname[44].var, "Hf2_02\0");
  strcpy(vname[45].var, "Hf3_02\0");
  strcpy(vname[46].var, "Hf4_02\0");
  strcpy(vname[47].var, "Hf5_02\0");
  numvars = numvars + 15;
  /* -1.0 for density is temporary!!!*/
  for (k = 0; k < nframes; k++) {
    input[numvars].frame[k] = -1.0;
  }
  strcpy(vname[48].var, "Den\0");

  interpolate(input, numinput, nframes, missingframe);
  areaoftraject =
      timeSeriesVars(input, max_number_of_frames, number_of_wavelet_levels,
                     numinput, nframes, missingframe);

  std::vector<std::string> column_labels =
      create_column_labels(vname, (numinput - 2));

  std::vector<double> wavelet_values;
  std::vector<double> statistics_values;

  for (int index = 2; index < numinput; ++index) {
    std::vector<double> temp_wavelet_values(
        input[index].vars,
        input[index].vars + ((number_of_wavelet_levels + 1) * 3));
    temp_wavelet_values.erase(temp_wavelet_values.begin() + 3,
                              temp_wavelet_values.begin() + 6);
    wavelet_values.insert(wavelet_values.end(), temp_wavelet_values.begin(),
                          temp_wavelet_values.end());

    std::vector<double> temp_stats_values(input[index].stats,
                                          input[index].stats + 3);
    statistics_values.insert(statistics_values.end(), temp_stats_values.begin(),
                             temp_stats_values.end());
  }

  std::vector<double> combined_values;
  wavelet_values.insert(wavelet_values.end(), statistics_values.begin(),
                        statistics_values.end());
  combined_values.insert(combined_values.end(), wavelet_values.begin(),
                         wavelet_values.end());

  Rcpp::DataFrame results;

  results[column_labels[0]] = Rcpp::CharacterVector::create(classlabel);
  results[column_labels[1]] = Rcpp::NumericVector::create(areaoftraject);

  for (int index = 2; index < column_labels.size(); ++index) {
    results[column_labels[index]] =
        Rcpp::NumericVector::create(combined_values[index - 2]);
  }

  for (int index = 0; index < (number_of_wavelet_levels); ++index) {
    results.erase(results.findName("Trac_l" + std::to_string(index) + "_des"));
  }

  free(missingframe);

  for (j = 0; j < numinput; j++) {
    free(input[j].frame);
    free(input[j].stats);
    free(input[j].vars);
  }
  free(input);

  for (j = 0; j < max_number_of_frames; j++) {
    free(boundary[j].xpix);
    free(boundary[j].ypix);
  }
  free(boundary);

  for (j = 0; j < max_number_of_frames; j++) {
    free(object[j].xpix);
    free(object[j].ypix);
    free(object[j].intensity);
    free(object[j].image);
    free(object[j].mask);
    free(object[j].lev0Pix);
    free(object[j].lev1Pix);
    free(object[j].lev2Pix);
    free(object[j].cooc01);
    free(object[j].cooc12);
    free(object[j].cooc02);
  }
  free(object);

  free(vname);

  return results;
}

/*------------------------------------------------------------------------------
 * Procedure:   create_column_labels
 * Description: Create a vector of column labels for returning a data.frame of
 *              results to R.
 *
 *              Subtraction of 2 from numinputs is to match the way in which the
 *              vname data is currently stored. Addition of 1 to
 *              number_of_wavelet_levels is to account for level 0 of the
 *              transform. Subtraction of 4 from the total is to account for
 *              removal of 'descent' of TrackLength.
 *----------------------------------------------------------------------------*/
std::vector<std::string> create_column_labels(NAMES *vname,
                                              int number_of_variables) {
  std::vector<std::string> variable_names;
  std::vector<std::string> column_labels = {"class", "trajarea"};
  std::vector<std::string> wavelet_columns;
  std::vector<std::string> statistics_columns;

  for (int index = 2; index < (number_of_variables + 2); ++index) {
    variable_names.push_back(std::string(vname[index].var));
  }

  std::vector<std::string> wavelet_suffixes = {"_asc", "_des", "_max"};
  std::vector<std::string> wavelet_levels = {"_l0", "_l1", "_l2", "_l3"};

  for (const auto &variable_name : variable_names) {
    for (const auto &wavelet_level : wavelet_levels) {
      for (const auto &wavelet_suffix : wavelet_suffixes) {
        //        if (variable_name == "Trac" && wavelet_suffix == "_des")
        //          continue;
        wavelet_columns.push_back(variable_name + wavelet_level +
                                  wavelet_suffix);
      }
    }
  }

  std::vector<std::string> statistic_suffixes = {"_mean", "_std", "_skew"};

  for (const auto &variable_name : variable_names) {
    for (const auto &statistic_suffix : statistic_suffixes) {
      statistics_columns.push_back(variable_name + statistic_suffix);
    }
  }

  wavelet_columns.insert(wavelet_columns.end(), statistics_columns.begin(),
                         statistics_columns.end());
  column_labels.insert(column_labels.end(), wavelet_columns.begin(),
                       wavelet_columns.end());

  return column_labels;
}
/*------------------------------------------------------------------------------
 * Procedure: readLine
 * Description: reads past a line
 *----------------------------------------------------------------------------*/
int readLine(FILE *fp) {
  int ch;

  ch = fgetc(fp);

  while (ch != 10) {
    ch = fgetc(fp);
  }

  return ch;
}

/*------------------------------------------------------------------------------
 * Procedure:  minBox
 * Description:    calculates minimal rectangular box that the cell will fit in
 *----------------------------------------------------------------------------*/
void minBox(BOUND *boundary, INVARS *input, int inputnum, int nframes,
            int *missingframe) {
  int j, k, l;
  int x1, x2, y1, y2;
  double dist, maxdist = 0.0;
  int keepx1 = 0, keepx2 = 0, keepy1 = 0, keepy2 = 0;
  double xlength, ylength;
  double alpha, rotx, roty, maxx, minx, maxy, miny, max;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      /* find maximum distance between two boundary points */
      maxdist = 0.0;
      for (j = 0; j < boundary[k].blength; j++) {

        for (l = j + 1; l < boundary[k].blength; l++) {
          x1 = boundary[k].xpix[j];
          y1 = boundary[k].ypix[j];
          x2 = boundary[k].xpix[l];
          y2 = boundary[k].ypix[l];
          dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
          if (dist > maxdist) {
            maxdist = dist;
            keepx1 = x1;
            keepx2 = x2;
            keepy1 = y1;
            keepy2 = y2;
          }
        }
      }
      alpha = atan2((double)(keepy1 - keepy2), (double)(keepx1 - keepx2));

      /* rotating points by -alpha makes keepx1 -keepx2 lie along x-axis */
      /* get rotated y values and find max-min for ylength */
      maxx = 0.0;
      minx = 5000.0;
      maxy = 0.0;
      miny = 5000.0;
      /* double rotx */
      for (j = 0; j < boundary[k].blength; j++) {
        roty = keepy1 - sin(alpha) * (float)(boundary[k].xpix[j] - keepx1) +
               cos(alpha) * (float)(boundary[k].ypix[j] - keepy1);
        rotx = keepx1 + cos(alpha) * (float)(boundary[k].xpix[j] - keepx1) +
               sin(alpha) * (float)(boundary[k].ypix[j] - keepy1);
        if (roty < miny)
          miny = roty;
        if (roty > maxy)
          maxy = roty;
        if (rotx < minx)
          minx = rotx;
        if (rotx > maxx)
          maxx = rotx;
        /*   if (k == 0) printf("%d %f %f\n",j, rotx, roty);*/
      }
      xlength = maxx - minx;
      ylength = maxy - miny;
      input[inputnum].frame[k] =
          (float)(xlength * ylength) / input[13].frame[k];
      max = xlength;
      if (ylength > max)
        max = ylength;
      input[inputnum + 1].frame[k] = max / (xlength + ylength);
    } else {
      input[inputnum].frame[k] = -1.0;
      input[inputnum + 1].frame[k] = -1.0;
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  getVar
 * Description: calcultes the variance/mean of an array
 *----------------------------------------------------------------------------*/
double getVar(double *array, int num) {
  int i;
  double var, sum1, sum2, fnum;

  sum1 = 0.0;
  sum2 = 0.0;
  fnum = (float)num;
  for (i = 0; i < num; i++) {
    sum1 += array[i];
    sum2 += array[i] * array[i];
  }
  sum1 /= fnum;
  sum2 /= fnum;
  var = sum2 - (sum1 * sum1);
  var = var / sum1;
  return (var);
}

/*------------------------------------------------------------------------------
 * Procedure:  varFromCentre
 * Description: finds variance of the didtance from the boundary to the centre
 *----------------------------------------------------------------------------*/
void varFromCentre(BOUND *boundary, INVARS *input, int maximum_boundary_length,
                   int inputnum, int nframes, int *missingframe) {
  int j, k, x, y, icx, icy;

  double *dist = NULL;
  dist = (double *)malloc(maximum_boundary_length * sizeof(double));

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      for (j = 0; j < boundary[k].blength; j++) {
        x = boundary[k].xpix[j];
        y = boundary[k].ypix[j];
        icx = input[0].frame[k];
        icy = input[1].frame[k];
        dist[j] = sqrt((x - icx) * (x - icx) + (y - icy) * (y - icy));
      }
      input[inputnum].frame[k] = getVar(dist, boundary[k].blength);
    } else
      input[inputnum].frame[k] = -1.0;
  }
  free(dist);
}

/*------------------------------------------------------------------------------
 * Procedure:  curvature
 * Description:    calculates the curvature of the boundary
 *----------------------------------------------------------------------------*/
void curvature(BOUND *boundary, INVARS *input, int inputnum, int nframes,
               int *missingframe) {
  int j, k, x, y, n, xmgap, ymgap, xpgap, ypgap;
  int i1, i2, i3;
  double d1, d2, d3;
  double curv, newval;
  int gap = 4;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      n = boundary[k].blength;
      curv = 0.0;
      for (j = 0; j < n; j++) {
        x = boundary[k].xpix[j];
        y = boundary[k].ypix[j];
        xmgap = boundary[k].xpix[(j - gap + n) % n];
        ymgap = boundary[k].ypix[(j - gap + n) % n];
        xpgap = boundary[k].xpix[(j + gap + n) % n];
        ypgap = boundary[k].ypix[(j + gap + n) % n];
        i1 = (x - xmgap) * (x - xmgap) + (y - ymgap) * (y - ymgap);
        i2 = (x - xpgap) * (x - xpgap) + (y - ypgap) * (y - ypgap);
        i3 = (xmgap - xpgap) * (xmgap - xpgap) +
             (ymgap - ypgap) * (ymgap - ypgap);
        d1 = sqrt((float)i1);
        d2 = sqrt((float)i2);
        d3 = sqrt((float)i3);
        newval = d1 + d2 - d3;
        curv = curv + newval;
      }
      input[inputnum].frame[k] = curv / (float)n;
    } else
      input[inputnum].frame[k] = -1.0;
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  atob
 * Description:    calculates the ratio of the area to the boundary squared
 *----------------------------------------------------------------------------*/
void atob(BOUND *boundary, INVARS *input, int inputnum, int nframes,
          int *missingframe) {
  int k;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      input[inputnum].frame[k] =
          input[inputnum - 1].frame[k] /
          (float)(boundary[k].blength * boundary[k].blength);
    } else
      input[inputnum].frame[k] = -1.0;
  }
}

/*------------------------------------------------------------------------------
 * Procedure: pointttolinedist
 * Description:  calculates the distance between the point (x0, y0) and the line
 *----------------------------------------------------------------------------*/
double pointttolinedist(int x0, int y0, int x1, int y1, int x2, int y2) {
  double fnum, fden, dist;
  int numer, denom;

  numer = (y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1;
  denom = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
  fnum = fabs((float)numer);
  fden = sqrt((float)denom);
  dist = fnum / fden;

  return (dist);
}

/*------------------------------------------------------------------------------
 * Procedure: polygon
 * Description: extracts a set of points from the cell boundary (xArray and
 *              yArray) that form a polygon where the sides are within thresh
 *              of the boundary
 *----------------------------------------------------------------------------*/
void polygon(int *xArray, int *yArray, int num, double thresh, int *xPoints,
             int *yPoints, int *nPoints) {
  int k, i, xs, ys, x, y, numpoints, indkeep, n;
  int x0, y0, x1, y1, x2, y2;
  int alldone;
  double dist, distMax;

  xs = xArray[0];
  ys = yArray[0];
  distMax = 0.0;
  indkeep = 0;
  for (k = 0; k < num; k++) {
    x = xArray[k];
    y = yArray[k];
    dist = sqrt((double)((x - xs) * (x - xs) + (y - ys) * (y - ys)));
    if (dist > distMax) {
      distMax = dist;
      indkeep = k;
    }
  }

  int *pointArray = NULL;
  pointArray = (int *)malloc(num * sizeof(int));
  int *tempArray = NULL;
  tempArray = (int *)malloc(num * sizeof(int));
  int *done = NULL;
  done = (int *)malloc(num * sizeof(int));

  numpoints = 2;
  pointArray[0] = 0;
  pointArray[1] = indkeep;

  alldone = 0;
  while (alldone == 0) {
    alldone = 1;
    n = 0;
    tempArray[n] = 0;
    n++;
    for (k = 1; k < numpoints; k++) {
      x1 = xArray[pointArray[k - 1]];
      y1 = yArray[pointArray[k - 1]];
      x2 = xArray[pointArray[k]];
      y2 = yArray[pointArray[k]];
      distMax = 0.0;
      for (i = pointArray[k - 1] + 1; i < pointArray[k] - 1; i++) {
        x0 = xArray[i];
        y0 = yArray[i];
        dist = pointttolinedist(x0, y0, x1, y1, x2, y2);
        if (dist > distMax) {
          distMax = dist;
          indkeep = i;
        }
      }
      if (distMax > thresh) {
        tempArray[n] = indkeep;
        n++;
        alldone = 0;
      }
      tempArray[n] = pointArray[k];
      n++;
    }
    done[numpoints] = 0;
    x1 = xArray[pointArray[numpoints - 1]];
    y1 = yArray[pointArray[numpoints - 1]];
    x2 = xArray[pointArray[0]];
    y2 = yArray[pointArray[0]];
    distMax = 0.0;
    for (i = pointArray[numpoints - 1] + 1; i < num; i++) {
      x0 = xArray[i];
      y0 = yArray[i];
      dist = pointttolinedist(x0, y0, x1, y1, x2, y2);
      if (dist > distMax) {
        distMax = dist;
        indkeep = i;
      }
    }
    if (distMax > thresh) {
      tempArray[n] = indkeep;
      n++;
      alldone = 0;
    }
    numpoints = n;
    for (i = 0; i < numpoints; i++) {
      pointArray[i] = tempArray[i];
    }
  }

  for (i = 0; i < numpoints; i++) {
    xPoints[i] = xArray[pointArray[i]];
    yPoints[i] = yArray[pointArray[i]];
  }
  xPoints[i] = xPoints[0];
  yPoints[i] = yPoints[0];

  (*nPoints) = numpoints;

  free(pointArray);
  free(tempArray);
  free(done);
}

/*------------------------------------------------------------------------------
 * Procedure: polyclass
 * Description:  calculates variables from the polygonal estimates of the cells
 *----------------------------------------------------------------------------*/
void polyClass(BOUND *boundary, INVARS *input, int inputnum, int nframes,
               int *missingframe) {
  int j, k, x0, y0, x1, y1, x2, y2, num;
  double dist, distMax;
  double angle, angleMin;
  double meana, vara, meand, vard;
  double Asq, Bsq, Csq, a, b;
  double thresh = 2.0;
  int nPoints;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      num = boundary[k].blength;
      int *xArray = NULL;
      xArray = (int *)malloc(num * sizeof(int));
      int *yArray = NULL;
      yArray = (int *)malloc(num * sizeof(int));
      for (j = 0; j < num; j++) {
        xArray[j] = boundary[k].xpix[j];
        yArray[j] = boundary[k].ypix[j];
      }
      int *xPoints = NULL;
      xPoints = (int *)malloc((num + 2) * sizeof(int));
      int *yPoints = NULL;
      yPoints = (int *)malloc((num + 2) * sizeof(int));
      polygon(xArray, yArray, num, thresh, xPoints, yPoints, &nPoints);

      distMax = 0.0;
      angleMin = 2.0 * PI;
      meana = 0.0;
      vara = 0.0;
      meand = 0.0;
      vard = 0.0;
      xPoints[nPoints] = xPoints[0];
      yPoints[nPoints] = yPoints[0];
      xPoints[nPoints + 1] = xPoints[1];
      yPoints[nPoints + 1] = yPoints[2];
      for (j = 1; j < nPoints + 1; j++) {
        x0 = xPoints[j];
        y0 = yPoints[j];
        x1 = xPoints[j - 1];
        y1 = yPoints[j - 1];
        x2 = xPoints[j + 1];
        y2 = yPoints[j + 1];
        Asq = (double)((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
        Bsq = (double)((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));
        Csq = (double)((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        a = sqrt(Asq);
        b = sqrt(Bsq);
        dist = a;
        if (dist > distMax)
          distMax = dist;
        meand = meand + dist;
        vard = vard + dist * dist;
        if (fabs((Asq + Bsq - Csq) / (2.0 * a * b) - 1.0) > 0.001) {
          angle = acos((Asq + Bsq - Csq) / (2.0 * a * b));
          if (angle < angleMin)
            angleMin = angle;
          meana = meana + angle;
          vara = vara + angle * angle;
        }
      }
      meand /= (float)(nPoints + 1);
      vard /= (float)(nPoints + 1);
      meana /= (float)(nPoints + 1);
      vara /= (float)(nPoints + 1);
      input[inputnum].frame[k] = distMax;
      input[inputnum + 1].frame[k] = angleMin;
      input[inputnum + 2].frame[k] = vara - (meana * meana);
      input[inputnum + 3].frame[k] = vard - (meand * meand);
      free(xArray);
      free(yArray);
      free(xPoints);
      free(yPoints);
    } else {
      input[inputnum].frame[k] = -1;
      input[inputnum + 1].frame[k] = -1;
      input[inputnum + 2].frame[k] = -1;
      input[inputnum + 3].frame[k] = -1;
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  fSort
 * Description:    Sorts an array of integers into increasing order
 *----------------------------------------------------------------------------*/
void intSort(int *array, int nRef) {
  int nInt, ii, i, j, nf;
  int tf;
  nInt = nRef;
  if (nInt == 1)
    return;

  do {
    nInt = nInt / 2;
    if ((nInt % 2) == 0) {
      nInt = nInt - 1;
    }
    for (ii = 0; ii < nRef - nInt; ii++) {
      i = ii;
      j = i + nInt;
      if (array[i] > array[j]) {
        tf = array[j];
        do {
          array[j] = array[i];
          j = i;
          i = i - nInt;
        } while ((i >= 0) && (array[i] > tf));
        array[j] = tf;
      }
    }
  } while (nInt > 1);
}

/*------------------------------------------------------------------------------
 * Procedure:  textureVariables
 * Description:
 *----------------------------------------------------------------------------*/
void textureVariables(AREA *object, BOUND *boundary, INVARS *input,
                      int inputnum, int nframes, int *missingframe) {
  int width, height, c, k, a, term1, term2;
  int b, xp, xc, yp, yc, Indicator, floorh, num, inputvar;
  double maxx, maxy, r, d, h, q, polyarea;
  float quantile;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      inputvar = inputnum;
      term1 = 0;
      term2 = 0;
      Indicator = 0;
      d = 0;

      width = object[k].width;
      height = object[k].height;
      maxx = (float)(width - 1);
      maxy = (float)(height - 1);
      if (maxx < maxy)
        r = maxx * 0.25;
      else
        r = maxy * 0.25;

      int *pointstoevaluatex = NULL;
      pointstoevaluatex = (int *)malloc(object[k].npix * sizeof(int));
      int *pointstoevaluatey = NULL;
      pointstoevaluatey = (int *)malloc(object[k].npix * sizeof(int));

      int *sortintensities = NULL;
      sortintensities = (int *)malloc(object[k].npix * sizeof(int));

      for (c = 0; c < object[k].npix; c++) {
        sortintensities[c] = object[k].intensity[c];
      }
      intSort(sortintensities, object[k].npix);

      for (q = 0.1; q < 0.95; q += 0.1) {
        Indicator = 0;
        h = ((object[k].npix) - 1) * q + 1;
        floorh = floor(h);
        quantile = sortintensities[floorh - 1] +
                   (h - floorh) *
                       (sortintensities[floorh] - sortintensities[floorh - 1]);

        num = 0;
        for (c = 0; c < object[k].npix; c++) {
          if (object[k].intensity[c] >= quantile) {
            pointstoevaluatex[num] = object[k].xpix[c];
            pointstoevaluatey[num] = object[k].ypix[c];
            num++;
          }
        }
        term1 = 0;
        term2 = 0;
        polyarea = 0;
        for (a = 0; a < (boundary[k].blength - 1); a++) {
          term1 = term1 + (boundary[k].xpix[a] * boundary[k].ypix[a + 1]);
          term2 = term2 + (boundary[k].xpix[a + 1] * boundary[k].ypix[a]);
        }
        term1 = term1 + (boundary[k].xpix[boundary[k].blength - 1] *
                         boundary[k].ypix[0]);
        term2 = term2 + (boundary[k].xpix[0] *
                         boundary[k].ypix[boundary[k].blength - 1]);
        polyarea =
            0.5 * abs(term1 + (boundary[k].xpix[0] * boundary[k].ypix[0]) -
                      term2 - (boundary[k].xpix[0] * boundary[k].ypix[0]));

        for (a = 0; a < num; a++) {
          for (b = 0; b < num; b++) {
            if (a != b) {
              xp = pointstoevaluatex[b];
              yp = pointstoevaluatey[b];
              xc = pointstoevaluatex[a];
              yc = pointstoevaluatey[a];
              d = sqrt(((xp - xc) * (xp - xc)) + (yp - yc) * (yp - yc));
              if (d < r) {
                Indicator++;
              }
            }
          }
        }
        double Kemp = (polyarea / (num * (num - 1))) * Indicator;
        double Ktheo = PI * r * r;
        input[inputvar].frame[k] = Kemp - Ktheo;
        inputvar = inputvar + 1;
      }

      free(sortintensities);
      free(pointstoevaluatex);
      free(pointstoevaluatey);
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure: FOstats
 * Description: calculates mean, stdev and skewness from an array of values
 *----------------------------------------------------------------------------*/
void FOstats(double *array, int num, double *stats) {
  int j;
  double m1, m2, m3, en, enm1, enm2;

  en = (float)(num);
  enm1 = en - 1.0;
  enm2 = en - 2.0;

  stats[0] = 0.0;
  stats[1] = 0.0;
  stats[2] = 0.0;
  for (j = 0; j < num; j++) {
    stats[0] += array[j];
  }
  /* stat 0 is mean */
  stats[0] /= en;

  /* get moments */
  m1 = stats[0];
  m2 = 0.0;
  m3 = 0.0;
  for (j = 0; j < num; j++) {
    m2 += (array[j] - m1) * (array[j] - m1);
    m3 += (array[j] - m1) * (array[j] - m1) * (array[j] - m1);
  }
  /* stat 1 is standard deviation */
  stats[1] = sqrt(m2 / enm1);
  /* stat 2 is skewness */
  if ((m2 == 0.0) || (enm2 <= 0.0))
    stats[2] = 0.0;
  else
    stats[2] = (m3 * en * sqrt(enm1)) / (enm2 * m2 * sqrt(m2));
}

/*------------------------------------------------------------------------------
 * Procedure:  firstOrderOriginal
 * Description: Extracts first order features from pixel intensity histogram of
 *              each cell in each frame
 *----------------------------------------------------------------------------*/
void firstOrderOriginal(AREA *object, INVARS *input, int maximum_cell_area,
                        int inputnum, int nframes, int *missingframe) {
  int c, k;

  double *stats = NULL;
  stats = (double *)malloc(3 * sizeof(double));

  double *array = NULL;
  array = (double *)malloc(maximum_cell_area * sizeof(double));

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      for (c = 0; c < object[k].npix; c++) {
        array[c] = (float)object[k].intensity[c];
      }
      FOstats(array, object[k].npix, stats);

      input[inputnum].frame[k] = stats[0];
      input[inputnum + 1].frame[k] = stats[1];
      input[inputnum + 2].frame[k] = stats[2];
    } else {
      input[inputnum].frame[k] = -1.0;
      input[inputnum + 1].frame[k] = -1.0;
      input[inputnum + 2].frame[k] = -1.0;
    }
  }

  free(stats);
  free(array);
}

/*------------------------------------------------------------------------------
 * Procedure:  reScale
 * Description: real-valued input values are re-scaled to have a minimum equal
 *              to 0 and a maximum equal to 255
 *----------------------------------------------------------------------------*/
void reScale(int num, double *inputArray, int *mask) {
  int ind;
  double max, min, scale;
  min = 1000.0;
  max = -1000.0;

  for (ind = 0; ind < num; ind++) {
    if ((mask[ind] == 1) && (inputArray[ind] < min))
      min = inputArray[ind];
    if (inputArray[ind] > max)
      max = inputArray[ind];
  }

  if (max != min) {
    scale = 255.0 / (max - min);
    for (ind = 0; ind < num; ind++) {
      if (mask[ind] == 1) {
        inputArray[ind] = (inputArray[ind] - min) * scale;
      }
    }
  }

  for (ind = 0; ind < num; ind++) {
    if ((mask[ind] == 1) && (inputArray[ind] < min))
      min = inputArray[ind];
    if (inputArray[ind] > max)
      max = inputArray[ind];
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  getCoocMatrix
 * Description: creates co-occurrence matrices images on two levels
 *----------------------------------------------------------------------------*/
void getCoocMatrix(double *cooc, int cooccurrence_levels, double *bigimage,
                   int *bigmask, double *smallimage, int *smallmask, int w,
                   int h, int numlevs) {
  int i, j, ls, lb, ind, inds, indb, x, y, nc, twos;
  float n;

  nc = cooccurrence_levels;
  n = (float)nc;
  twos = pow(2, numlevs);
  /* rescale the two images */
  reScale(w * h, smallimage, smallmask);
  reScale(twos * twos * w * h, bigimage, bigmask);

  for (j = 0; j < cooccurrence_levels * cooccurrence_levels; j++) {
    cooc[j] = 0;
  }

  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      inds = i + j * w;
      for (x = 0; x < twos; x++) {
        for (y = 0; y < twos; y++) {
          indb = twos * i + x + twos * w * (twos * j + y);
          if (bigmask[indb] != 0) {
            ls = (int)(n * smallimage[inds] / 256.0);
            lb = (int)(n * bigimage[indb] / 256.0);
            ind = lb + ls * cooccurrence_levels;
            cooc[ind]++;
          }
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure: code for Daubechies' 2-coefficient (Haar) wavelet from Numerical
 *            Recipes
 *----------------------------------------------------------------------------*/
void daub2(double *a, int n, int isign) {
  double *wksp = NULL;
  int nh, i, j;
  double D0 = 0.70710678;
  double D1 = 0.70710678;

  wksp = (double *)malloc(n * sizeof(double));

  if (n < 4)
    return;
  nh = n / 2;
  if (isign >= 0) {
    i = 0;
    for (j = 0; j < n - 1; j += 2) {
      wksp[i] = D0 * a[j] + D1 * a[j + 1];
      wksp[i + nh] = D1 * a[j] - D0 * a[j + 1];
      i++;
    }
  } else {
    j = 0;
    for (i = 0; i < nh - 1; i++) {
      wksp[j] = D0 * a[i] + D1 * a[i + nh];
      wksp[j + 1] = D1 * a[i] - D0 * a[i + nh];
      j = j + 2;
    }
  }
  for (i = 0; i < n; i++) {
    a[i] = wksp[i];
  }
  free(wksp);
}

/*------------------------------------------------------------------------------
 * Procedure:    wavelet transform
 * Description:  performs 3-level wavelet transform
 *----------------------------------------------------------------------------*/
void waveTran(double *inputArray, int number_of_wavelet_levels, int arraynum,
              double *outputArray, int *detlength) {
  int x, k, length;
  double root2 = sqrt(2.0);
  int ww;

  detlength[0] = 0;
  length = arraynum;
  ww = 0;
  if (length % 2 != 0)
    ww = 1;

  double *xVector = NULL;
  xVector = (double *)malloc((length + ww) * sizeof(double));
  for (x = 0; x < length; x++) {
    xVector[x] = inputArray[x];
  }

  for (k = 0; k < number_of_wavelet_levels; k++) {
    ww = 0;
    if (length % 2 != 0)
      ww = 1;

    /* do transform */
    if (ww == 1)
      xVector[length] = xVector[length - 1];
    length = length + ww;

    daub2(xVector, length, 1);

    /* rescale for next level */
    for (x = 0; x < length - ww; x++) {
      xVector[x] /= root2;
    }

    length = length / 2;
    /* output detail coefficients */
    for (x = 0; x < length - ww; x++) {
      outputArray[x + detlength[k]] = inputArray[x + length];
    }
    detlength[k + 1] = detlength[k] + length - ww;
  }

  free(xVector);
}

/*------------------------------------------------------------------------------
 * Procedure:	1-level 2D wavelet transform
 * Description:    does 1-level x and y wavelet transform
 *----------------------------------------------------------------------------*/
void waveTran2D(int width, int height, double *inputImage, double *levImage) {
  int x, y, ind, yind, jnd;
  int xLength = width;
  int yLength = height;

  double *xVector = NULL;
  double *yVector = NULL;
  xVector = (double *)malloc(width * sizeof(double));
  yVector = (double *)malloc(height * sizeof(double));
  double *tempImage = NULL;
  tempImage = (double *)malloc(width * height * sizeof(double));

  /* do one level x transform */
  for (y = 0; y < yLength; y++) {
    yind = width * y;
    for (x = 0; x < xLength; x++) {
      ind = x + yind;
      xVector[x] = inputImage[ind];
    }
    daub2(xVector, xLength, 1);
    for (x = 0; x < xLength; x++) {
      ind = x + yind;
      tempImage[ind] = xVector[x] / sqrt(2.0);
    }
  }
  /* do one level y transform */
  for (x = 0; x < xLength; x++) {
    for (y = 0; y < yLength; y++) {
      ind = x + width * y;
      yVector[y] = tempImage[ind];
    }
    daub2(yVector, yLength, 1);
    for (y = 0; y < yLength; y++) {
      ind = x + width * y;
      tempImage[ind] = yVector[y] / sqrt(2.0);
    }
  }

  xLength = xLength / 2;
  yLength = yLength / 2;
  /* output smoothed image */
  for (x = 0; x < xLength; x++) {
    for (y = 0; y < yLength; y++) {
      ind = x + width * y;
      jnd = x + (width / 2) * y;
      levImage[jnd] = tempImage[ind];
    }
  }
  free(xVector);
  free(yVector);
  free(tempImage);
}

/*------------------------------------------------------------------------------
 * Procedure:  shrinkmask
 * Description: produces a mask at a particular level from the one at the
 *              previous level new mask is only 0 if all 4 pixels in previous
 *              level are 0
 *----------------------------------------------------------------------------*/
void shrinkmask(int width, int height, int *mask, int *newmask) {
  int i, j, inds, x, y, indb, n;

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      inds = i + j * width;
      newmask[inds] = 1;
      n = 0;
      for (x = 0; x < 2; x++) {
        for (y = 0; y < 2; y++) {
          indb = 2 * i + x + 2 * width * (2 * j + y);
          if (mask[indb] == 0)
            n++;
        }
      }
      if (n == 4)
        newmask[inds] = 0;
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  coocur
 * Description:
 *----------------------------------------------------------------------------*/
void cooccur(AREA *object, int cooccurrence_levels, int nframes,
             int *missingframe) {
  int j, k, x, y, ind, jnd;
  int width, height, newwidth, newheight;

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      width = object[k].width;
      height = object[k].height;
      /* check the image size is correct for a 2-level transform */
      newwidth = (width / 4) * 4;
      if (width - newwidth > 0)
        newwidth += 4;
      newheight = (height / 4) * 4;
      if (height - newheight > 0)
        newheight += 4;

      double *newimage = NULL;
      newimage = (double *)malloc(newwidth * newheight * sizeof(double));
      int *newmask = NULL;
      newmask = (int *)malloc(newwidth * newheight * sizeof(int));
      double *lev1image = NULL;
      lev1image =
          (double *)malloc((newwidth / 2) * (newheight / 2) * sizeof(double));
      int *lev1mask = NULL;
      lev1mask = (int *)malloc((newwidth / 2) * (newheight / 2) * sizeof(int));
      double *lev2image = NULL;
      lev2image =
          (double *)malloc((newwidth / 4) * (newheight / 4) * sizeof(double));
      int *lev2mask = NULL;
      lev2mask = (int *)malloc((newwidth / 4) * (newheight / 4) * sizeof(int));

      /* transfer the image to the double array, extending if necessary */
      for (j = 0; j < newwidth * newheight; j++) {
        newimage[j] = 0.0;
        newmask[j] = 0;
      }
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          ind = x + y * width;
          jnd = x + y * newwidth;
          newimage[jnd] = (float)object[k].image[ind];
          newmask[jnd] = object[k].mask[ind];
        }
      }
      /* store non-mask pixels as real numbers */
      object[k].lev0num = 0;
      for (x = 0; x < newwidth * newheight; x++) {
        if (newmask[x] != 0) {
          object[k].lev0Pix[object[k].lev0num] = newimage[x];
          object[k].lev0num++;
        }
      }
      /* do both level wavelet transforms and store non-mask pixels (as real
       * numbers) */
      waveTran2D(newwidth, newheight, newimage, lev1image);
      shrinkmask(newwidth / 2, newheight / 2, newmask, lev1mask);
      object[k].lev1num = 0;
      for (x = 0; x < (newwidth / 2) * (newheight / 2); x++) {
        if (lev1mask[x] != 0) {
          object[k].lev1Pix[object[k].lev1num] = lev1image[x];
          object[k].lev1num++;
        }
      }
      waveTran2D(newwidth / 2, newheight / 2, lev1image, lev2image);
      shrinkmask(newwidth / 4, newheight / 4, lev1mask, lev2mask);
      object[k].lev2num = 0;
      for (x = 0; x < (newwidth / 4) * (newheight / 4); x++) {
        if (lev2mask[x] != 0) {
          object[k].lev2Pix[object[k].lev2num] = lev2image[x];
          object[k].lev2num++;
        }
      }

      /* calculate co-occurrence matrix between image and first level wavelet
       * approximation */
      double *cooc01 = NULL;
      cooc01 = (double *)malloc(cooccurrence_levels * cooccurrence_levels *
                                sizeof(double));
      getCoocMatrix(cooc01, cooccurrence_levels, newimage, newmask, lev1image,
                    lev1mask, newwidth / 2, newheight / 2, 1);

      /* calculate co-occurrence matrix between first and second level wavelet
       * approximations */
      double *cooc12 = NULL;
      cooc12 = (double *)malloc(cooccurrence_levels * cooccurrence_levels *
                                sizeof(double));
      getCoocMatrix(cooc12, cooccurrence_levels, lev1image, lev1mask, lev2image,
                    lev2mask, newwidth / 4, newheight / 4, 1);

      /* calculate co-occurrence matrix between image and second level wavelet
       * approximation */
      double *cooc02 = NULL;
      cooc02 = (double *)malloc(cooccurrence_levels * cooccurrence_levels *
                                sizeof(double));
      getCoocMatrix(cooc02, cooccurrence_levels, newimage, newmask, lev2image,
                    lev2mask, newwidth / 4, newheight / 4, 2);

      for (x = 0; x < cooccurrence_levels * cooccurrence_levels; x++) {
        object[k].cooc01[x] = cooc01[x];
        object[k].cooc12[x] = cooc12[x];
        object[k].cooc02[x] = cooc02[x];
      }
      free(newimage);
      free(newmask);
      free(lev1image);
      free(lev1mask);
      free(lev2image);
      free(lev2mask);
      free(cooc01);
      free(cooc12);
      free(cooc02);
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:  haralick
 * Description: calculates Haralick features from a co-occurrence matrix
 *----------------------------------------------------------------------------*/
void haralick(double *cooc, int cooccurrence_levels, int num, double *hf) {
  int a, b, c;
  double energy = 0.0;
  double contrast = 0.0;
  double homogeneity = 0.0;
  double correlation = 0.0;
  double entropy = 0.0;
  double rowmean = 0.0;
  double colmean = 0.0;
  double sdevx = 0.0;
  double sdevy = 0.0;
  double corrsum = 0.0;

  a = 0;
  for (c = 0; c < (cooccurrence_levels * cooccurrence_levels); c++) {
    energy = energy + ((cooc[c] / num) * (cooc[c] / num));
  }

  for (b = 0; b < cooccurrence_levels; b++) {
    for (c = 0; c < cooccurrence_levels; c++) {
      if (c == 0 && b > 0) {
        a = a + 9;
      }
      contrast = contrast + ((b - c) * (b - c) * cooc[b + c + a] / num);
      homogeneity =
          homogeneity + (cooc[b + c + a] / num) / (1 + (b - c) * (b - c));
      rowmean = rowmean + (b + 1) * (cooc[b + c + a] / num);
      colmean = colmean + (c + 1) * (cooc[b + c + a] / num);
      if (cooc[b + c + a] != 0) {
        entropy =
            entropy + ((cooc[b + c + a] / num) * (log(cooc[b + c + a] / num)));
      }
    }
  }

  a = 0;
  for (b = 0; b < cooccurrence_levels; b++) {
    for (c = 0; c < cooccurrence_levels; c++) {
      if (c == 0 && b > 0) {
        a = a + 9;
      }
      sdevx = sdevx + (((b + 1) - rowmean) * ((b + 1) - rowmean) *
                       cooc[b + c + a] / num);
      sdevy = sdevy + (((c + 1) - colmean) * ((c + 1) - colmean) *
                       cooc[b + c + a] / num);
      corrsum = corrsum + (((b + 1) - rowmean) * ((c + 1) - colmean) *
                           cooc[b + c + a] / num);
    }
  }
  sdevx = sqrt(sdevx);
  sdevy = sqrt(sdevy);
  correlation = corrsum / (sdevx * sdevy);
  entropy = entropy * -1.0;

  hf[0] = energy;
  hf[1] = contrast;
  hf[2] = homogeneity;
  hf[3] = correlation;
  hf[4] = entropy;
}

/*------------------------------------------------------------------------------
 * Procedure:  cooccurVariables
 * Description: extract Haralick features from co-occurrence matrices between
 *              different wavelet levels
 *----------------------------------------------------------------------------*/
void cooccurVariables(AREA *object, int cooccurrence_levels, INVARS *input,
                      int inputnum, int nframes, int *missingframe) {
  int k;

  double *hf = NULL;
  hf = (double *)malloc(5 * sizeof(double));

  for (k = 0; k < nframes; k++) {
    if (missingframe[k] != 1) {
      haralick(object[k].cooc01, cooccurrence_levels, object[k].lev0num, hf);
      input[inputnum].frame[k] = hf[0];
      input[inputnum + 1].frame[k] = hf[1];
      input[inputnum + 2].frame[k] = hf[2];
      input[inputnum + 3].frame[k] = hf[3];
      input[inputnum + 4].frame[k] = hf[4];

      haralick(object[k].cooc12, cooccurrence_levels, object[k].lev1num, hf);
      input[inputnum + 5].frame[k] = hf[0];
      input[inputnum + 6].frame[k] = hf[1];
      input[inputnum + 7].frame[k] = hf[2];
      input[inputnum + 8].frame[k] = hf[3];
      input[inputnum + 9].frame[k] = hf[4];

      haralick(object[k].cooc02, cooccurrence_levels, object[k].lev0num, hf);
      input[inputnum + 10].frame[k] = hf[0];
      input[inputnum + 11].frame[k] = hf[1];
      input[inputnum + 12].frame[k] = hf[2];
      input[inputnum + 13].frame[k] = hf[3];
      input[inputnum + 14].frame[k] = hf[4];
    } else {
      input[inputnum].frame[k] = -1.0;
      input[inputnum + 1].frame[k] = -1.0;
      input[inputnum + 2].frame[k] = -1.0;
      input[inputnum + 3].frame[k] = -1.0;
      input[inputnum + 4].frame[k] = -1.0;
      input[inputnum + 5].frame[k] = -1.0;
      input[inputnum + 6].frame[k] = -1.0;
      input[inputnum + 7].frame[k] = -1.0;
      input[inputnum + 8].frame[k] = -1.0;
      input[inputnum + 9].frame[k] = -1.0;
      input[inputnum + 10].frame[k] = -1.0;
      input[inputnum + 11].frame[k] = -1.0;
      input[inputnum + 12].frame[k] = -1.0;
      input[inputnum + 13].frame[k] = -1.0;
      input[inputnum + 14].frame[k] = -1.0;
    }
  }

  free(hf);
}

/*------------------------------------------------------------------------------
 * Procedure:  interpolate
 * Description:  gives values to missing frames by interpolation
 *----------------------------------------------------------------------------*/
void interpolate(INVARS *input, int numinput, int nframes, int *missingframe) {
  int j, k, r, m;
  double value;

  for (k = 1; k < nframes; k++) {
    r = 0;
    if (missingframe[k] == 1) {
      r = 1;
      while (missingframe[k + r] == 1)
        r++;
    }

    if (r > 0) {
      for (j = 0; j < numinput; j++) {
        value =
            (input[j].frame[k + r] - input[j].frame[k - 1]) / (float)(r + 1);
        for (m = 0; m < r; m++) {
          input[j].frame[k + m] = input[j].frame[k - 1] + (m + 1) * value;
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:    wavevars
 * Description:  calculate variables from wavelet coefficients
 *----------------------------------------------------------------------------*/
void wavevars(double *inputArray, int start, int end, double *vars, int lev) {
  int k, n;
  int length = end - start;
  double diff, ascent, descent, max;

  double *tempArray = NULL;
  tempArray = (double *)malloc(length * sizeof(double));

  n = 0;
  for (k = start; k < end; k++) {
    tempArray[n] = inputArray[k];
    n++;
  }

  /* calculate total ascent and descent */
  ascent = 0.0;
  descent = 0.0;
  max = 0.0;
  for (k = 1; k < length; k++) {
    diff = tempArray[k] - tempArray[k - 1];
    if (diff > 0)
      ascent += diff;
    else
      descent += diff;

    if (fabs(tempArray[k]) > max)
      max = fabs(tempArray[k]);
  }
  vars[3 * (lev + 1)] = ascent / (float)length;
  vars[3 * (lev + 1) + 1] = descent / (float)length;
  vars[3 * (lev + 1) + 2] = max / (float)length;
  free(tempArray);
}

/*------------------------------------------------------------------------------
 * Procedure: summarystats
 * Description: calculates various stats from cell variables
 *----------------------------------------------------------------------------*/
void summarystats(INVARS *input, int numinput, int nframes, int *missingframe) {
  int j, k, nummissing, notmissing;
  double m1, m2, m3, en, enm1, enm2;

  nummissing = 0;
  for (j = 0; j < nframes; j++) {
    if (missingframe[j] == 1)
      nummissing++;
  }
  notmissing = nframes - nummissing;
  if (notmissing > 5) {
    en = (float)(notmissing);
    enm1 = en - 1.0;
    enm2 = en - 2.0;
    for (k = 0; k < numinput; k++) {
      /* get means */
      input[k].stats[0] = 0.0;
      for (j = 0; j < nframes; j++) {
        if (missingframe[j] != 1) {
          input[k].stats[0] += input[k].frame[j];
        }
      }
      /* stat 0 is mean */
      input[k].stats[0] /= en;

      /* get moments */
      m1 = input[k].stats[0];
      m2 = 0.0;
      m3 = 0.0;
      for (j = 0; j < nframes; j++) {
        if (missingframe[j] != 1) {
          m2 += (input[k].frame[j] - m1) * (input[k].frame[j] - m1);
          m3 += (input[k].frame[j] - m1) * (input[k].frame[j] - m1) *
                (input[k].frame[j] - m1);
        }
      }
      /* stat 1 is standard deviation */
      input[k].stats[1] = sqrt(m2 / enm1);
      /* stat 2 is skewness */
      if (m2 == 0.0)
        input[k].stats[2] = 0.0;
      else
        input[k].stats[2] = (m3 * en * sqrt(enm1)) / (enm2 * m2 * sqrt(m2));
    }
  }
}

/*------------------------------------------------------------------------------
 * Procedure:    timeSeriesVars
 * Description:  calculates variables from time series
 *----------------------------------------------------------------------------*/
double timeSeriesVars(INVARS *input, int max_number_of_frames,
                      int number_of_wavelet_levels, int numinput, int nframes,
                      int *missingframe) {
  int j, k;
  double diff, areaoftraject;
  double minX, maxX, minY, maxY, area;

  double *wave = NULL;
  wave = (double *)malloc((max_number_of_frames) * sizeof(double));
  int *wl = NULL;
  wl = (int *)malloc((number_of_wavelet_levels + 1) * sizeof(int));

  /* calculate total ascent and descent */
  for (j = 0; j < numinput; j++) {
    input[j].vars[2] = input[j].frame[0];
    for (k = 1; k < nframes; k++) {
      diff = input[j].frame[k] - input[j].frame[k - 1];
      if (diff > 0)
        input[j].vars[0] += diff;
      else
        input[j].vars[1] += diff;
      if (input[j].frame[k] > input[j].vars[2])
        input[j].vars[2] = input[j].frame[k];
    }
    input[j].vars[0] /= (float)nframes;
    input[j].vars[1] /= (float)nframes;
    /* this doesn't need dividing by nframes - just rescaling the maximum */
    input[j].vars[2] /= (float)nframes;
  }

  /* calculate wavelet detail coefficients */
  for (j = 0; j < numinput; j++) {
    waveTran(input[j].frame, number_of_wavelet_levels, nframes, wave, wl);
    for (k = 0; k < number_of_wavelet_levels; k++) {
      wavevars(wave, wl[k], wl[k + 1], input[j].vars, k);
    }
  }

  /* calculate summary statistics */
  summarystats(input, numinput, nframes, missingframe);

  /* calculate trajectory area */
  minX = input[0].frame[0];
  maxX = 0.0;
  minY = input[1].frame[0];
  maxY = 0.0;

  for (j = 0; j < nframes; j++) {
    if (input[0].frame[j] < minX) {
      minX = input[0].frame[j];
    }
    if (input[1].frame[j] < minY) {
      minY = input[1].frame[j];
    }
    if (input[0].frame[j] > maxX) {
      maxX = input[0].frame[j];
    }
    if (input[1].frame[j] > maxY) {
      maxY = input[1].frame[j];
    }
  }
  area = (float)((maxX - minX) * (maxY - minY));
  areaoftraject = area / (float)nframes;

  free(wave);
  free(wl);
  return (areaoftraject);
}

/*------------------------------------------------------------------------------
 * Procedure: writedata
 * Description: write out all values for particular variables
 *----------------------------------------------------------------------------*/
void writedata(INVARS *input, int number_of_wavelet_levels, int numinput,
               int nstats, double areaoftraject, NAMES *vname,
               char *classlabel) {
  int i, j, k;
  char varname[100];
  char levname[100];

  FILE *fp = NULL;
  fp = fopen("outputdata.txt", "w");

  NAMES *str = NULL;
  str = (NAMES *)malloc(3 * sizeof(NAMES));
  strcpy(str[0].var, "_asc\0");
  strcpy(str[1].var, "_des\0");
  strcpy(str[2].var, "_max\0");

  NAMES *lev = NULL;
  lev = (NAMES *)malloc(3 * sizeof(NAMES));
  strcpy(lev[0].var, "_l1\0");
  strcpy(lev[1].var, "_l2\0");
  strcpy(lev[2].var, "_l3\0");

  NAMES *stat = NULL;
  stat = (NAMES *)malloc(nstats * sizeof(NAMES));
  strcpy(stat[0].var, "_mean\0");
  strcpy(stat[1].var, "_std\0");
  strcpy(stat[2].var, "_skew\0");

  fprintf(fp, "class trajarea ");
  for (j = 2; j < numinput; j++) {
    for (k = 0; k < 3; k++) {
      strcpy(varname, vname[j].var);
      strcat(varname, str[k].var);
      if ((j != 9) || ((j == 9) && (k != 1)))
        fprintf(fp, "%s ", varname);
    }
    for (k = 0; k < 3; k++) {
      strcpy(varname, vname[j].var);
      strcat(varname, lev[k].var);
      for (i = 0; i < 3; i++) {
        strcpy(levname, varname);
        strcat(levname, str[i].var);
        if ((j != 9) || ((j == 9) && (i != 1)))
          fprintf(fp, "%s ", levname);
      }
    }
  }
  for (j = 2; j < numinput; j++) {
    for (k = 0; k < nstats; k++) {
      strcpy(varname, vname[j].var);
      strcat(varname, stat[k].var);
      fprintf(fp, "%s ", varname);
    }
  }
  fprintf(fp, "\n");

  fprintf(fp, "%s %f ", classlabel, areaoftraject);
  for (j = 2; j < numinput; j++) {
    for (k = 0; k < 3; k++) {
      if ((j != 9) || ((j == 9) && (k != 1)))
        fprintf(fp, "%f ", input[j].vars[k]);
    }
    for (k = 6; k < 3 * (number_of_wavelet_levels + 1); k++) {
      if ((j != 9) || ((j == 9) && (k % 3 != 1)))
        fprintf(fp, "%f ", input[j].vars[k]);
    }
  }
  for (j = 2; j < numinput; j++) {
    for (k = 0; k < nstats; k++) {
      fprintf(fp, "%f ", input[j].stats[k]);
    }
  }
  fprintf(fp, "\n");

  free(stat);
  free(lev);
  free(str);
  fclose(fp);
}
