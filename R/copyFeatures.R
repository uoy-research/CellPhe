#' Copy metadata and cell-frame features from an existing PhaseFocus or Trackmate table
#'
#' Loads the frame and cell IDs along with the filename used to refer to each ROI.
#' For PhaseFocus generated data, volume and sphericity features are also extracted.
#' Only cells that are tracked for a minimum of \code{minframes} are included.
#'
#' @param file The filepath to a CSV file containing features output by PhaseFocus or Trackmate software.
#' @param minframes The minimum number of frames a cell must be tracked for to
#' be included in the output features.
#' @param source The name of the software that produced the metadata file, either 'Phase' or 'Trackmate'
#' are currently supported.
#' @return A dataframe with 1 row corresponding to 1 cell tracked in 1 frame
#' with the following columns:
#' \itemize{
#'   \item{\code{FrameID}: the numeric FrameID}
#'   \item{\code{CellID}: the numeric CellID}
#'   \item{\code{ROI_filename}: the label used to refer to this ROI}
#'   \item{\code{Volume}: a real-valued number}
#'   \item{\code{Sphericity}: a real-valued number}
#' }
#' @export
copyFeatures = function(file,
                        minframes,
                        source = c("Phase", "Trackmate")) {
  source <- match.arg(source)
  if (source == "Phase") {
    full_ft <-
      utils::read.csv(
        file,
        header = TRUE,
        skip = 1,
        stringsAsFactors = FALSE
      )
    full_ft$ROI_filename <-
      sprintf("%s-%s", full_ft$Frame, full_ft$Tracking.ID)
    out <-
      full_ft[, c("Frame",
                  "Tracking.ID",
                  "ROI_filename",
                  "Volume..µm..",
                  "Sphericity...")]
    colnames(out) <-
      c("FrameID",
        "CellID",
        "ROI_filename",
        "Volume",
        "Sphericity")
  } else if (source == "Trackmate") {
    full_ft <- utils::read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    # Lines 2-4 in the raw file contain additional header information and can be safely discarded
    out <- full_ft[-(1:3), c("FRAME", "TRACK_ID", "LABEL")]
    colnames(out) <- c("FrameID", "CellID", "ROI_filename")
    out$FrameID <- as.integer(out$FrameID) + 1  # Convert from 0-indexed to 1-indexed
  }
  
  # Restrict to cells which are in minimum number of frames
  out_sub <- out |>
    dplyr::mutate(CellID = as.integer(CellID),
                  FrameID = as.integer(FrameID)) |>
    dplyr::group_by(CellID) |>
    dplyr::filter(dplyr::n() >= minframes) |>
    dplyr::ungroup() |>
    dplyr::arrange(CellID, FrameID)
}