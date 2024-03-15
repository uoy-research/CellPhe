#' Export CellPhe frame-specific features to a SQLite database that can be read by CellProfiler Analyst
#'
#' @description CellProfiler Analyst can read the frame-level features from CellPhe for data visualisation
#'   and exploration.
#' @param features A feature table output by \code{extractFeatures}.
#' @param image_folder The folder where the images are stored.
#' @param image_prefix The code assumes that images are labelled <prefix>_<x>.tif, where \code{prefix} is
#'   \code{image_prefix} and \code{x} is the frameID 0-padded to 4 digits.
#' @param db_filename The name of the SQLite database to save to.
#' @param properties_filename The name of the CellProfiler Analyst properties file to save to.
#' @return Doesn't return anything, but instead creates a SQLite database and a plain text properties file.
#'
#' @export
export_to_cellprofiler_analyst <- function(features,
                                           image_folder,
                                           image_prefix,
                                           db_filename="cellphe.db",
                                           properties_filename="cellphe.properties") {
  
  # The per_image table needs
  #  - image id
  #  - image path (1 column per channel)
  #  - image filename (1 column per channel)
  #  - Image_Group_Index (can be same as id for now)
  #  - Image_Group_Number (can all 1s)
  # We only have 1 channel
  per_image <- features |>
                  dplyr::distinct(FrameID) |>
                  dplyr::mutate(
                    path_channel_1 = image_folder,
                    filename_channel_1 = sprintf("%s-%04d.tif", image_prefix, FrameID),
                    Image_Group_Index = FrameID,
                    Image_Group_Number = 1
                  )
  
  # the per_object table requires the following columns:
  #  - image id (FrameID)
  #  - object id (CellID)
  #  - centroid x (xpos)
  #  - centroid y (ypos)
  #  - ... Any number of feature columns
  # Drop ROI filename as it doesn't fit into any of these categories
  features <- features[, setdiff(colnames(features), "ROI_filename")]
  cell_fields <- setdiff(colnames(features), c("FrameID", "CellID", "xpos", "ypos"))
  
  # We also need 3 more tables:
  #  - FramePer_Relationships: mapping objects in sequential frames
  frame_per_relationships <- features |>
    dplyr::distinct(FrameID, CellID) |>
    dplyr::inner_join(features |> dplyr::distinct(FrameID, CellID),
               dplyr::join_by(CellID == CellID, FrameID < FrameID), suffix = c("1", "2")) |>
    dplyr::group_by(CellID, FrameID1) |>
    dplyr::top_n(-1, FrameID2) |>
    dplyr::mutate(CellID2 = CellID, relationship_type_id = 1) |>
    dplyr::select(
      relationship_type_id,
      image_number1 = FrameID1, 
      object_number1 = CellID,
      image_number2 = FrameID2, 
      object_number2 = CellID2
    )
  #  - FramePer_RelationshipTypes: Possible relationship types
  frame_per_relationship_types <- dplyr::tribble(
    ~relationship_type_id, ~module_number, ~relationship, ~object_name1, ~object_name2,
    # TODO make cell a constant as use it in properties if this works
    1, 1, "Parent", "cell", "cell"
  )
  
  #  - FramePer_RelationshipsView: A view combining the two tables above
  
  # Create DB connection
  con <- DBI::dbConnect(RSQLite::SQLite(), db_filename)
  
  # Setup tables
  DBI::dbCreateTable(con, 
                     "frames", 
          
                                c(
                       "FrameID"="INTEGER", 
                       "path_channel_1"="TEXT",
                       "filename_channel_1"="TEXT",
                       "Image_Group_Number"="INTEGER",
                       "Image_Group_Index"="INTEGER"
                       )
                     )
  DBI::dbCreateTable(con, 
                     "cells",
                     c(
                       "FrameID"="INTEGER",
                       "CellID"="INTEGER",
                       "xpos"="REAL",
                       "ypos"="REAL",
                       setNames(
                         rep("REAL", length(cell_fields)),
                         cell_fields
                       )
                     ))
  
  DBI::dbCreateTable(
    con,
    "FramePer_Relationships",
    c(
      "relationship_type_id"="INTEGER",
      "image_number1"="INTEGER",
      "object_number1"="INTEGER",
      "image_number2"="INTEGER",
      "object_number2"="INTEGER"
    )
  )
  
  DBI::dbCreateTable(
    con,
    "FramePer_RelationshipTypes",
    c(
      "relationship_type_id"="INTEGER",
      "module_number"="INTEGER",
      "relationship"="TEXT",
      "object_name1"="TEXT",
      "object_name2"="TEXT"
    )
  )
  
  DBI::dbExecute(
    con,
    "
    CREATE VIEW FramePer_RelationshipsView AS 
    SELECT module_number, relationship, object_name1, object_name2, 
           image_number1, object_number1, image_number2, object_number2
    FROM FramePer_Relationships fpr
    INNER JOIN FramePer_RelationshipTypes fpt
    ON fpr.relationship_type_id = fpt.relationship_type_id;
  ")
  
  # Populate tables
  DBI::dbAppendTable(con, "frames", per_image)
  DBI::dbAppendTable(con, "cells", features)
  DBI::dbAppendTable(con, "FramePer_Relationships", frame_per_relationships)
  DBI::dbAppendTable(con, "FramePer_RelationshipTypes", frame_per_relationship_types)
  
  DBI::dbDisconnect(con)
  
  write_properties(properties_filename, db_filename)
}

write_properties <- function(properties_filename, db_filename) {
  raw <- "
  db_type      =  sqlite
  db_sqlite_file  =  %s
  
  image_table   =  frames
  object_table  =  cells
  
  image_id    =  FrameID
  object_id   =  CellID
  
  cell_x_loc  =  xpos
  cell_y_loc  =  ypos
  
  image_path_cols  =  path_channel_1
  image_file_cols  =  filename_channel_1
  
  image_names   =  channel1
  image_channel_colors  =  red
  
  image_channel_blend_modes = add
  
  # ======== Dynamic Groups ========
  # OPTIONAL
  # Here you can define ways of grouping your image data, by linking column(s)
  # that identify unique images (the image-key) to a unique group of columns the
  # (group-key). Note that the group-key columns may come from other tables, so 
  # long as the tables have a common key.
  # 
  # Example: With the 'Well' group defined below, Classifier will allow you to
  #   fetch cells from images from a particular well by providing you with a well
  #   dropdown in the user interface. It will also allow you to group your data
  #   by each unique well value when scoring.
  #
  # Example 2: Also note the 'Plate_and_Well' group. This group specifies unique
  #   pairs of plate and well values. Since well values such as 'A01' are likely
  #   to NOT be unique across multiple plates, this will provide a way to refer
  #   to cells from, plate X, well A01, rather than just any well named 'A01'.
  #  
  # FORMAT:
  #   group_XXX  =  MySQL select statement that returns image-key columns followed by group-key columns. XXX will be the name of the group.
  # EXAMPLE GROUPS:
  #   group_SQL_Well        =  SELECT TableNumber, ImageNumber, well FROM Per_Image_Table
  #   group_SQL_Plate_and_Well  =  SELECT Per_Image_Table.TableNumber, Per_Image_Table.ImageNumber, Well_ID_Table.Plate, Per_Image_Table.well FROM Per_Image_Table, WELL_ID_Table WHERE Per_Image_Table.well=Well_ID_Table.well
  #   group_SQL_Treatment   =  SELECT Per_Image_Table.TableNumber, Per_Image_Table.ImageNumber, Well_ID_Table.treatment FROM Per_Image_Table, Well_ID_Table WHERE Per_Image_Table.well=Well_ID_Table.well
  
  #group_SQL_YourGroupName  =  
    
    
  # ======== Image Filters ========
  # OPTIONAL
  # Here you can define image filters to let you fetch or score objects from a 
  # subset of the images in your experiment.
  #
  # Example: With the CDKs filter defined below, Classifier will provide an extra
  #   option to fetch cells from CDKs... that is, images who's corresponding gene
  #   entry starts with CDK.
  #
  # FORMAT:
  #   filter_SQL_XXX  =  MySQL select statement that returns image-key columns for images you wish to filter out. XXX will be the name of the filter.
  # EXAMPLE FILTERS:
  #   filter_SQL_EMPTY  =  SELECT TableNumber, ImageNumber FROM CPA_per_image, Well_ID_Table WHERE CPA_per_image.well=Well_ID_Table.well AND Well_ID_Table.Gene='EMPTY'
  #   filter_SQL_CDKs   =  SELECT TableNumber, ImageNumber FROM CPA_per_image, Well_ID_Table WHERE CPA_per_image.well=Well_ID_Table.well AND Well_ID_Table.Gene REGEXP 'CDK.*'
  
  #filter_SQL_YourFilterName  =  
    
    
  # ======== Meta data ======== 
  # What are your objects called? (e.g. cells, worms, etc.)
  # This is used to provide the correct syntax for the GUI.
  # FORMAT:  object_name  =  singular name, plural name
  object_name  =  cell, cells
  
  classifier_ignore_columns  =  xpos, ypos, meta_.* 
  
  image_tile_size  =  50
  
  #image_width  = 
  #image_height = 
    
  # OPTIONAL
  # Image Gallery can use a different tile size (in pixels) to create thumbnails for images
  # If not set, it will be the same as image_tile_size
  
  #image_size = 
  
  classification_type = object
    
  #training_set  =  
    
    
  # ======== Area Based Scoring ========
  # OPTIONAL
  # You may specify a column in your per-object table which will be summed and 
  # reported in place of object-counts when scoring.  The typical use for this
  # is to report the areas of objects on a per-image or per-group basis. 
  
  #area_scoring_column  =  
    
    
  # ======== Output Per-Object Classes ========
  # OPTIONAL
  # Here you can specify a MySQL table in your Database where you would like 
  # Classifier to write out class information for each object in the 
  # object_table 
  
  #class_table  =  
    
    
  # ======== Check Tables ========
  # OPTIONAL
  # [yes/no]  You can ask CPA to check your tables for anomalies such as
  # orphaned objects or missing column indices.  Default is on.
  # This check is run when Classifier starts and may take up to a minute if
  # your object_table is extremely large.
  
  check_tables = yes
  
  
  # ======== Force BioFormats ========
  # OPTIONAL
  # [yes/no]  By default, CPA will try to use the imageio library to load images
  # which are in supported formats, then fall back to using the older BioFormats
  # loader if something goes wrong. ImageIO is faster but some unusual file
  # compression formats can cause errors when loading. This option forces CPA to
  # always use the BioFormats reader. Try this if images aren't displayed correctly.
  
  force_bioformats = no
  
  
  # ======== Use Legacy Fetcher ========
  # OPTIONAL
  # [yes/no]  In CPA 3.0 the object fetching system has been revised to be more
  # efficient. In the vast majority of cases it should be faster than the previous
  # versions. However, some complex object filters can still cause problems. If you
  # encounter slowdowns this setting allows you to switch back to the old method of
  # fetching and randomisation.
  
  use_legacy_fetcher = no
  
  # ======== Process as 3D (visualize a different z position per object) ========
  # OPTIONAL
  # [yes/no]  In 3D datasets, this optionally displays in CPA classifier a separate
  # z slice for each object depending on that object's center position in z. Useful
  # for classifying cells from 3D data.
  
  process_3D = no
  "
  processed <- sprintf(raw, db_filename)
  write(processed, properties_filename)
}
