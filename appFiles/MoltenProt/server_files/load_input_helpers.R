source_python("moltenprot_shiny.py")

## Count folders in the current directory
count_folders <- function(dir) { return(length(list.files(dir))) }

getFileNameExtension <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split='/')[[1]]
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l]
  # the extention must be the suffix of a non-empty name
  return(ext)
}

## Get color palette
get_colors <- function(n) {

    if (n <= 9) {
      return(global_palette_9[1:n])
    } else if (n <= 40) {
      return(global_palette_40[1:n])
    } else if (n <= 96) {
      return(global_palette_96[1:n])
    } else if (n <= 364) {
      return(global_palette_364[1:n])
    } else {
      return(rep(global_palette_364,10)[1:n])
    }
}


## Get include and conditions vectors from capillary versus condition tables

get_include_vector <- function(table1,table2,table3,table4,
                               tot_cond,row_per_table,maxConditions) {

  if (tot_cond <= (maxConditions-row_per_table*1)) {table4 <- NULL}
  if (tot_cond <= (maxConditions-row_per_table*2)) {table3 <- NULL}
  if (tot_cond <= (maxConditions-row_per_table*3)) {table2 <- NULL}

  DF1 <- hot_to_r(table1)
  include_vector    <- DF1$Include
  conditions_vector <- DF1$Condition
  series_vector     <- DF1$Series
  file_vector <- DF1$File

  color_vector        <- NULL
  color_given_by_user <- "Color" %in% colnames(DF1)


  # Find if there is a column named "Color" in the first table
  if (color_given_by_user) color_vector <- DF1$Color

  if (!(is.null(table2))) {
    DF2 <- hot_to_r(table2)
    include_vector    <- c(include_vector,DF2$Include)
    conditions_vector <- c(conditions_vector,DF2$Condition)
    series_vector     <- c(series_vector,DF2$Series)
    file_vector  <- c(file_vector,DF2$File)

    if (color_given_by_user) color_vector <- c(color_vector,DF2$Color)

  }
  
  if (!(is.null(table3))) {
    DF3 <- hot_to_r(table3) 
    include_vector    <- c(include_vector,DF3$Include)
    conditions_vector <- c(conditions_vector,DF3$Condition)
    series_vector     <- c(series_vector,DF3$Series)
    file_vector  <- c(file_vector,DF3$File)
    if (color_given_by_user) color_vector <- c(color_vector,DF3$Color)
  }
  
  if (!(is.null(table4))) {
    DF4 <- hot_to_r(table4) 
    include_vector    <- c(include_vector,DF4$Include)
    conditions_vector <- c(conditions_vector,DF4$Condition)
    series_vector     <- c(series_vector,DF4$Series)
    file_vector  <- c(file_vector,DF4$File)
    if (color_given_by_user) color_vector <- c(color_vector,DF4$Color)

  }
  
  return(list(
    "include_vector"=include_vector,
    "conditions_vector"=conditions_vector,
    "series_vector"=series_vector,
    "color_vector"=color_vector,
    "file_vector"=file_vector))
}

## Constraint the median filter value between 0 and 6

get_median_filter <- function(median_value) {
  
  median_filter <- median_value
  
  if (median_value < 0) {median_filter <- 0}
  if (median_value > 6) {median_filter <- 6}
  
  return(median_filter)
}

## Splits a vector into a list of n elements 
split_vec <- function(vector,chunck_n) {

  vector <- as.vector(vector)
  vector <- na.omit(vector)

  sels    <- list()
  chuncks <- ceiling( length(vector) /chunck_n )
  
  for (i in 1:chuncks) {
    idx <- seq(1+chunck_n*(i-1),chunck_n*i)
    sels[[i]] <- na.omit(vector[idx])
  }
  
  return(sels)
}


## Get the 4 Tables data that will store the conditions names, series and include information
get_table_data <- function(conditions,files, n_rows_conditions_table, colors=NULL, series=NULL, include=NULL) {

  data_lst <- list()

  total_cond <- length(conditions)
  d1max      <- min(n_rows_conditions_table,total_cond)

  # use default colors if not provided
  if (is.null(colors)) {
    colors <- get_colors(total_cond)
  }

  # use default series if not provided - all set to "A"
  if (is.null(series)) {
    series <- rep("A",total_cond)
  }

  # select all if not provided
  if (is.null(include)) {
    include <- rep(TRUE,total_cond)
  }

  data1 <- data.frame(Condition=as.character(conditions[1:d1max]),
                      Series=as.character(series[1:d1max]),
                      Include=as.logical(include[1:d1max]),
                      Color=as.character(colors[1:d1max]),
                      File=as.character(files[1:d1max])
                      )

  data_lst[[1]] <- data1

  if (total_cond > n_rows_conditions_table ) {

    d2max <- min(n_rows_conditions_table*2,total_cond)
    nreps <- d2max - n_rows_conditions_table

    data2 <- data.frame(
                        Condition=as.character(conditions[(n_rows_conditions_table+1):d2max]),
                        Series=as.character(series[(n_rows_conditions_table+1):d2max]),
                        Include=as.logical(include[(n_rows_conditions_table+1):d2max]),
                        Color=as.character(colors[(n_rows_conditions_table+1):d2max]),
                        File=as.character(files[(n_rows_conditions_table+1):d2max])
                        )

    data_lst[[2]] <- data2

  }

  if (total_cond > n_rows_conditions_table*2 ) {

    d3max <- min(n_rows_conditions_table*3,total_cond)
    nreps <- d3max - (n_rows_conditions_table*2)

    data3 <- data.frame(
                        Condition=as.character(conditions[(n_rows_conditions_table*2+1):d3max]),
                        Series=as.character(series[(n_rows_conditions_table*2+1):d3max]),
                        Include=as.logical(include[(n_rows_conditions_table*2+1):d3max]),
                        Color=as.character(colors[(n_rows_conditions_table*2+1):d3max]),
                        File=as.character(files[(n_rows_conditions_table*2+1):d3max])
                        )

    data_lst[[3]] <- data3

  }

  if (total_cond > n_rows_conditions_table*3 ) {

    d4max <- min(n_rows_conditions_table*4,total_cond)
    nreps <- d4max - (n_rows_conditions_table*3)

    data4 <- data.frame(
                        Condition=as.character(conditions[(n_rows_conditions_table*3+1):d4max]),
                        Series=as.character(series[(n_rows_conditions_table*3+1):d4max]),
                        Include=as.logical(include[(n_rows_conditions_table*3+1):d4max]),
                        Color=as.character(colors[(n_rows_conditions_table*3+1):d4max]),
                        File=as.character(files[(n_rows_conditions_table*3+1):d4max])
                        )

    data_lst[[4]] <- data4

  }

  return(data_lst)
}

## Get the 4 renderRHandsontable Tables
get_renderRHandsontable_list <- function(
    conditions,
    files,
    n_rows_conditions_table,
    colors=NULL,
    series=NULL,
    include=NULL,
    hide_color_column=FALSE) {

    tables_data <- get_table_data(conditions,files,n_rows_conditions_table,colors,series,include)

    tables <- lapply(seq_along(tables_data),function(i){

        df <- tables_data[[i]]

        if (hide_color_column) {

            df_no_color <- df[, names(df) != "Color", drop = FALSE]

            table <- renderRHandsontable({
                rhandsontable(df_no_color, # Remove the color column
                    rowHeaders = NULL,
                    maxRows=n_rows_conditions_table
                ) %>%
                hot_col("Condition",allowInvalid = TRUE,renderer = myrenderer) %>%
                hot_col(c("Series"),renderer = myrenderer) %>%
                hot_col("Include",renderer = myrendererBoolean) %>%
                hot_col("File",readOnly = TRUE)
            })

        } else {

            color_cells <- data.frame(col = which(names(df) == "Color"), row = 1:nrow(df))

            table <- renderRHandsontable({
                rhandsontable(df,
                    rowHeaders = NULL,
                    col_highlight = color_cells$col - 1,
                    row_highlight = color_cells$row - 1,
                    maxRows=n_rows_conditions_table
                ) %>%
                hot_col("Condition",allowInvalid = TRUE,renderer = myrenderer) %>%
                hot_col(c("Series","Color"),renderer = myrenderer) %>%
                hot_col("Include",renderer = myrendererBoolean) %>%
                hot_col("File",readOnly = TRUE)
            })

        }

        return(table)

    })

    return(tables)

}



## Divide each column of the fluorescence matrix by the median of the first 2 degrees
normalize_matrix_by_initial_value <- function(fluo_matrix,temp_vector) {

  npoints <- length(temp_vector[temp_vector < min(temp_vector)+2])

  initial_vales <- apply(fluo_matrix[1:npoints,,drop = FALSE], 2, median)

  fluo_matrix_norm   <- t(t(fluo_matrix) / initial_vales)

  return(fluo_matrix_norm)
}

## Normalize each column of the fluorescence matrix by the maximum and minimum value
normalize_matrix_max_min <- function(fluo_matrix) {

  max_values <- apply(fluo_matrix, 2, max)
  min_values <- apply(fluo_matrix, 2, min)

  delta <- max_values - min_values

  fluo_matrix_norm   <- t(t(fluo_matrix) - min_values)
  fluo_matrix_norm   <- t(t(fluo_matrix_norm) / delta)

  return(fluo_matrix_norm)
}

## Normalize each column of the fluorescence matrix such that the area under the curve is 1
normalize_matrix_area <- function(fluo_matrix,temps) {

  trapz      <- apply(fluo_matrix, 2, FUN = function(y) trapz(temps,y))
  fluo_matrix_norm   <- t(t(fluo_matrix) / trapz)

  return(fluo_matrix_norm)
}

## Normalize fluo matrix according to selection

normalize_fluo_matrix_by_option <- function(option,fluo,temps) {
  
  if (option == "Divide_by_init_value") {
    fluo <- normalize_matrix_by_initial_value(fluo,temps)}
  
  if (option == "MC_Normalization") {
    fluo <- normalize_matrix_max_min(fluo)}
  
  if (option == "area_Normalization") {
    fluo <- normalize_matrix_area(fluo,temps)}
  
  return(fluo)
}
