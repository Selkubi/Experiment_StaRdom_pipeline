library(tidyverse)

absorbance_Read=function (absorbance_path, order = TRUE, recursive = TRUE, dec = NULL, 
          sep = NULL, verbose = FALSE, cores = parallel::detectCores(logical = FALSE), 
          ...) 
{
  if (dir.exists(absorbance_path)) {
    abs_data <- list.files(absorbance_path, full.names = TRUE, 
                           recursive = recursive, no.. = TRUE, include.dirs = FALSE, 
                           pattern = "*.txt|*.dx", ignore.case = TRUE)
    abs_data <- abs_data[!file.info(abs_data)$isdir]
  }
  else if (file.exists(absorbance_path)) {
    abs_data <- absorbance_path
  }
  else stop("Absorbance data was not found!")
  if (length(abs_data) < 1) 
    stop("No valid files found in absorbance_path!")
  cl <- makeCluster(min(cores, length(abs_data)), type = "PSOCK")
  clusterExport(cl, c("dec", "sep", "verbose"), envir = environment())
  clusterEvalQ(cl, require(dplyr))
  clusterEvalQ(cl, require(stringr))
  abs_data <- parLapply(cl, abs_data, function(tab) {
    tryCatch({
      rawdata <- readLines(tab)
      data <- rawdata %>% sapply(str_remove, pattern = "([^0-9]*$)")
      first_number <- min(which((substr(data, 1, 1) %>% 
                                   grepl("[0-9]", .))))
      last_number <- max(which((substr(data, 1, 1) %>% 
                                  grepl("[0-9]", .))))
      if (is.null(sep) | is.null(dec)) {
        nsepdec <- data[first_number] %>% str_extract_all("[^-0-9eE]") %>% 
          unlist()
        example_number <- data[first_number] %>% str_extract("([-]?[0-9]+[.,]?[0-9]+[eE]?[-0-9]+)$")
        if (is.null(dec) & length(nsepdec) > 1) 
          dec <- example_number %>% str_replace("([-0-9eE]+)([.,]?)([-0-9eE]*)", 
                                                "\\2")
        if (is.null(sep)) 
          sep <- gsub(pattern = dec, replacement = "", 
                      x = data[first_number], fixed = TRUE) %>% 
          str_extract(paste0("[^-0-9eE", dec, "]"))
        if (verbose) 
          warning("processing", tab, ": using", sep, 
                  "as field separator and", dec, "as decimal separator.", 
                  fill = TRUE)
      }
      data <- str_split(data, sep)
      table <- data[(first_number):last_number] %>% unlist() %>% 
        matrix(ncol = length(data[[first_number]]), byrow = TRUE) %>% 
        data.frame(stringsAsFactors = FALSE) %>% mutate_all(gsub, 
                                                            pattern = ifelse(dec != "", dec, "."), replacement = ".", 
                                                            fixed = TRUE)
      table <- table %>% mutate_all(as.numeric)
      attr(table, "location") <- rep(tab, ncol(table) - 
                                       1)
      if (ncol(table) == 2) {
        samples <- tab %>% basename() %>% str_replace_all(regex(".txt$|.dx$", 
                                                                ignore_case = TRUE), "")
      }
      else {
        samples <- rawdata[[1]] %>% str_split(sep) %>% 
          unlist() %>% matrix(ncol = length(.), byrow = TRUE) %>% 
          data.frame(stringsAsFactors = FALSE) %>% .[-1]
      }
      table <- table %>% setNames(c("wavelength", samples))
    }, error = function(err) {
      stop("Error while reading ", tab, ": ", err)
    })
  })
  stopCluster(cl)
  locations <- lapply(abs_data, function(tab) {
    attr(tab, "location")
  }) %>% unlist()
  if (length(abs_data) == 1) 
    abs_data <- abs_data[[1]] %>% as.data.frame()
  else abs_data <- abs_data %>% list_join(by = "wavelength")
  if (order) 
    abs_data <- abs_data %>% arrange(wavelength)
  attr(abs_data, "location") <- locations
  abs_data
}


import_fluoromax4_reverse <- function(file) {
  dat <- fread(file)
  dat <- dat[-2,]
  ex <- as.matrix(as.numeric(na.omit(dat[1, -1])))
  em <- as.matrix(as.numeric(unlist(na.omit(dat[-1,1]))))

  x <- dat[, -1]
  x <- x[-1, ]
  x <- matrix(as.numeric(unlist(x, use.names = FALSE)), ncol = 41, byrow = FALSE)

  l <- list(
    file = file,
    x = x,
    ex = ex,
    em = em
  )
  
  return(l)
}

