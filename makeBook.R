

  library("bookdown")
  library("tidyverse")
  library("fs")
  
  file_delete("_main.Rmd")
  
  render_book("rmd")
  
  outDir <- paste0("S:/DEHTemp/",basename(here::here()),"/")
  
  dir_create(outDir)
  
  folders <- dir_info(recursive = TRUE) %>%
    dplyr::filter(grepl("_book|out",path)
                  , type == "directory"
                  ) %>%
    dplyr::mutate(toDir = paste0(outDir,path)) %>%
    dplyr::pull(toDir)
  
  #file_delete(folders)
  dir_create(folders)
  
  files <- dir_info(recursive = TRUE) %>%
    dplyr::filter(grepl("_book",path)
                  , type != "directory"
                  ) %>%
    dplyr::mutate(toDir = paste0(outDir,path))
  
  #file_delete(files$path)
  
  walk2(files$path
        , files$toDir
        , file_copy
        , overwrite = TRUE
        )
  
  file_delete("rPackages.bib")
  