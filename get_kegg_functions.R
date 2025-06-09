# Defining function to get kegg_ids ----
# df                  dataframe with molecular formulas
# id_column           column with unique metabolite IDS
# mol_formula_column  column with molecular formulas
# match_names         wheter to try to match kegg compounds by name and formula
# names_column        if match_names is TRUE, column with metabolite names

get_keggs_id <- function(df,
                         id_column,
                         molformula_column,
                         names_column = NULL,
                         match_names = FALSE){
  
  if(match_names){
    if(is.null(names_column)){
      stop('You must define a column with names for matching')
    }
  }
  
  ids <- df[[id_column]]
  formulas <- df[[molformula_column]]
  
  if(is.null(names_column)){
    des <- rep(NA, nrow(df))
  } else des <- df[[names_column]]
  
  kegg_info <- pmap(
    list(ids, formulas, des), 
    function(id, formula, des){
      
      print(paste0('Working on ', formula))
      
      cpd_ids <- names(keggFind('compound', formula, 'formula'))
      
      if(!is.null(cpd_ids)){
        searches <- map(cpd_ids, function(k_id){
          cpd_info <- keggGet(k_id)
          
          df <- tibble(
            KEGG_id = k_id,
            KEGG_name = paste(cpd_info[[1]]$NAME, collapse = ''),
            KEGG_formula = cpd_info[[1]]$FORMULA,
            KEGG_pathway = ifelse(!is.null(cpd_info[[1]]$PATHWAY), paste(cpd_info[[1]]$PATHWAY, collapse = ';'), NA),
            KEGG_module = ifelse(!is.null(cpd_info[[1]]$MODULE), paste(cpd_info[[1]]$MODULE, collapse = ';'), NA),
            KEGG_brite = ifelse(!is.null(cpd_info[[1]]$BRITE), paste(cpd_info[[1]]$BRITE, collapse = ';'), NA),
            KEGG_enzyme = ifelse(!is.null(cpd_info[[1]]$ENZYME), paste(cpd_info[[1]]$ENZYME, collapse = ';'), NA),
            KEGG_reaction = ifelse(!is.null(cpd_info[[1]]$REACTION), paste(cpd_info[[1]]$REACTION, collapse = ';'), NA)
          )
        })
        
        cpd_df <- do.call(rbind, searches) %>% 
          mutate(id = id,
                 Formula = formula,
                 Name = des)
      } else {
        
        cpd_df <- tibble(
          KEGG_id = NA,
          KEGG_name = NA,
          KEGG_formula = NA,
          KEGG_pathway = NA,
          KEGG_module = NA,
          KEGG_brite = NA,
          KEGG_enzyme = NA,
          KEGG_reaction = NA
        ) %>% 
          mutate(id = id,
                 Formula = formula,
                 Name = des)
      }
      
      return(cpd_df)
    })
  
  kegg_info <- do.call(rbind, kegg_info) %>% 
    filter(!is.na(KEGG_id),
           Formula == KEGG_formula)
  
  if(match_names){
    kegg_info <- kegg_info %>% 
      mutate(regex_description = str_replace_all(Name, '\\{', '\\('),
             regex_description = str_replace_all(regex_description, '\\}', '\\)'),
             regex_description = str_replace_all(regex_description, '\\[', '\\('),
             regex_description = str_replace_all(regex_description, '\\]', '\\)'),
             regex_description = str_remove_all(regex_description, '\\(\\+\\)-')) %>% 
      filter(str_detect(KEGG_name, regex_description))
    
  } else {
    
    if(is.null(names_column)){
      kegg_info <- kegg_info %>% 
        select(-Name)
    } else {
      kegg_info
    }  
    
    
  }
  
  return(kegg_info)
  
}

# Defining function to get pathways of a specific organism ----

get_paths_org <- function(org_code){
  
  path_list <- keggLink(org_code, 'pathway')
  path_list <- unique(names(path_list))
  
  path_df <- map(path_list, function(path){
    
    info <- keggGet(path)
    
    df <- tibble(PathwayID = path,
                 Pathway_name = info[[1]]$NAME)
    
  })
  
  path_df <- do.call(rbind, path_df)
  
}


# # Using the functions ----
# 
# data <- read_csv('name_for_kegg.csv') %>% 
#   mutate(id = paste0('Metabolite_', n():1),
#          Formula = str_remove_all(Formula, ' '))
# 
# # with name matching
# 
# kegg_df_names_match <- get_keggs_id(data[11:20,],
#                                     id_column = 'id',
#                                     molformula_column = 'Formula',
#                                     names_column = 'Name',
#                                     match_names = FALSE)


get_kegg_by_name <- function(compound_names){
  
  kegg_info <- map(compound_names, function(name){
      
      print(paste0('Working on ', name))
      
      cpd_ids <- names(keggFind('compound', name))
      
      if(!is.null(cpd_ids)){
        searches <- map(cpd_ids, function(k_id){
          cpd_info <- keggGet(k_id)
          
          df <- tibble(
            query = name,
            KEGG_id = k_id,
            KEGG_name = paste(cpd_info[[1]]$NAME, collapse = ''),
            KEGG_formula = ifelse(!is.null(cpd_info[[1]]$FORMULA), paste(cpd_info[[1]]$FORMULA, collapse = ';'), NA),
            KEGG_pathway = ifelse(!is.null(cpd_info[[1]]$PATHWAY), paste(cpd_info[[1]]$PATHWAY, collapse = ';'), NA),
            KEGG_module = ifelse(!is.null(cpd_info[[1]]$MODULE), paste(cpd_info[[1]]$MODULE, collapse = ';'), NA),
            KEGG_brite = ifelse(!is.null(cpd_info[[1]]$BRITE), paste(cpd_info[[1]]$BRITE, collapse = ';'), NA),
            KEGG_enzyme = ifelse(!is.null(cpd_info[[1]]$ENZYME), paste(cpd_info[[1]]$ENZYME, collapse = ';'), NA),
            KEGG_reaction = ifelse(!is.null(cpd_info[[1]]$REACTION), paste(cpd_info[[1]]$REACTION, collapse = ';'), NA)
          )
        })
        
        cpd_df <- reduce(searches, rbind)
        
      } else {
        
        cpd_df <- tibble(
          query = name,
          KEGG_id = NA,
          KEGG_name = NA,
          KEGG_formula = NA,
          KEGG_pathway = NA,
          KEGG_module = NA,
          KEGG_brite = NA,
          KEGG_enzyme = NA,
          KEGG_reaction = NA
        )
      }
      
      return(cpd_df)
    })
  
  final <- reduce(kegg_info, rbind)
}


