library(tidyverse)
library(DBI)
library(MSnbase)
library(readxl)

# Load files ----

## Link file
compound_link <- read_xlsx('data/CD_results_palsa_1.5e6_old_al_node.xlsx') %>% 
  select(Name, `Calc. MW`, `m/z`, `RT [min]`, `Reference Ion`) %>% 
  mutate(MW = formatC(round(`Calc. MW`, 5), digits = 5, format = 'f'),
         RT = formatC(round(`RT [min]`, 3), digits = 3, format = 'f'),
         FeatureID = paste0('MW_', MW, '@RT_', RT),
         Name = ifelse(is.na(Name), FeatureID, Name))

## Database
msdb <- dbConnect(RSQLite::SQLite(), 'data/palsa_1.5e6_old_al_node_fix_names.db')

### Acquire the necessary information from the database
#### Compound information:
compound_info <- dbGetQuery(msdb, 'SELECT CompoundId, Formula, Name 
                            FROM CompoundTable;')
#### Spectra InformationL:
spectra_info <- dbGetQuery(msdb, 'SELECT SpectrumId, CompoundId, RetentionTime, PrecursorMass, ScanNumber, RawFileURL 
                           FROM SpectrumTable;')
#### Close the connection to the database:
dbDisconnect(msdb)

# Data wrangling ----

## Dataframe with all information
#### Rawfile URL is modified to match file paths in my laptop
total_info <- left_join(compound_info, spectra_info) 

#### Creating a dataframe with the system calls
total_nest <- total_info %>% 
  group_by(RawFileURL) %>%
  nest() %>% 
  mutate(scans = map_chr(data, function(x){
    sc <- x %>%
      drop_na() %>%
      dplyr::pull(ScanNumber) %>%
      unique() %>%
      paste0(., collapse = ',')
    return(sc) }),
    new_path = str_replace(RawFileURL, '.*Palsa\\\\', 
                           'Z:/PhD_stuff/moira_data/Palsa/'),
    call = paste0("C:/Users/Chris/Desktop/UofA/ThermoRawFileParser1.4.2/ThermoRawFileParser.exe query -i=", 
                  new_path,
                  " -n=\"",
                  scans,
                  "\" -s")) 

## Extracting spectra using the ThermoRawFileParser
spectra <- map_df(total_nest$call, function(x){
  specframe <- jsonlite::fromJSON(system(x, intern = T)) %>%
    mutate(ScanNumber = map_chr(attributes, function(x) x %>% filter(name == 'scan number') %>% pull(value))) %>%
    modify_at('ScanNumber', as.integer)
})

## Combining spectra data with compounds table
spectra_linked <- spectra %>%
  left_join(total_info, relationship = "many-to-many") %>% 
  mutate(Name = str_trim(Name)) %>% 
  left_join(compound_link, relationship = 'many-to-many') %>%
  filter(!is.na(FeatureID))

## Get consensus spectra
spectra_combined <- spectra_linked %>%
  group_by(CompoundId) %>%
  nest() %>% 
  mutate(spec_list = purrr::map(data, function(s) {
    map2(s$mzs, s$intensities, function(x,y) new('Spectrum2', mz = x, intensity = y))}))  %>%
  mutate(m_spec = purrr::map(spec_list, MSpectra)) %>%
  mutate(con_spec = purrr::map(m_spec, function(x){
    mz <- numeric(0)
    prop <- 0.5
    while(identical(mz, numeric(0)) && prop > 0.3){
      s <- consensusSpectrum(x, minProp = prop, ppm = 1)
      prop <- prop - 0.05
    }
    
    return(x)
  }))  %>%
  mutate(mz_list = purrr::map(con_spec, mz)) %>%
  mutate(int_list = purrr::map(con_spec, intensity)) %>%
  #dplyr::select(-spec_list, -m_spec, -con_spec) %>%
  mutate(reduced = purrr::map(data, function(x) x[1,])) %>%
  dplyr::select(-data) %>%
  unnest(reduced)

## Build the MGF file
mgf_output <- list()
for(i in 1:nrow(spectra_combined)){
  
  #Extract out the feature name
  feature_name <- spectra_combined[i,]$FeatureID
  
  #Pull mz
  mz <- as.numeric(spectra_combined[i,]$`m/z`)
  #Pull rt 
  rt <- as.numeric(spectra_combined[i,]$`RT [min]`)
  #Grab the actual fragment data
  ions <- data.frame(mz = spectra_combined$mzs[[i]],
                     int = spectra_combined$intensities[[i]]) %>%
    magrittr::set_colnames(NULL)
  
  #Figure out charge
  raw_charge <- spectra_combined$attributes[[i]] %>%
    filter(name == 'scan polarity') %>%
    pull(value)
  
  #Set charge
  charge <- ifelse(raw_charge == 'positive scan', '1+', '1-')
  
  #Grab the adduct
  # adduct <- spectra_combined[i,]$`Reference Ion`
  
  #Format the fragment data
  ions_collapsed <- apply(ions, 1, function(x) paste0(x, collapse = ' '))
  
  #Build the mgf file
  mgf_format <- list(c('BEGIN IONS',
                       paste0('FEATURE_ID=', feature_name),
                       paste0('SPECTRUMID=', feature_name),
                       'MSLEVEL=2',
                       paste0('PEPMASS=', mz),
                       paste0('CHARGE=', charge),
                       paste0('TITLE=', feature_name),
                       paste0('FILENAME=', feature_name),
                       paste0('RTINSECONDS=', rt*60),
                       #paste0('ION=', adduct),
                       ions_collapsed,
                       'END IONS',
                       ''))
  #Ammend the mgf file
  mgf_output <- c(mgf_output, mgf_format)
  
}
## Turn into a vector to be written out
final_mgf <- unlist(mgf_output)
## Save the results
write_lines(final_mgf, 'output_tables/palsa_ms2_spectra.mgf')


