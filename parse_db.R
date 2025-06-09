library(tidyverse)
library(DBI)
library(MSnbase)
library(readxl)

db_file <- 'data/palsa/palsa_labeled/palsa_labeled_snt_3_int_2e6.db'
output_file <- 'data/palsa_labeled.mgf'

raw_files <- tibble(full_path = list.files(path = 'data/mzML', pattern = '.raw', full.names = TRUE)) %>% 
  mutate(file = basename(full_path))

compound_file <- read_xlsx('data/palsa_labeled_snt_3_int_2e6.xlsx') %>% 
  # mutate(Name = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Name),
  #        Formula = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Formula)) %>% 
  arrange(desc(`Calc. MW`)) %>% 
  mutate(FeatureID = paste0('Feature',formatC(n():0001, 
                                              width = 4, 
                                              flag = '0'))) %>% 
  select(FeatureID, Name, Formula, `Calc. MW`, `RT [min]`) %>% 
  mutate(MW = formatC(round(`Calc. MW`, 5), digits = 5, format = 'f'),
         RT = formatC(round(`RT [min]`, 3), digits = 3, format = 'f'),
         Name = ifelse(is.na(Name), paste0('MW_', MW, '@RT_', RT), Name),
         `m/z` = `Calc. MW` - 1.00784) %>% 
  select(-`Calc. MW`, -`RT [min]`)

#Open the database
msdb <- dbConnect(RSQLite::SQLite(), db_file)

#Aquire the necessary information from the database
#Compound information:
compound_info <- dbGetQuery(msdb, 'SELECT CompoundId, Formula, Name 
                            FROM CompoundTable;')
#Spectra InformationL:
spectra_info <- dbGetQuery(msdb, 'SELECT SpectrumId, CompoundId, RetentionTime, PrecursorMass, ScanNumber, RawFileURL 
                           FROM SpectrumTable;')
#Close the connection to the database:
dbDisconnect(msdb)

#Write out the results into a dataframe
total_info <- left_join(compound_info, spectra_info)

subset_compounds <- total_info


total_nest <- subset_compounds %>%
  group_by(RawFileURL) %>%
  nest() %>%
  mutate(final_name = basename(RawFileURL)) %>%
  mutate(scans = map_chr(data, function(x){
    sc <- x %>%
      drop_na() %>%
      dplyr::pull(ScanNumber) %>%
      unique() %>%
      paste0(., collapse = ',')
    return(sc) })) %>%
  inner_join(raw_files, by = c('final_name' = 'file')) %>% 
  mutate(call = paste0("../ThermoRawFileParser1.4.2/ThermoRawFileParser.exe query -i=", 
                       full_path,
                       " -n=\"",
                       scans,
                       "\" -s"))

spectra <- map_df(total_nest$call, function(x){
  print(x)
  specframe <- jsonlite::fromJSON(system(x, intern = T)) %>%
    mutate(ScanNumber = map_chr(attributes, function(x) x %>% filter(name == 'scan number') %>% pull(value))) %>%
    modify_at('ScanNumber', as.integer)
})


spectra_linked <- spectra  %>% 
  left_join(subset_compounds, relationship = 'many-to-many') %>%
  left_join(compound_file, by = c('Name', 'Formula'), relationship = 'many-to-many')

spectra_combined <- spectra_linked %>%
  group_by(CompoundId) %>%
  mutate(med_pepmass = median(PrecursorMass),
         med_rt = median(RetentionTime)) %>% 
  nest() %>% 
  mutate(spec_list = purrr::map(data, function(s) {
    map2(s$mzs, s$intensities, function(x,y) new('Spectrum2', mz = x, intensity = y))}))  %>%
  mutate(m_spec = purrr::map(spec_list, MSpectra)) %>%
  mutate(con_spec = purrr::map(m_spec, function(x){
    mz <- numeric(0)
    prop <- 0.5
    while(identical(mz, numeric(0)) && prop > 0.1){
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

#Step 5: Build the MGF file
mgf_output <- list()
for(i in 1:nrow(spectra_combined)){
  
  #Extract out the feature name
  feature_name <- spectra_combined[i,]$FeatureID
  
  #Pull mz
  mz <- spectra_combined[i,]$med_pepmass
  #Pull rt 
  rt <- as.numeric(spectra_combined[i,]$med_rt)
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
#Turn into a vector to be written out
final_mgf <- unlist(mgf_output)
#Writeout the results
final_file <- file(output_file)
writeLines(final_mgf, final_file)
close(final_file)

