library(tidyverse)
library(classyfireR)
source('functions_cdis_exploration.R')

# Load data ----

cd_results_table <- read_xlsx('data/CD_results_palsa_1.5e6_old_al_node.xlsx') %>%
  # mutate(Name = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Name),
  #        Formula = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Formula)) %>% 
  mutate(MW = formatC(round(`Calc. MW`, 5), digits = 5, format = 'f'),
         RT = formatC(round(`RT [min]`, 3), digits = 3, format = 'f'),
         FeatureID = paste0('MW_', MW, '@RT_', RT),
         Name = ifelse(is.na(Name), FeatureID, Name)) %>%
  select(FeatureID, Name, Formula, MW, mz = `m/z`, contains('Area:'), 
         contains('Labeling Status:'),
         contains('Rel. Exchange'),
         contains('Peak Rating')) %>% 
  #filter(if_any(contains('Peak Rating'), ~.x > 5)) %>% 
  rename_with(function(.x){
    y <- str_remove(.x, 'Saleska_') %>% 
      str_remove('_26Jan.*')
    return(y)
  }) %>% 
  select(-contains('Peak Rating'))

# Split formula column into elemental counts
lcms_names <- cd_results_table %>% 
  select(FeatureID, Name, MW, mz, Formula)

# write_csv(lcms_names, 'output_tables/lcms_annotated_name.csv')

# The mgf file provide the order of the spectra for matching
mgf_file <- read_lines('output_tables/palsa_ms2_spectra.mgf')
mgf_names <- str_remove(mgf_file[which(str_detect(mgf_file, 'FEATURE_ID'))], 'FEATURE_ID=')

gnps_annot <- read_tsv('data/MOLECULAR-LIBRARYSEARCH-V2-a153089b-download_all_identifications-main.tsv') %>% 
  mutate(spectra_number = as.numeric(str_remove(FileScanUniqueID, 'temp.*mgf')),
         FeatureID = mgf_names[spectra_number]) %>% 
  filter(str_detect(IonMode, 'negative|Negative'),
         str_detect(Adduct, 'M-H'),
         !str_detect(Adduct, '2'))

write_csv(gnps_annot, 'output_tables/gnps_annot.csv')

canopus_data <- read_tsv('data/canopus_compound_summary.tsv') %>% 
  select(featureId, molecularFormula, contains('ClassyFire')) %>% 
  select(-contains('most specific'), -contains('all classifications')) %>% 
  rename_with(~str_remove(.x, 'ClassyFire#'))

# Adding GNPS annotations

lcms_annotation_w_gnps <- lcms_names %>% 
  left_join(select(gnps_annot, FeatureID, GNPS_Compound_Name = Compound_Name)) %>% 
  mutate(final_name = case_when(!str_detect(Name, '^MW_') ~ Name,
                                str_detect(Name, '^MW_') & !is.na(GNPS_Compound_Name) ~ GNPS_Compound_Name,
                                TRUE ~ Name),
         annot_source = case_when(!str_detect(Name, '^MW_') ~ 'Compound Discoverer',
                                  str_detect(Name, '^MW_') & !is.na(GNPS_Compound_Name) ~ 'GNPS',
                                  TRUE ~ 'No annotation'))

write_csv(lcms_annotation_w_gnps, 'output_tables/lcms_annotation_to_curate.csv')

# Curation include search and addition of InChI, SMILES and change of formula for those with GNPS annotation

  
lcms_annotation_curated <- read_csv('output_tables/lcms_annotation_curated.csv') %>% 
  mutate(Formula = str_remove_all(Formula, ' '))

# Get classification from Classyfire

inchikeys <- set_names(
  lcms_annotation_curated$InChIKey,
  nm = lcms_annotation_curated$FeatureID
)[!is.na(lcms_annotation_curated$InChIKey)]

classification_list <- map(inchikeys, get_classification)

saveRDS(classification_list, 'output_tables/classification_list.rds')

classification_list <- read_rds('output_tables/classification_list.rds')

classyfire_df <- imap(classification_list, function(cf_obj, name){
  
  temp <- tibble(Level = c('kingdom', 'superclass', 'class', 'subclass', 'level 5'))
  
  if(!is.null(cf_obj)){
    df <- classification(cf_obj) %>% 
      right_join(temp, by = 'Level') %>% 
      mutate(FeatureID = name) %>% 
      select(-CHEMONT) %>% 
      pivot_wider(names_from = Level, values_from = Classification) %>% 
      select(FeatureID, kingdom, superclass, class, subclass, `level 5`)
  } else {
    df <- tibble(FeatureID = name,
                 kingdom = NA,
                 superclass = NA,
                 class = NA,
                 subclass = NA,
                 `level 5` = NA)
  }
  return(df)
}) %>% 
  reduce(rbind) %>% 
  inner_join(select(lcms_annotation_curated, FeatureID, Formula)) %>% 
  mutate(Classification_from = 'Classyfire from InChI')

# Arranging data from SIRIUS

canopus_class <- canopus_data %>% 
  select(-contains('Probability'), -contains('probability')) %>% 
  pivot_longer(!c(featureId, molecularFormula), names_to = 'Level', values_to = 'Classification') 

canopus_probability <- canopus_data %>% 
  select(featureId, molecularFormula, contains('Probability'), contains('probability')) %>% 
  pivot_longer(!c(featureId, molecularFormula), names_to = 'Level', values_to = 'probability') %>% 
  mutate(Level = str_remove(Level, ' Probability'),
         Level = str_remove(Level, ' probability'))

canopus_final <- canopus_class %>% 
  inner_join(canopus_probability) %>% 
  filter(probability > 0.75) %>% 
  select(-probability) %>% 
  pivot_wider(names_from = Level, values_from = Classification) %>% 
  mutate(kingdom = 'Organic compounds') %>% 
  select(FeatureID = featureId, Formula = molecularFormula, kingdom, superclass, class, subclass, `level 5`) %>%
  filter(!(FeatureID %in% classyfire_df$FeatureID)) %>% 
  mutate(Classification_from = 'CANOPUS')

classes_final <- rbind(classyfire_df, canopus_final) %>% 
  rename(annot_formula = Formula)

lcms_annotation_class <- lcms_annotation_curated %>% 
  left_join(classes_final, by = c('FeatureID')) %>% 
  mutate(Formula = ifelse(Classification_from == 'CANOPUS', annot_formula, Formula)) %>% 
  select(-annot_formula)
  
lcms_annotation_final <- lcms_annotation_class %>% 
  mutate(Formula = gsub('(\\d*)', '\\1 ', Formula),
         Formula = str_trim(Formula)) %>% 
    separate_formula(.) %>% 
    calc_ratios_n_idxs(.)

write_csv(lcms_annotation_final, 'output_tables/lcms_annotation_final.csv')
