# -------------------------------------------------------------------------
plot_dwt_clusters <- function(dwt_res, nrow = NULL, ncol = NULL){

  clusts <- tibble(FeatureID = names(dwt_res@datalist),
                   cluster = paste0('Cluster ', dwt_res@cluster))
  
  feature_dwt_clusts <- as.data.frame(dwt_res@datalist) %>% 
    rowid_to_column(var = 'time') %>% 
    pivot_longer(!time, names_to = 'FeatureID', values_to = 'value') %>% 
    mutate(label_status = case_when(FeatureID %in% label_lists$T0_LO &
                                      FeatureID %in% label_inc_all ~ 'Label in PL and LO',
                                    FeatureID %in% label_lists$T0_LO ~ 'Label only in LO',
                                    FeatureID %in% label_inc_all ~ 'Label only in PL',
                                    TRUE ~ 'Unlabeled'),
           time = time -1,
           time = paste0('T', time)) %>% 
    left_join(clusts, by = 'FeatureID')
  
  rename_vec <- set_names(1:length(unique(dwt_res@cluster)),
                          paste0('Cluster ', sort(unique(dwt_res@cluster))))
  
  dwt_centroids <- as.data.frame(dwt_res@centroids) %>% 
    rename(rename_vec) %>% 
    rowid_to_column(var = 'time') %>% 
    mutate(time = time -1,
           time = paste0('T', time)) %>% 
    pivot_longer(!time, names_to = 'cluster', values_to = 'value')
  
  #label_status_color <- label_status_color[unique(feature_dwt_clusts$label_status)]
  
  feature_dwt_plot <- ggplot(feature_dwt_clusts) +
    geom_line(aes(x = time,
                  y = value,
                  color = label_status,
                  group = FeatureID),
              linewidth = 0.8) +
    geom_line(data = dwt_centroids,
              aes(y = value,
                  x = time,
                  group = 1),
              linetype = 'dashed',
              color = 'black',
              linewidth = 1) +
    facet_wrap(~cluster, nrow = nrow, ncol = ncol) +
    scale_color_manual(values = c(get_palette('d3', 3), 'gray')) +
    theme_bw() 
  
  return(feature_dwt_plot)
}