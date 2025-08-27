
# drivers on bird species

gam_bird_PLS <- function(bird_data,pressure_data,
                             pressure_name = c("d_impervious","d_treedensity","d_agri",
                                               "d_tempsrping","tempsrping","d_tempsrpingvar","d_precspring","precspring",
                                               "d_shannon","shannon","drymatter","protectedarea_perc",
                                               "eulandsystem_farmland_low","eulandsystem_farmland_medium","eulandsystem_farmland_high",
                                               "eulandsystem_forest_lowmedium","eulandsystem_forest_high","milieu_cat"),
                             min_site_number_per_species = 60,
                             min_occurence_species=200,
                             family="quasipoisson"){
  
  species_press_data_year <- merge(bird_data, pressure_data[which(pressure_data$siteID %in% unique(bird_data$siteID) & pressure_data$year %in% unique(bird_data$year)),], by =c("siteID","year"), all.x=TRUE)
  
  poisson_df <- na.omit(species_press_data_year[,c("siteID","count","year","area_sampled_m2","scheme_code","Long_LAEA","Lat_LAEA",
                                                   pressure_name,"PLS")])
  
  poisson_df$year <- scale(poisson_df$year)
  
  if(length(table(poisson_df$area_sampled_m2)) > length(unique(poisson_df$scheme_code))){
    one_scheme_time_area <- 0 
    poisson_df$area_sampled_m2 <- scale(poisson_df$area_sampled_m2)
  }else{
    one_scheme_time_area <- 1
  }
  
  poisson_df$count_scale_all <- poisson_df$count
  
  if(length(pressure_name) > 1){
    formula_gam <- "count_scale_all ~ year + year:d_impervious + year:d_treedensity +
    year:eulandsystem_forest_lowmedium + year:eulandsystem_forest_high +
    year:d_agri + year:eulandsystem_farmland_low + year:eulandsystem_farmland_medium + year:eulandsystem_farmland_high +
    year:d_tempsrping + year:d_tempsrpingvar + year:d_precspring + year:d_shannon + year:protectedarea_perc +
    milieu_cat + tempsrping + precspring + shannon + drymatter"
  }else{
    formula_gam <- paste("count_scale_all ~", paste(pressure_name,sep="", collapse = " + "))
  }
  
  col_names <- c("(Intercept)","year","milieu_catopenland","milieu_catothers","milieu_caturban",
                 "tempsrping","precspring","shannon","drymatter","year:d_impervious","year:d_treedensity",
                 "year:eulandsystem_forest_lowmedium","year:eulandsystem_forest_high","year:d_agri",
                 "year:eulandsystem_farmland_low","year:eulandsystem_farmland_medium","year:eulandsystem_farmland_high",
                 "year:d_tempsrping","year:d_tempsrpingvar","year:d_precspring","year:d_shannon","year:protectedarea_perc")
  
  if(nrow(poisson_df) >= min_occurence_species){
    
    ### global poisson model (gamm too resource consumming over the whole Europe)
    
    if(length(unique(poisson_df$scheme_code)) > 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2:scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) == 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) > 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) == 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    
   
    if(global_mod$converged){
      
      unique_poisson_df <- distinct(poisson_df, Long_LAEA, Lat_LAEA,.keep_all = TRUE)
      
      if_fail <- rbind(matrix(NA,(nrow=length(col_names)+1),ncol=4),c(0,rep(0,3)))
      row.names(if_fail) <- c(col_names,"dev_exp","n_obs")
      
      result_all_site <- daply(unique_poisson_df,.(PLS),.fun=purrr::possibly(otherwise=if_fail,
                                                                             .f=function(x,min_site_number_per_species,poisson_df){
                                                                               
                                                                               if(nrow(x) >= min_site_number_per_species){
                                                                                 
                                                                                 poisson_df_i <- poisson_df[which(poisson_df$PLS == unique(x$PLS)),]

                                                                                 if(length(table(poisson_df_i$area_sampled_m2)) > length(unique(poisson_df_i$scheme_code))){
                                                                                   one_scheme_time_area <- 0 
                                                                                   poisson_df_i$area_sampled_m2 <- scale(poisson_df_i$area_sampled_m2)
                                                                                 }else{
                                                                                   one_scheme_time_area <- 1
                                                                                 }
                                                                                 
                                                                                 if(length(unique(poisson_df_i$scheme_code)) > 1 && one_scheme_time_area == 0){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2:scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("scheme_code|area_sampled_m2",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$scheme_code)) == 1 && one_scheme_time_area == 0){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("area_sampled_m2",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$scheme_code)) > 1 && one_scheme_time_area == 1){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("scheme_code",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$scheme_code)) == 1 && one_scheme_time_area == 1){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("no_",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 
                                                                                 if(nrow(result_i) == length(col_names)){
                                                                                   result_site <- result_i
                                                                                 }else{
                                                                                   row_to_add <- matrix(NA,nrow=length(which(!(col_names %in% row.names(result_i)))), ncol=1)
                                                                                   row.names(row_to_add) <- col_names[which(!(col_names %in% row.names(result_i)))]
                                                                                   result_i_complet <- merge(result_i,row_to_add,by="row.names",all=TRUE)
                                                                                   result_i_complet <- result_i_complet[match(col_names, result_i_complet$Row.names),]
                                                                                   result_i_complet <- as.matrix(result_i_complet[2:5])
                                                                                   result_site <- result_i_complet
                                                                                 }
                                                                                 result_site <- rbind(result_site,c(dev_exp,rep(0,3)))
                                                                                 result_site <- rbind(result_site,c(n_obs,rep(0,3)))
                                                                                 
                                                                               }else{
                                                                                 n_obs <- nrow(poisson_df[which(poisson_df$PLS == unique(x$PLS)),])
                                                                                 result_site <- matrix(NA,(nrow=length(col_names)+1),ncol=4)
                                                                                 result_site <- rbind(result_site,c(n_obs,rep(0,3)))
                                                                               }
                                                                               
                                                                               row.names(result_site) <- c(col_names,"dev_exp","n_obs")
                                                                               
                                                                               return(result_site)
                                                                             }),
                               min_site_number_per_species=min_site_number_per_species,poisson_df=poisson_df,
                               .progress="none")
      
      if(!is.na(dim(result_all_site)[3])){
        result_all_site <- aperm(result_all_site, c(2,3,1))
        
        if(dim(result_all_site)[3] > 1){
          res.poisson_df <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),1,])))
          res.poisson_sd <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)),2,])))
          res.poisson_pval <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),4,])))
        }
        if(dim(result_all_site)[3] == 1){
          res.poisson_df <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),1])))
          res.poisson_sd <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)),2])))
          res.poisson_pval <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),4])))
        }
      }
      
      if(is.na(dim(result_all_site)[3])){
        res.poisson_df <- data.frame(matrix(NA,nrow=1,ncol=(length(col_names)+2)))
        res.poisson_sd <- data.frame(matrix(NA,nrow=1,ncol=(length(col_names))))
        res.poisson_pval <- matrix(1,nrow=1,ncol=(length(col_names)+2))
      }
      
      
      
      res.poisson_df[res.poisson_pval > 0.05] <- NA
      res.poisson_sd[res.poisson_pval[,1:length(col_names)] > 0.05] <- NA 
      
      if(is.na(dim(result_all_site)[3])){
        
        res.poisson_df$PLS <- NA
        
      }else{
        
        res.poisson_df$PLS <- gsub("X","",row.names(res.poisson_df))
        
      }
      
      names(res.poisson_sd) <- c(paste0(col_names,"_sd"))
      
      res.poisson_df <- cbind(res.poisson_df,res.poisson_sd)
      
      
      
      global_mod_coef <- summary(global_mod)$p.table[grep("scheme_code|area_sampled_m2|time_effort|no_",row.names(summary(global_mod)$p.table),invert = TRUE),]
      
      if(nrow(global_mod_coef) < length(col_names)){
        row_to_add <- matrix(NA,nrow=length(which(!(col_names %in% row.names(global_mod_coef)))), ncol=1)
        row.names(row_to_add) <- col_names[which(!(col_names %in% row.names(global_mod_coef)))]
        global_mod_coef_complet <- merge(global_mod_coef,row_to_add,by="row.names",all=TRUE)
        global_mod_coef_complet <- global_mod_coef_complet[match(col_names, global_mod_coef_complet$Row.names),]
        global_mod_coef_complet <- as.matrix(global_mod_coef_complet[2:5])
        global_mod_coef <- global_mod_coef_complet
      }
      
      global_mod_coef <- rbind(global_mod_coef,c(summary(global_mod)$dev.expl,rep(0,3)),c(summary(global_mod)$n,rep(0,3)))
      
      global_mod_coef1 <- global_mod_coef[,1]
      global_mod_coef1[which(global_mod_coef[,4] > 0.05)] <- NA
      global_mod_df <- data.frame(t(global_mod_coef1))
      names(global_mod_df) <- c(col_names,"dev_exp","n_obs")
      global_mod_df$PLS <- "europe"
      
      global_mod_coef2 <- global_mod_coef[1:length(col_names),2]
      global_mod_coef2[which(global_mod_coef[,4] > 0.05)] <- NA
      global_mod_sd <- data.frame(t(global_mod_coef2))
      names(global_mod_sd) <- c(paste0(col_names,"_sd"))
      
      global_mod_df <- cbind(global_mod_df,global_mod_sd)
      
      res.poisson_df <- rbind(res.poisson_df,global_mod_df)
      
    }else{
      res.poisson_df <- data.frame(t(rep(NA,(2*length(col_names)+3))))
      names(res.poisson_df) <- c(col_names,"dev_exp","n_obs","PLS",paste0(col_names,"_sd"))
    }
    
  }else{
    res.poisson_df <- data.frame(t(rep(NA,(2*length(col_names)+3))))
    names(res.poisson_df) <- c(col_names,"dev_exp","n_obs","PLS",paste0(col_names,"_sd"))
  }
  
  return(res.poisson_df)
}


# drivers on butterfly species


gam_butterfly_PLS <- function(butterfly_data,pressure_data,
                              pressure_name = c("d_impervious","d_treedensity","d_agri",
                                                "d_tempsrping","tempsrping","d_tempsrpingvar","d_precspring","precspring",
                                                "d_shannon","shannon","drymatter","protectedarea_perc",
                                                "eulandsystem_farmland_low","eulandsystem_farmland_medium","eulandsystem_farmland_high",
                                                "eulandsystem_forest_lowmedium","eulandsystem_forest_high","milieu_cat"),
                              min_site_number_per_species = 60,
                              min_occurence_species=200,
                              family="quasipoisson"){
  
  species_press_data_year <- merge(butterfly_data, pressure_data[which(pressure_data$transect_id %in% unique(butterfly_data$transect_id) & pressure_data$year %in% unique(butterfly_data$year)),], by =c("transect_id","year"), all.x=TRUE)
  
  poisson_df <- na.omit(species_press_data_year[,c("transect_id","count_corrected","year","transect_length","bms_id","Long_LAEA","Lat_LAEA",
                                                   pressure_name,"PLS")])
  
  poisson_df$year <- scale(poisson_df$year)
  
  if(length(table(poisson_df$transect_length)) > length(unique(poisson_df$bms_id))){
    one_scheme_time_area <- 0 
    poisson_df$transect_length <- scale(poisson_df$transect_length)
  }else{
    one_scheme_time_area <- 1
  }
  
  poisson_df$count_corrected_scale_all <- poisson_df$count_corrected
  
  if(length(pressure_name) > 1){
    formula_gam <- "count_corrected_scale_all ~ year + year:d_impervious + year:d_treedensity +
    year:eulandsystem_forest_lowmedium + year:eulandsystem_forest_high +
    year:d_agri + year:eulandsystem_farmland_low + year:eulandsystem_farmland_medium + year:eulandsystem_farmland_high +
    year:d_tempsrping + year:d_tempsrpingvar + year:d_precspring + year:d_shannon + year:protectedarea_perc +
    milieu_cat + tempsrping + precspring + shannon + drymatter"
  }else{
    formula_gam <- paste("count_corrected_scale_all ~", paste(pressure_name,sep="", collapse = " + "))
  }
  
  col_names <- c("(Intercept)","year","milieu_catopenland","milieu_catothers","milieu_caturban",
                 "tempsrping","precspring","shannon","drymatter","year:d_impervious","year:d_treedensity",
                 "year:eulandsystem_forest_lowmedium","year:eulandsystem_forest_high","year:d_agri",
                 "year:eulandsystem_farmland_low","year:eulandsystem_farmland_medium","year:eulandsystem_farmland_high",
                 "year:d_tempsrping","year:d_tempsrpingvar","year:d_precspring","year:d_shannon","year:protectedarea_perc")
  
  if(nrow(poisson_df) >= min_occurence_species){
    
    ### global poisson model (gamm too resource consumming over the whole Europe)
    
    if(length(unique(poisson_df$bms_id)) > 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length:bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) == 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) > 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) == 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    
    
    if(global_mod$converged){
      
      unique_poisson_df <- distinct(poisson_df, Long_LAEA, Lat_LAEA,.keep_all = TRUE)
      
      if_fail <- rbind(matrix(NA,(nrow=length(col_names)+1),ncol=4),c(0,rep(0,3)))
      row.names(if_fail) <- c(col_names,"dev_exp","n_obs")
      
      result_all_site <- daply(unique_poisson_df,.(PLS),.fun=purrr::possibly(otherwise=if_fail,
                                                                             .f=function(x,min_site_number_per_species,poisson_df){
                                                                               
                                                                               if(nrow(x) >= min_site_number_per_species){
                                                                                 
                                                                                 poisson_df_i <- poisson_df[which(poisson_df$PLS == unique(x$PLS)),]

                                                                                 if(length(table(poisson_df_i$transect_length)) > length(unique(poisson_df_i$bms_id))){
                                                                                   one_scheme_time_area <- 0 
                                                                                   poisson_df_i$transect_length <- scale(poisson_df_i$transect_length)
                                                                                 }else{
                                                                                   one_scheme_time_area <- 1
                                                                                 }
                                                                                 
                                                                                 if(length(unique(poisson_df_i$bms_id)) > 1 && one_scheme_time_area == 0){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length:bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("bms_id|transect_length",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$bms_id)) == 1 && one_scheme_time_area == 0){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("transect_length",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$bms_id)) > 1 && one_scheme_time_area == 1){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("bms_id",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 if(length(unique(poisson_df_i$bms_id)) == 1 && one_scheme_time_area == 1){
                                                                                   res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                         family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   result_i <- summary(res.poisson_i$gam)$p.table
                                                                                   dev_exp <- summary(res.poisson_i$gam)$r.sq
                                                                                   n_obs <- summary(res.poisson_i$gam)$n
                                                                                   result_i <- as.matrix(result_i[grep("no_",row.names(result_i),invert = TRUE),])
                                                                                 }
                                                                                 
                                                                                 if(nrow(result_i) == length(col_names)){
                                                                                   result_site <- result_i
                                                                                 }else{
                                                                                   row_to_add <- matrix(NA,nrow=length(which(!(col_names %in% row.names(result_i)))), ncol=1)
                                                                                   row.names(row_to_add) <- col_names[which(!(col_names %in% row.names(result_i)))]
                                                                                   result_i_complet <- merge(result_i,row_to_add,by="row.names",all=TRUE)
                                                                                   result_i_complet <- result_i_complet[match(col_names, result_i_complet$Row.names),]
                                                                                   result_i_complet <- as.matrix(result_i_complet[2:5])
                                                                                   result_site <- result_i_complet
                                                                                 }
                                                                                 result_site <- rbind(result_site,c(dev_exp,rep(0,3)))
                                                                                 result_site <- rbind(result_site,c(n_obs,rep(0,3)))
                                                                                 
                                                                               }else{
                                                                                 n_obs <- nrow(poisson_df[which(poisson_df$PLS == unique(x$PLS)),])
                                                                                 result_site <- matrix(NA,(nrow=length(col_names)+1),ncol=4)
                                                                                 result_site <- rbind(result_site,c(n_obs,rep(0,3)))
                                                                               }
                                                                               
                                                                               row.names(result_site) <- c(col_names,"dev_exp","n_obs")
                                                                               
                                                                               return(result_site)
                                                                             }),
                               min_site_number_per_species=min_site_number_per_species,poisson_df=poisson_df,
                               .progress="none")
      
      if(!is.na(dim(result_all_site)[3])){
        result_all_site <- aperm(result_all_site, c(2,3,1))
        
        if(dim(result_all_site)[3] > 1){
          res.poisson_df <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),1,])))
          res.poisson_sd <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)),2,])))
          res.poisson_pval <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),4,])))
        }
        if(dim(result_all_site)[3] == 1){
          res.poisson_df <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),1])))
          res.poisson_sd <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)),2])))
          res.poisson_pval <- as.data.frame(t(data.frame(result_all_site[1:(length(col_names)+2),4])))
        }
      }
      
      if(is.na(dim(result_all_site)[3])){
        res.poisson_df <- data.frame(matrix(NA,nrow=1,ncol=(length(col_names)+2)))
        res.poisson_sd <- data.frame(matrix(NA,nrow=1,ncol=(length(col_names))))
        res.poisson_pval <- matrix(1,nrow=1,ncol=(length(col_names)+2))
      }
      
      
      
      res.poisson_df[res.poisson_pval > 0.05] <- NA
      res.poisson_sd[res.poisson_pval[,1:length(col_names)] > 0.05] <- NA
      
      if(is.na(dim(result_all_site)[3])){
        
        res.poisson_df$PLS <- NA
        
      }else{
        
        res.poisson_df$PLS <- gsub("X","",row.names(res.poisson_df))
        
      }
      
      names(res.poisson_sd) <- c(paste0(col_names,"_sd"))
      
      res.poisson_df <- cbind(res.poisson_df,res.poisson_sd)
      
      global_mod_coef <- summary(global_mod)$p.table[grep("bms_id|transect_length|time_effort|no_",row.names(summary(global_mod)$p.table),invert = TRUE),]
      
      if(nrow(global_mod_coef) < length(col_names)){
        row_to_add <- matrix(NA,nrow=length(which(!(col_names %in% row.names(global_mod_coef)))), ncol=1)
        row.names(row_to_add) <- col_names[which(!(col_names %in% row.names(global_mod_coef)))]
        global_mod_coef_complet <- merge(global_mod_coef,row_to_add,by="row.names",all=TRUE)
        global_mod_coef_complet <- global_mod_coef_complet[match(col_names, global_mod_coef_complet$Row.names),]
        global_mod_coef_complet <- as.matrix(global_mod_coef_complet[2:5])
        global_mod_coef <- global_mod_coef_complet
      }
      
      global_mod_coef <- rbind(global_mod_coef,c(summary(global_mod)$dev.expl,rep(0,3)),c(summary(global_mod)$n,rep(0,3)))
      
      global_mod_coef1 <- global_mod_coef[,1]
      global_mod_coef1[which(global_mod_coef[,4] > 0.05)] <- NA
      global_mod_df <- data.frame(t(global_mod_coef1))
      names(global_mod_df) <- c(col_names,"dev_exp","n_obs")
      global_mod_df$PLS <- "europe"
      
      global_mod_coef2 <- global_mod_coef[1:length(col_names),2]
      global_mod_coef2[which(global_mod_coef[,4] > 0.05)] <- NA
      global_mod_sd <- data.frame(t(global_mod_coef2))
      names(global_mod_sd) <- c(paste0(col_names,"_sd"))
      
      global_mod_df <- cbind(global_mod_df,global_mod_sd)
      
      res.poisson_df <- rbind(res.poisson_df,global_mod_df)
      
    }else{
      res.poisson_df <- data.frame(t(rep(NA,(2*length(col_names)+3))))
      names(res.poisson_df) <- c(col_names,"dev_exp","n_obs","PLS",paste0(col_names,"_sd"))
    }
    
  }else{
    res.poisson_df <- data.frame(t(rep(NA,(2*length(col_names)+3))))
    names(res.poisson_df) <- c(col_names,"dev_exp","n_obs","PLS",paste0(col_names,"_sd"))
  }
  
  return(res.poisson_df)
}

# trend in each scenario (core function)

predict_trend <- function(mod,
                          pressure_data_unscale,
                          poisson_df_unscale,
                          poisson_df,
                          lulc_pls_short,
                          climate_pls,
                          pa_pls_short,
                          PLS,
                          nb_rep=1000,
                          pressure_remove = NULL){
  
  mod_coef <- summary(mod)$p.table[grep("year",row.names(summary(mod)$p.table)),]
  
  if(!is.null(pressure_remove)){
    mod_coef[pressure_remove,c("Estimate","Std. Error")] <- 0
    pressure_removed <- c("year","d_impervious","d_tempsrping","d_tempsrpingvar","d_precspring",
                          "d_shannon","protectedarea_perc","d_treedensity","eulandsystem_forest_lowmedium","eulandsystem_forest_high",
                          "d_agri","eulandsystem_farmland_low","eulandsystem_farmland_medium","eulandsystem_farmland_high")[pressure_remove]
  }else{pressure_removed <- "none"}
  
  year_si <- sd(na.omit(pressure_data_unscale$year))
  d_impervious_si <- sd(na.omit(pressure_data_unscale$d_impervious))
  d_tempspring_si <- sd(na.omit(pressure_data_unscale$d_tempsrping))
  d_tempspringvar_si <- sd(na.omit(pressure_data_unscale$d_tempsrpingvar))
  d_precspring_si <- sd(na.omit(pressure_data_unscale$d_precspring))
  d_shannon_si <- sd(na.omit(pressure_data_unscale$d_shannon))
  protectedarea_perc_si <- sd(na.omit(pressure_data_unscale$protectedarea_perc))
  d_treedensity_si <- sd(na.omit(pressure_data_unscale$d_treedensity))
  d_agri_si <- sd(na.omit(pressure_data_unscale$d_agri))
  eulandsystem_farmland_low_si <- sd(na.omit(pressure_data_unscale$eulandsystem_farmland_low))
  eulandsystem_farmland_medium_si <- sd(na.omit(pressure_data_unscale$eulandsystem_farmland_medium))
  eulandsystem_farmland_high_si <- sd(na.omit(pressure_data_unscale$eulandsystem_farmland_high))
  eulandsystem_forest_lowmedium_si <- sd(na.omit(pressure_data_unscale$eulandsystem_forest_lowmedium))
  eulandsystem_forest_high_si <- sd(na.omit(pressure_data_unscale$eulandsystem_forest_high))
  
  year_mu <- mean(na.omit(pressure_data_unscale$year))
  d_impervious_mu <- mean(na.omit(pressure_data_unscale$d_impervious))
  d_tempspring_mu <- mean(na.omit(pressure_data_unscale$d_tempsrping))
  d_tempspringvar_mu <- mean(na.omit(pressure_data_unscale$d_tempsrpingvar))
  d_precspring_mu <- mean(na.omit(pressure_data_unscale$d_precspring))
  d_shannon_mu <- mean(na.omit(pressure_data_unscale$d_shannon))
  protectedarea_perc_mu <- mean(na.omit(pressure_data_unscale$protectedarea_perc))
  d_treedensity_mu <- mean(na.omit(pressure_data_unscale$d_treedensity))
  d_agri_mu <- mean(na.omit(pressure_data_unscale$d_agri))
  eulandsystem_farmland_low_mu <- mean(na.omit(pressure_data_unscale$eulandsystem_farmland_low))
  eulandsystem_farmland_medium_mu <- mean(na.omit(pressure_data_unscale$eulandsystem_farmland_medium))
  eulandsystem_farmland_high_mu <- mean(na.omit(pressure_data_unscale$eulandsystem_farmland_high))
  eulandsystem_forest_lowmedium_mu <- mean(na.omit(pressure_data_unscale$eulandsystem_forest_lowmedium))
  eulandsystem_forest_high_mu <- mean(na.omit(pressure_data_unscale$eulandsystem_forest_high))
  
  nb_rep <- 1000
  
  beta1_past <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*mean(poisson_df$d_impervious) +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*mean(poisson_df$d_tempsrping) +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*mean(poisson_df$d_tempsrpingvar) +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*mean(poisson_df$d_precspring) +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*mean(poisson_df$d_shannon) +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*mean(poisson_df$protectedarea_perc) +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*mean(poisson_df$d_treedensity) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*mean(poisson_df$eulandsystem_forest_lowmedium) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*mean(poisson_df$eulandsystem_forest_high) +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*mean(poisson_df$d_agri) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_low) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_medium) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_high)
  
  beta1_past_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*mean(poisson_df$d_impervious) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*mean(poisson_df$d_tempsrping) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*mean(poisson_df$d_tempsrpingvar) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*mean(poisson_df$d_precspring) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*mean(poisson_df$d_shannon) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*mean(poisson_df$protectedarea_perc) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*mean(poisson_df$d_treedensity) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_lowmedium) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_high) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*mean(poisson_df$d_agri) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_low) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_medium) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_high)
  
  
  d_impervious_ssp1 <- mean(poisson_df_unscale$impervious_2018)*
    (lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]-1)/(2050-2018) 
  d_shannon_ssp1 <- mean(poisson_df_unscale$shannon_2018)*
    (lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]-1)/(2050-2018) 
  d_tempspring_ssp1 <- mean(poisson_df_unscale$tempspring_2020)*
    (climate_pls$mean_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$mean_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_tempspringvar_ssp1 <- mean(poisson_df_unscale$tempspringvar_2020)*
    (climate_pls$var_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$var_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_precspring_ssp1 <- mean(poisson_df_unscale$precspring_2020)*
    (climate_pls$sum_p_4_5[which(climate_pls$PLS==PLS)]/climate_pls$sum_p_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_agri_ssp1 <- mean(poisson_df_unscale$agri_2018)*
    (sum(lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])-1)/(2050-2018)
  agri_low_ssp1 <- mean(poisson_df_unscale$eulandsystem_farmland_low)*lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]
  agri_medium_ssp1 <- mean(poisson_df_unscale$eulandsystem_farmland_medium)*lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]
  agri_high_ssp1 <- mean(poisson_df_unscale$eulandsystem_farmland_high)*lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]
  d_treedensity_ssp1 <- mean(poisson_df_unscale$treedensity_2018)*
    (sum(lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])-1)/(2050-2018)
  forest_lowmedium_ssp1 <- mean(poisson_df_unscale$eulandsystem_forest_lowmedium)*lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]
  forest_high_ssp1 <- mean(poisson_df_unscale$eulandsystem_forest_high)*lulc_pls_short$ssp1[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]
  protectedarea_perc_ssp1 <- mean(poisson_df_unscale$protectedarea_perc)*pa_pls_short$ssp1[which(pa_pls_short$PLS==PLS)]/pa_pls_short$initial[which(pa_pls_short$PLS==PLS)]
  
  beta1_SSP1 <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*(d_impervious_ssp1 - d_impervious_mu)/d_impervious_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*(d_shannon_ssp1 - d_shannon_mu)/d_shannon_si +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*(d_treedensity_ssp1 - d_treedensity_mu)/d_treedensity_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_ssp1 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_ssp1 - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*(d_agri_ssp1 - d_agri_mu)/d_agri_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_ssp1 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_ssp1 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_ssp1 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_SSP1_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*(d_impervious_ssp1 - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*(d_shannon_ssp1 - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_ssp1 - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_ssp1 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_ssp1 +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*(d_agri_ssp1 - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_ssp1 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_ssp1 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_ssp1 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  
  beta1_BAU <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*mean(poisson_df$d_impervious) +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*mean(poisson_df$d_shannon) +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*mean(poisson_df$d_treedensity) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*mean(poisson_df$eulandsystem_forest_lowmedium) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*mean(poisson_df$eulandsystem_forest_high) +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*mean(poisson_df$d_agri) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_low) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_medium) +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_high)
  
  beta1_BAU_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*mean(poisson_df$d_impervious) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*mean(poisson_df$d_shannon) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*mean(poisson_df$d_treedensity) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_lowmedium) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_high) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*mean(poisson_df$d_agri) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_low) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_medium) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_high)
  
  
  
  d_impervious_ssp3 <- mean(poisson_df_unscale$impervious_2018)*
    (lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]-1)/(2050-2018) 
  d_shannon_ssp3 <- mean(poisson_df_unscale$shannon_2018)*
    (lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]-1)/(2050-2018) 
  d_tempspring_ssp3 <- mean(poisson_df_unscale$tempspring_2020)*
    (climate_pls$mean_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$mean_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_tempspringvar_ssp3 <- mean(poisson_df_unscale$tempspringvar_2020)*
    (climate_pls$var_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$var_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_precspring_ssp3 <- mean(poisson_df_unscale$precspring_2020)*
    (climate_pls$sum_p_4_5[which(climate_pls$PLS==PLS)]/climate_pls$sum_p_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_agri_ssp3 <- mean(poisson_df_unscale$agri_2018)*
    (sum(lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])-1)/(2050-2018)
  agri_low_ssp3 <- mean(poisson_df_unscale$eulandsystem_farmland_low)*lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]
  agri_medium_ssp3 <- mean(poisson_df_unscale$eulandsystem_farmland_medium)*lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]
  agri_high_ssp3 <- mean(poisson_df_unscale$eulandsystem_farmland_high)*lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]
  d_treedensity_ssp3 <- mean(poisson_df_unscale$treedensity_2018)*
    (sum(lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])-1)/(2050-2018)
  forest_lowmedium_ssp3 <- mean(poisson_df_unscale$eulandsystem_forest_lowmedium)*lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]
  forest_high_ssp3 <- mean(poisson_df_unscale$eulandsystem_forest_high)*lulc_pls_short$ssp3[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]
  protectedarea_perc_ssp3 <- mean(poisson_df_unscale$protectedarea_perc)*pa_pls_short$ssp3[which(pa_pls_short$PLS==PLS)]/pa_pls_short$initial[which(pa_pls_short$PLS==PLS)]
  
  beta1_SSP3 <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*(d_impervious_ssp3 - d_impervious_mu)/d_impervious_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp3 - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp3 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp3 - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*(d_shannon_ssp3 - d_shannon_mu)/d_shannon_si +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp3 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*(d_treedensity_ssp3 - d_treedensity_mu)/d_treedensity_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_ssp3 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_ssp3 - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*(d_agri_ssp3 - d_agri_mu)/d_agri_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_ssp3 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_ssp3 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_ssp3 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_SSP3_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*(d_impervious_ssp3 - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp3 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp3 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp3 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*(d_shannon_ssp3 - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp3 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_ssp3 - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_ssp3 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_ssp3 +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*(d_agri_ssp3 - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_ssp3 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_ssp3 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_ssp3 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  d_impervious_nac <- mean(poisson_df_unscale$impervious_2018)*
    (lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]-1)/(2050-2018) 
  d_shannon_nac <- mean(poisson_df_unscale$shannon_2018)*
    (lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]-1)/(2050-2018) 
  d_tempspring_nac <- mean(poisson_df_unscale$tempspring_2020)*
    (climate_pls$mean_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$mean_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_tempspringvar_nac <- mean(poisson_df_unscale$tempspringvar_2020)*
    (climate_pls$var_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$var_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_precspring_nac <- mean(poisson_df_unscale$precspring_2020)*
    (climate_pls$sum_p_4_5[which(climate_pls$PLS==PLS)]/climate_pls$sum_p_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_agri_nac <- mean(poisson_df_unscale$agri_2018)*
    (sum(lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])-1)/(2050-2018)
  agri_low_nac <- mean(poisson_df_unscale$eulandsystem_farmland_low)*lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]
  agri_medium_nac <- mean(poisson_df_unscale$eulandsystem_farmland_medium)*lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]
  agri_high_nac <- mean(poisson_df_unscale$eulandsystem_farmland_high)*lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]
  d_treedensity_nac <- mean(poisson_df_unscale$treedensity_2018)*
    (sum(lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])-1)/(2050-2018)
  forest_lowmedium_nac <- mean(poisson_df_unscale$eulandsystem_forest_lowmedium)*lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]
  forest_high_nac <- mean(poisson_df_unscale$eulandsystem_forest_high)*lulc_pls_short$nac[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]
  protectedarea_perc_nac <- mean(poisson_df_unscale$protectedarea_perc)*pa_pls_short$nac[which(pa_pls_short$PLS==PLS)]/pa_pls_short$initial[which(pa_pls_short$PLS==PLS)]
  
  beta1_nac <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*(d_impervious_nac - d_impervious_mu)/d_impervious_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nac - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nac - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_nac - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*(d_shannon_nac - d_shannon_mu)/d_shannon_si +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nac - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nac - d_treedensity_mu)/d_treedensity_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nac - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nac - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*(d_agri_nac - d_agri_mu)/d_agri_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nac - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nac - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nac - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nac_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*(d_impervious_nac - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nac - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nac - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_nac - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*(d_shannon_nac - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nac - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nac - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nac - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nac +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*(d_agri_nac - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nac - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nac - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nac - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  d_impervious_nfn <- mean(poisson_df_unscale$impervious_2018)*
    (lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]-1)/(2050-2018) 
  d_shannon_nfn <- mean(poisson_df_unscale$shannon_2018)*
    (lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]-1)/(2050-2018) 
  d_tempspring_nfn <- mean(poisson_df_unscale$tempspring_2020)*
    (climate_pls$mean_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$mean_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_tempspringvar_nfn <- mean(poisson_df_unscale$tempspringvar_2020)*
    (climate_pls$var_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$var_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_precspring_nfn <- mean(poisson_df_unscale$precspring_2020)*
    (climate_pls$sum_p_4_5[which(climate_pls$PLS==PLS)]/climate_pls$sum_p_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_agri_nfn <- mean(poisson_df_unscale$agri_2018)*
    (sum(lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])-1)/(2050-2018)
  agri_low_nfn <- mean(poisson_df_unscale$eulandsystem_farmland_low)*lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]
  agri_medium_nfn <- mean(poisson_df_unscale$eulandsystem_farmland_medium)*lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]
  agri_high_nfn <- mean(poisson_df_unscale$eulandsystem_farmland_high)*lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]
  d_treedensity_nfn <- mean(poisson_df_unscale$treedensity_2018)*
    (sum(lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])-1)/(2050-2018)
  forest_lowmedium_nfn <- mean(poisson_df_unscale$eulandsystem_forest_lowmedium)*lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]
  forest_high_nfn <- mean(poisson_df_unscale$eulandsystem_forest_high)*lulc_pls_short$nfn[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]
  protectedarea_perc_nfn <- mean(poisson_df_unscale$protectedarea_perc)*pa_pls_short$nfn[which(pa_pls_short$PLS==PLS)]/pa_pls_short$initial[which(pa_pls_short$PLS==PLS)]
  
  beta1_nfn <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*(d_impervious_nfn - d_impervious_mu)/d_impervious_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nfn - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nfn - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_nfn - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*(d_shannon_nfn - d_shannon_mu)/d_shannon_si +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nfn - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nfn - d_treedensity_mu)/d_treedensity_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nfn - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nfn - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*(d_agri_nfn - d_agri_mu)/d_agri_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nfn - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nfn - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nfn - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfn_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*(d_impervious_nfn - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nfn - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nfn - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_nfn - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*(d_shannon_nfn - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nfn - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nfn - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nfn - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nfn +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*(d_agri_nfn - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nfn - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nfn - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nfn - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  d_impervious_nfs <- mean(poisson_df_unscale$impervious_2018)*
    (lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("urban"))]-1)/(2050-2018) 
  d_shannon_nfs <- mean(poisson_df_unscale$shannon_2018)*
    (lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("landscape_div"))]-1)/(2050-2018) 
  d_tempspring_nfs <- mean(poisson_df_unscale$tempspring_2020)*
    (climate_pls$mean_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$mean_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_tempspringvar_nfs <- mean(poisson_df_unscale$tempspringvar_2020)*
    (climate_pls$var_t_4_5[which(climate_pls$PLS==PLS)]/climate_pls$var_t_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_precspring_nfs <- mean(poisson_df_unscale$precspring_2020)*
    (climate_pls$sum_p_4_5[which(climate_pls$PLS==PLS)]/climate_pls$sum_p_2016[which(climate_pls$PLS==PLS)]-1)/(2050-2018)
  d_agri_nfs <- mean(poisson_df_unscale$agri_2018)*
    (sum(lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low","farmland_medium","farmland_high"))])-1)/(2050-2018)
  agri_low_nfs <- mean(poisson_df_unscale$eulandsystem_farmland_low)*lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_low"))]
  agri_medium_nfs <- mean(poisson_df_unscale$eulandsystem_farmland_medium)*lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_medium"))]
  agri_high_nfs <- mean(poisson_df_unscale$eulandsystem_farmland_high)*lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("farmland_high"))]
  d_treedensity_nfs <- mean(poisson_df_unscale$treedensity_2018)*
    (sum(lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])/sum(lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium","forest_high"))])-1)/(2050-2018)
  forest_lowmedium_nfs <- mean(poisson_df_unscale$eulandsystem_forest_lowmedium)*lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_lowmedium"))]
  forest_high_nfs <- mean(poisson_df_unscale$eulandsystem_forest_high)*lulc_pls_short$nfs[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]/lulc_pls_short$initial[which(lulc_pls_short$PLS==PLS & lulc_pls_short$variable %in% c("forest_high"))]
  protectedarea_perc_nfs <- mean(poisson_df_unscale$protectedarea_perc)*pa_pls_short$nfs[which(pa_pls_short$PLS==PLS)]/pa_pls_short$initial[which(pa_pls_short$PLS==PLS)]
  
  beta1_nfs <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"]*(d_impervious_nfs - d_impervious_mu)/d_impervious_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nfs - d_tempspring_mu)/d_tempspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nfs - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"]*(d_precspring_nfs - d_precspring_mu)/d_precspring_si +
    mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"]*(d_shannon_nfs - d_shannon_mu)/d_shannon_si +
    mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nfs - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nfs - d_treedensity_mu)/d_treedensity_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nfs - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nfs - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"]*(d_agri_nfs - d_agri_mu)/d_agri_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nfs - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nfs - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nfs - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfs_sample <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_impervious"),"Std. Error"])*(d_impervious_nfs - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nfs - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nfs - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_precspring"),"Std. Error"])*(d_precspring_nfs - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_shannon"),"Std. Error"])*(d_shannon_nfs - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nfs - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nfs - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nfs - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nfs +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:d_agri"),"Std. Error"])*(d_agri_nfs - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nfs - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nfs - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nfs - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  
  mod_coef_signif <- mod_coef
  mod_coef_signif[which(mod_coef_signif[,c("Pr(>|t|)")] > 0.05),c("Estimate","Std. Error")] <- 0
  
  beta1_past_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*mean(poisson_df$d_impervious) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*mean(poisson_df$d_tempsrping) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*mean(poisson_df$d_tempsrpingvar) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*mean(poisson_df$d_precspring) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*mean(poisson_df$d_shannon) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*mean(poisson_df$protectedarea_perc) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*mean(poisson_df$d_treedensity) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*mean(poisson_df$eulandsystem_forest_lowmedium) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*mean(poisson_df$eulandsystem_forest_high) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*mean(poisson_df$d_agri) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_low) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_medium) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_high)
  
  beta1_past_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*mean(poisson_df$d_impervious) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*mean(poisson_df$d_tempsrping) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*mean(poisson_df$d_tempsrpingvar) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*mean(poisson_df$d_precspring) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*mean(poisson_df$d_shannon) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*mean(poisson_df$protectedarea_perc) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*mean(poisson_df$d_treedensity) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_lowmedium) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_high) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*mean(poisson_df$d_agri) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_low) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_medium) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_high)
  
  beta1_SSP1_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*(d_impervious_ssp1 - d_impervious_mu)/d_impervious_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*(d_shannon_ssp1 - d_shannon_mu)/d_shannon_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*(d_treedensity_ssp1 - d_treedensity_mu)/d_treedensity_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_ssp1 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_ssp1 - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*(d_agri_ssp1 - d_agri_mu)/d_agri_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_ssp1 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_ssp1 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_ssp1 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_SSP1_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*(d_impervious_ssp1 - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*(d_shannon_ssp1 - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_ssp1 - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_ssp1 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_ssp1 +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*(d_agri_ssp1 - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_ssp1 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_ssp1 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_ssp1 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  
  beta1_BAU_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*mean(poisson_df$d_impervious) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*mean(poisson_df$d_shannon) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*mean(poisson_df$d_treedensity) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*mean(poisson_df$eulandsystem_forest_lowmedium) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*mean(poisson_df$eulandsystem_forest_high) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*mean(poisson_df$d_agri) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_low) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_medium) +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*mean(poisson_df$eulandsystem_farmland_high)
  
  beta1_BAU_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*mean(poisson_df$d_impervious) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp1 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp1 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp1 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*mean(poisson_df$d_shannon) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp1 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*mean(poisson_df$d_treedensity) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_lowmedium) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*mean(poisson_df$eulandsystem_forest_high) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*mean(poisson_df$d_agri) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_low) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_medium) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*mean(poisson_df$eulandsystem_farmland_high)
  
  
  beta1_SSP3_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*(d_impervious_ssp3 - d_impervious_mu)/d_impervious_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_ssp3 - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_ssp3 - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_ssp3 - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*(d_shannon_ssp3 - d_shannon_mu)/d_shannon_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_ssp3 - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*(d_treedensity_ssp3 - d_treedensity_mu)/d_treedensity_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_ssp3 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_ssp3 - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*(d_agri_ssp3 - d_agri_mu)/d_agri_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_ssp3 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_ssp3 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_ssp3 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_SSP3_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*(d_impervious_ssp3 - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_ssp3 - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_ssp3 - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_ssp3 - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*(d_shannon_ssp3 - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_ssp3 - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_ssp3 - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_ssp3 - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_ssp3 +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*(d_agri_ssp3 - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_ssp3 - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_ssp3 - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_ssp3 - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nac_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*(d_impervious_nac - d_impervious_mu)/d_impervious_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nac - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nac - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_nac - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*(d_shannon_nac - d_shannon_mu)/d_shannon_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nac - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nac - d_treedensity_mu)/d_treedensity_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nac - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nac - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*(d_agri_nac - d_agri_mu)/d_agri_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nac - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nac - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nac - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nac_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*(d_impervious_nac - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nac - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nac - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_nac - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*(d_shannon_nac - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nac - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nac - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nac - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nac +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*(d_agri_nac - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nac - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nac - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nac - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfn_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*(d_impervious_nfn - d_impervious_mu)/d_impervious_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nfn - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nfn - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_nfn - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*(d_shannon_nfn - d_shannon_mu)/d_shannon_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nfn - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nfn - d_treedensity_mu)/d_treedensity_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nfn - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nfn - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*(d_agri_nfn - d_agri_mu)/d_agri_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nfn - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nfn - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nfn - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfn_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*(d_impervious_nfn - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nfn - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nfn - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_nfn - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*(d_shannon_nfn - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nfn - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nfn - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nfn - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nfn +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*(d_agri_nfn - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nfn - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nfn - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nfn - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfs_signif <- mod_coef[which(row.names(mod_coef)=="year"),"Estimate"] +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"]*(d_impervious_nfs - d_impervious_mu)/d_impervious_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"]*(d_tempspring_nfs - d_tempspring_mu)/d_tempspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"]*(d_tempspringvar_nfs - d_tempspringvar_mu)/d_tempspringvar_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"]*(d_precspring_nfs - d_precspring_mu)/d_precspring_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"]*(d_shannon_nfs - d_shannon_mu)/d_shannon_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"]*(protectedarea_perc_nfs - protectedarea_perc_mu)/protectedarea_perc_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"]*(d_treedensity_nfs - d_treedensity_mu)/d_treedensity_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"]*(forest_lowmedium_nfs - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"]*(forest_high_nfs - eulandsystem_forest_high_mu)/eulandsystem_forest_high_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"]*(d_agri_nfs - d_agri_mu)/d_agri_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"]*(agri_low_nfs - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"]*(agri_medium_nfs - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"]*(agri_high_nfs - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  beta1_nfs_sample_signif <- rnorm(nb_rep,mod_coef[which(row.names(mod_coef)=="year"),"Estimate"], sd=mod_coef[which(row.names(mod_coef)=="year"),"Std. Error"]) +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_impervious"),"Std. Error"])*(d_impervious_nfs - d_impervious_mu)/d_impervious_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrping"),"Std. Error"])*(d_tempspring_nfs - d_tempspring_mu)/d_tempspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_tempsrpingvar"),"Std. Error"])*(d_tempspringvar_nfs - d_tempspringvar_mu)/d_tempspringvar_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_precspring"),"Std. Error"])*(d_precspring_nfs - d_precspring_mu)/d_precspring_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_shannon"),"Std. Error"])*(d_shannon_nfs - d_shannon_mu)/d_shannon_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:protectedarea_perc"),"Std. Error"])*(protectedarea_perc_nfs - protectedarea_perc_mu)/protectedarea_perc_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_treedensity"),"Std. Error"])*(d_treedensity_nfs - d_treedensity_mu)/d_treedensity_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_lowmedium"),"Std. Error"])*(forest_lowmedium_nfs - eulandsystem_forest_lowmedium_mu)/eulandsystem_forest_lowmedium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_forest_high"),"Std. Error"])*forest_high_nfs +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:d_agri"),"Std. Error"])*(d_agri_nfs - d_agri_mu)/d_agri_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_low"),"Std. Error"])*(agri_low_nfs - eulandsystem_farmland_low_mu)/eulandsystem_farmland_low_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_medium"),"Std. Error"])*(agri_medium_nfs - eulandsystem_farmland_medium_mu)/eulandsystem_farmland_medium_si +
    rnorm(nb_rep,mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Estimate"], sd=mod_coef_signif[which(row.names(mod_coef_signif)=="year:eulandsystem_farmland_high"),"Std. Error"])*(agri_high_nfs - eulandsystem_farmland_high_mu)/eulandsystem_farmland_high_si
  
  
  return(data.frame(intercept = summary(mod)$p.table[1,1], trend_past=beta1_past,sd_past=sd(beta1_past_sample), trend_BAU=beta1_BAU,sd_BAU=sd(beta1_BAU_sample),trend_SSP1=beta1_SSP1,sd_SSP1=sd(beta1_SSP1_sample),
                    trend_SSP3=beta1_SSP3,sd_SSP3=sd(beta1_SSP3_sample),trend_nac=beta1_nac,sd_nac=sd(beta1_nac_sample),
                    trend_nfn=beta1_nfn,sd_nfn=sd(beta1_nfn_sample),trend_nfs=beta1_nfs,sd_nfs=sd(beta1_nfs_sample),
                    trend_past_signif=beta1_past_signif,sd_past_signif=sd(beta1_past_sample_signif),
                    trend_BAU_signif=beta1_BAU_signif,sd_BAU_signif=sd(beta1_BAU_sample_signif),trend_SSP1_signif=beta1_SSP1_signif,sd_SSP1_signif=sd(beta1_SSP1_sample_signif),
                    trend_SSP3_signif=beta1_SSP3_signif,sd_SSP3_signif=sd(beta1_SSP3_sample_signif),trend_nac_signif=beta1_nac_signif,sd_nac_signif=sd(beta1_nac_sample_signif),
                    trend_nfn_signif=beta1_nfn_signif,sd_nfn_signif=sd(beta1_nfn_sample_signif),trend_nfs_signif=beta1_nfs_signif,sd_nfs_signif=sd(beta1_nfs_sample_signif),PLS=PLS, pressure_removed))
}

# trend by scenario for birds

predict_trend_bird <- function(bird_data,pressure_data,pressure_data_unscale,
                               lulc_pls_short,climate_pls,pa_pls_short,
                               pressure_name = c("d_impervious","d_treedensity","d_agri",
                                                 "d_tempsrping","tempsrping","d_tempsrpingvar","d_precspring","precspring",
                                                 "d_shannon","shannon","drymatter","protectedarea_perc",
                                                 "eulandsystem_farmland_low","eulandsystem_farmland_medium","eulandsystem_farmland_high",
                                                 "eulandsystem_forest_lowmedium","eulandsystem_forest_high","milieu_cat"),
                               min_site_number_per_species = 60,
                               min_occurence_species=200,
                               family="quasipoisson"){
  
  species_press_data_year <- merge(bird_data, pressure_data[which(pressure_data$siteID %in% unique(bird_data$siteID) & pressure_data$year %in% unique(bird_data$year)),], by =c("siteID","year"), all.x=TRUE)
  
  poisson_df <- na.omit(species_press_data_year[,c("siteID","count","year","area_sampled_m2","scheme_code","Long_LAEA","Lat_LAEA",
                                                   pressure_name,"PLS")])
  
  species_press_data_year_unscale <- merge(bird_data, pressure_data_unscale[which(pressure_data_unscale$siteID %in% unique(bird_data$siteID) & pressure_data_unscale$year %in% unique(bird_data$year)),], by =c("siteID","year"), all.x=TRUE)
  
  poisson_df_unscale <- na.omit(species_press_data_year_unscale[,c("siteID","count","year","area_sampled_m2","scheme_code","Long_LAEA","Lat_LAEA",
                                                                   pressure_name,"tempspring_2020","tempspringvar_2020","precspring_2020","agri_2018","shannon_2018","impervious_2018","treedensity_2018","PLS")])
  
  
  poisson_df$year <- scale(poisson_df$year)
  
  if(length(table(poisson_df$area_sampled_m2)) > length(unique(poisson_df$scheme_code))){
    one_scheme_time_area <- 0 
    poisson_df$area_sampled_m2 <- scale(poisson_df$area_sampled_m2)
  }else{
    one_scheme_time_area <- 1
  }
  
  poisson_df$count_scale_all <- poisson_df$count
  
  if(length(pressure_name) > 1){
    formula_gam <- "count_scale_all ~ year + year:d_impervious + year:d_treedensity +
    year:eulandsystem_forest_lowmedium + year:eulandsystem_forest_high +
    year:d_agri + year:eulandsystem_farmland_low + year:eulandsystem_farmland_medium + year:eulandsystem_farmland_high +
    year:d_tempsrping + year:d_tempsrpingvar + year:d_precspring + year:d_shannon + year:protectedarea_perc +
    milieu_cat + tempsrping + precspring + shannon + drymatter"
  }else{
    formula_gam <- paste("count_scale_all ~", paste(pressure_name,sep="", collapse = " + "))
  }
  
  col_names <- c("(Intercept)","year","milieu_catopenland","milieu_catothers","milieu_caturban",
                 "tempsrping","precspring","shannon","drymatter","year:d_impervious","year:d_treedensity",
                 "year:eulandsystem_forest_lowmedium","year:eulandsystem_forest_high","year:d_agri",
                 "year:eulandsystem_farmland_low","year:eulandsystem_farmland_medium","year:eulandsystem_farmland_high",
                 "year:d_tempsrping","year:d_tempsrpingvar","year:d_precspring","year:d_shannon","year:protectedarea_perc")
  
  
  if(nrow(poisson_df) >= min_occurence_species){
    
    ### global poisson model
    
    if(length(unique(poisson_df$scheme_code)) > 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2:scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) == 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) > 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$scheme_code)) == 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    
    if(global_mod$converged){
      
      predict_trend_europe <- predict_trend(mod=global_mod,
                                            pressure_data_unscale,
                                            poisson_df_unscale,
                                            poisson_df,
                                            lulc_pls_short,
                                            climate_pls,
                                            pa_pls_short,
                                            PLS="europe")
      
      for(i in 1:14){
        predict_trend_europe_temp <- predict_trend(mod=global_mod,
                                                   pressure_data_unscale,
                                                   poisson_df_unscale,
                                                   poisson_df,
                                                   lulc_pls_short,
                                                   climate_pls,
                                                   pa_pls_short,
                                                   PLS="europe",
                                                   pressure_remove = i)
        
        predict_trend_europe <- rbind(predict_trend_europe,predict_trend_europe_temp)
      }
      
      unique_poisson_df <- distinct(poisson_df, Long_LAEA, Lat_LAEA,.keep_all = TRUE)
      
      if_fail <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                            trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                            trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                            trend_past_signif=NA,sd_past_signif=NA,
                            trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                            trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                            trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                            PLS=NA,pressure_removed=NA)
      
      predict_trend_pls <- ddply(unique_poisson_df,.(PLS),.fun=purrr::possibly(otherwise=if_fail,
                                                                               .f=function(x,min_site_number_per_species,poisson_df){
                                                                                 
                                                                                 if(nrow(x) >= min_site_number_per_species){
                                                                                   
                                                                                   poisson_df_i <- poisson_df[which(poisson_df$PLS == unique(x$PLS)),]
                                                                                   poisson_df_unscale_i <- poisson_df_unscale[which(poisson_df_unscale$PLS == unique(x$PLS)),]
                                                                                   
                                                                                   if(length(table(poisson_df_i$area_sampled_m2)) > length(unique(poisson_df_i$scheme_code))){
                                                                                     one_scheme_time_area <- 0 
                                                                                     poisson_df_i$area_sampled_m2 <- scale(poisson_df_i$area_sampled_m2)
                                                                                   }else{
                                                                                     one_scheme_time_area <- 1
                                                                                   }
                                                                                   
                                                                                   if(length(unique(poisson_df_i$scheme_code)) > 1 && one_scheme_time_area == 0){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2:scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$scheme_code)) == 1 && one_scheme_time_area == 0){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$scheme_code)) > 1 && one_scheme_time_area == 1){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("scheme_code","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$scheme_code)) == 1 && one_scheme_time_area == 1){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(siteID=~1))
                                                                                   }
                                                                                   
                                                                                   predict_trend_i <- predict_trend(mod=res.poisson_i$gam,
                                                                                                                    pressure_data_unscale,
                                                                                                                    poisson_df_unscale_i,
                                                                                                                    poisson_df_i,
                                                                                                                    lulc_pls_short,
                                                                                                                    climate_pls,
                                                                                                                    pa_pls_short,
                                                                                                                    PLS=unique(x$PLS))
                                                                                   
                                                                                 }else{
                                                                                   predict_trend_i <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                                                                                                 trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                                                                                                 trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                                                                                                 trend_past_signif=NA,sd_past_signif=NA,
                                                                                                                 trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                                                                                                 trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                                                                                                 trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                                                                                                 PLS=unique(x$PLS),pressure_removed=NA)
                                                                                 }
                                                                                 
                                                                                 return(predict_trend_i)
                                                                               }),
                                 min_site_number_per_species=min_site_number_per_species,poisson_df=poisson_df,
                                 .progress="none")
      
      predict_trend_all <- rbind(predict_trend_pls,predict_trend_europe)
      
      
    }else{
      predict_trend_all <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                      trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                      trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                      trend_past_signif=NA,sd_past_signif=NA,
                                      trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                      trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                      trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                      PLS=NA,pressure_removed=NA)
    }
    
  }else{
    predict_trend_all <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                    trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                    trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                    trend_past_signif=NA,sd_past_signif=NA,
                                    trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                    trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                    trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                    PLS=NA,pressure_removed=NA)
  }
  
  return(predict_trend_all)
}




# trend by scenario for butterflies

predict_trend_butterfly <- function(butterfly_data,pressure_data,pressure_data_unscale,
                                    lulc_pls_short,climate_pls,pa_pls_short,
                                    pressure_name = c("d_impervious","d_treedensity","d_agri",
                                                      "d_tempsrping","tempsrping","d_tempsrpingvar","d_precspring","precspring",
                                                      "d_shannon","shannon","drymatter","protectedarea_perc",
                                                      "eulandsystem_farmland_low","eulandsystem_farmland_medium","eulandsystem_farmland_high",
                                                      "eulandsystem_forest_lowmedium","eulandsystem_forest_high","milieu_cat"),
                                    min_site_number_per_species = 60,
                                    min_occurence_species=200,
                                    family="quasipoisson"){
  
  species_press_data_year <- merge(butterfly_data, pressure_data[which(pressure_data$transect_id %in% unique(butterfly_data$transect_id) & pressure_data$year %in% unique(butterfly_data$year)),], by =c("transect_id","year"), all.x=TRUE)
  
  poisson_df <- na.omit(species_press_data_year[,c("transect_id","count_corrected","year","transect_length","bms_id","Long_LAEA","Lat_LAEA",
                                                   pressure_name,"PLS")])
  
  species_press_data_year_unscale <- merge(butterfly_data, pressure_data_unscale[which(pressure_data_unscale$transect_id %in% unique(butterfly_data$transect_id) & pressure_data_unscale$year %in% unique(butterfly_data$year)),], by =c("transect_id","year"), all.x=TRUE)
  
  poisson_df_unscale <- na.omit(species_press_data_year_unscale[,c("transect_id","count_corrected","year","transect_length","bms_id","Long_LAEA","Lat_LAEA",
                                                                   pressure_name,"tempspring_2020","tempspringvar_2020","precspring_2020","agri_2018","shannon_2018","impervious_2018","treedensity_2018","PLS")])
  
  
  poisson_df$year <- scale(poisson_df$year)
  
  if(length(table(poisson_df$transect_length)) > length(unique(poisson_df$bms_id))){
    one_scheme_time_area <- 0 
    poisson_df$transect_length <- scale(poisson_df$transect_length)
  }else{
    one_scheme_time_area <- 1
  }
  
  poisson_df$count_corrected_scale_all <- poisson_df$count_corrected
  
  if(length(pressure_name) > 1){
    formula_gam <- "count_corrected_scale_all ~ year + year:d_impervious + year:d_treedensity +
    year:eulandsystem_forest_lowmedium + year:eulandsystem_forest_high +
    year:d_agri + year:eulandsystem_farmland_low + year:eulandsystem_farmland_medium + year:eulandsystem_farmland_high +
    year:d_tempsrping + year:d_tempsrpingvar + year:d_precspring + year:d_shannon + year:protectedarea_perc +
    milieu_cat + tempsrping + precspring + shannon + drymatter"
  }else{
    formula_gam <- paste("count_corrected_scale_all ~", paste(pressure_name,sep="", collapse = " + "))
  }
  
  col_names <- c("(Intercept)","year","milieu_catopenland","milieu_catothers","milieu_caturban",
                 "tempsrping","precspring","shannon","drymatter","year:d_impervious","year:d_treedensity",
                 "year:eulandsystem_forest_lowmedium","year:eulandsystem_forest_high","year:d_agri",
                 "year:eulandsystem_farmland_low","year:eulandsystem_farmland_medium","year:eulandsystem_farmland_high",
                 "year:d_tempsrping","year:d_tempsrpingvar","year:d_precspring","year:d_shannon","year:protectedarea_perc")
  
  
  if(nrow(poisson_df) >= min_occurence_species){
    
    ### global poisson model
    
    if(length(unique(poisson_df$bms_id)) > 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length:bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) == 1 && one_scheme_time_area == 0){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) > 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    if(length(unique(poisson_df$bms_id)) == 1 && one_scheme_time_area == 1){
      global_mod <- gam(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                        family=family, data=poisson_df)
    }
    
    if(global_mod$converged){
      
      predict_trend_europe <- predict_trend(mod=global_mod,
                                            pressure_data_unscale,
                                            poisson_df_unscale,
                                            poisson_df,
                                            lulc_pls_short,
                                            climate_pls,
                                            pa_pls_short,
                                            PLS="europe")
      
      for(i in 1:14){
        predict_trend_europe_temp <- predict_trend(mod=global_mod,
                                                   pressure_data_unscale,
                                                   poisson_df_unscale,
                                                   poisson_df,
                                                   lulc_pls_short,
                                                   climate_pls,
                                                   pa_pls_short,
                                                   PLS="europe",
                                                   pressure_remove = i)
        
        predict_trend_europe <- rbind(predict_trend_europe,predict_trend_europe_temp)
      }
      
      unique_poisson_df <- distinct(poisson_df, Long_LAEA, Lat_LAEA,.keep_all = TRUE)
      
      if_fail <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                            trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                            trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                            trend_past_signif=NA,sd_past_signif=NA,
                            trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                            trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                            trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                            PLS=NA,pressure_removed=NA)
      
      predict_trend_pls <- ddply(unique_poisson_df,.(PLS),.fun=purrr::possibly(otherwise=if_fail,
                                                                               .f=function(x,min_site_number_per_species,poisson_df){
                                                                                 
                                                                                 if(nrow(x) >= min_site_number_per_species){
                                                                                   
                                                                                   poisson_df_i <- poisson_df[which(poisson_df$PLS == unique(x$PLS)),]
                                                                                   poisson_df_unscale_i <- poisson_df_unscale[which(poisson_df_unscale$PLS == unique(x$PLS)),]
                                                                                   
                                                                                   if(length(table(poisson_df_i$transect_length)) > length(unique(poisson_df_i$bms_id))){
                                                                                     one_scheme_time_area <- 0 
                                                                                     poisson_df_i$transect_length <- scale(poisson_df_i$transect_length)
                                                                                   }else{
                                                                                     one_scheme_time_area <- 1
                                                                                   }
                                                                                   
                                                                                   if(length(unique(poisson_df_i$bms_id)) > 1 && one_scheme_time_area == 0){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length:bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$bms_id)) == 1 && one_scheme_time_area == 0){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("transect_length","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$bms_id)) > 1 && one_scheme_time_area == 1){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("bms_id","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   }
                                                                                   if(length(unique(poisson_df_i$bms_id)) == 1 && one_scheme_time_area == 1){
                                                                                     res.poisson_i <- gamm(as.formula(paste(formula_gam,sep=" + ",paste(c("te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=4)"), collapse = " + "))),
                                                                                                           family=family, data=poisson_df_i,random=list(transect_id=~1))
                                                                                   }
                                                                                   
                                                                                   predict_trend_i <- predict_trend(mod=res.poisson_i$gam,
                                                                                                                    pressure_data_unscale,
                                                                                                                    poisson_df_unscale_i,
                                                                                                                    poisson_df_i,
                                                                                                                    lulc_pls_short,
                                                                                                                    climate_pls,
                                                                                                                    pa_pls_short,
                                                                                                                    PLS=unique(x$PLS))
                                                                                   
                                                                                 }else{
                                                                                   predict_trend_i <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                                                                                                 trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                                                                                                 trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                                                                                                 trend_past_signif=NA,sd_past_signif=NA,
                                                                                                                 trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                                                                                                 trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                                                                                                 trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                                                                                                 PLS=unique(x$PLS),pressure_removed=NA)
                                                                                 }
                                                                                 
                                                                                 return(predict_trend_i)
                                                                               }),
                                 min_site_number_per_species=min_site_number_per_species,poisson_df=poisson_df,
                                 .progress="none")
      
      predict_trend_all <- rbind(predict_trend_pls,predict_trend_europe)

      
    }else{
      predict_trend_all <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                      trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                      trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                      trend_past_signif=NA,sd_past_signif=NA,
                                      trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                      trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                      trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                      PLS=NA,pressure_removed=NA)
    }
    
  }else{
    predict_trend_all <- data.frame(intercept=NA,trend_past=NA,sd_past=NA,trend_BAU=NA,sd_BAU=NA,trend_SSP1=NA,sd_SSP1=NA,
                                    trend_SSP3=NA,sd_SSP3=NA,trend_nac=NA,sd_nac=NA,
                                    trend_nfn=NA,sd_nfn=NA,trend_nfs=NA,sd_nfs=NA,
                                    trend_past_signif=NA,sd_past_signif=NA,
                                    trend_BAU_signif=NA,sd_BAU_signif=NA,trend_SSP1_signif=NA,sd_SSP1_signif=NA,
                                    trend_SSP3_signif=NA,sd_SSP3_signif=NA,trend_nac_signif=NA,sd_nac_signif=NA,
                                    trend_nfn_signif=NA,sd_nfn_signif=NA,trend_nfs_signif=NA,sd_nfs_signif=NA,
                                    PLS=NA,pressure_removed=NA)
  }
  
  return(predict_trend_all)
}



overall_mean_sd_trend <- function(data){
  n <- length(na.omit(data$trend_BAU))
  mu_past <- exp(mean(data$trend_past,na.rm=TRUE))
  var_past <- (sum(data$sd_past^2 + data$trend_past^2, na.rm = TRUE))/n - mean(data$trend_past,na.rm=TRUE)^2
  sd_past <- sqrt(mu_past^2*var_past)
  se_past <- mu_past/sqrt(n)*sd(data$trend_past)
  mu_bau <- exp(mean(data$trend_BAU,na.rm=TRUE))
  var_bau <- (sum(data$sd_BAU^2 + data$trend_BAU^2, na.rm = TRUE))/n - mean(data$trend_BAU,na.rm=TRUE)^2
  sd_bau <- sqrt(mu_bau^2*var_bau)
  se_bau <- mu_bau/sqrt(n)*sd(data$trend_BAU)
  mu_ssp1 <- exp(mean(data$trend_SSP1,na.rm=TRUE))
  var_ssp1 <- (sum(data$sd_SSP1^2 + data$trend_SSP1^2, na.rm = TRUE))/n - mean(data$trend_SSP1,na.rm=TRUE)^2
  sd_ssp1 <- sqrt(mu_ssp1^2*var_ssp1)
  se_ssp1 <- mu_ssp1/sqrt(n)*sd(data$trend_SSP1)
  mu_ssp3 <- exp(mean(data$trend_SSP3,na.rm=TRUE))
  var_ssp3 <- (sum(data$sd_SSP3^2 + data$trend_SSP3^2, na.rm = TRUE))/n - mean(data$trend_SSP3,na.rm=TRUE)^2
  sd_ssp3 <- sqrt(mu_ssp3^2*var_ssp3)
  se_ssp3 <- mu_ssp3/sqrt(n)*sd(data$trend_SSP3)
  mu_nac <- exp(mean(data$trend_nac,na.rm=TRUE))
  var_nac <- (sum(data$sd_nac^2 + data$trend_nac^2, na.rm = TRUE))/n - mean(data$trend_nac,na.rm=TRUE)^2
  sd_nac <- sqrt(mu_nac^2*var_nac)
  se_nac <- mu_nac/sqrt(n)*sd(data$trend_nac)
  mu_nfn <- exp(mean(data$trend_nfn,na.rm=TRUE))
  var_nfn <- (sum(data$sd_nfn^2 + data$trend_nfn^2, na.rm = TRUE))/n - mean(data$trend_nfn,na.rm=TRUE)^2
  sd_nfn <- sqrt(mu_nfn^2*var_nfn)
  se_nfn <- mu_nfn/sqrt(n)*sd(data$trend_nfn)
  mu_nfs <- exp(mean(data$trend_nfs,na.rm=TRUE))
  var_nfs <- (sum(data$sd_nfs^2 + data$trend_nfs^2, na.rm = TRUE))/n - mean(data$trend_nfs,na.rm=TRUE)^2
  sd_nfs <- sqrt(mu_nfs^2*var_nfs)
  se_nfs <- mu_nfs/sqrt(n)*sd(data$trend_nfs)
  
  mu_past_signif <- exp(mean(data$trend_past_signif,na.rm=TRUE))
  var_past_signif <- (sum(data$sd_past_signif^2 + data$trend_past_signif^2, na.rm = TRUE))/n - mean(data$trend_past_signif,na.rm=TRUE)^2
  sd_past_signif <- sqrt(mu_past_signif^2*var_past_signif)
  se_past_signif <- mu_past_signif/sqrt(n)*sd(data$trend_past_signif)
  mu_bau_signif <- exp(mean(data$trend_BAU_signif,na.rm=TRUE))
  var_bau_signif <- (sum(data$sd_BAU_signif^2 + data$trend_BAU_signif^2, na.rm = TRUE))/n - mean(data$trend_BAU_signif,na.rm=TRUE)^2
  sd_bau_signif <- sqrt(mu_bau_signif^2*var_bau_signif)
  se_bau_signif <- mu_bau_signif/sqrt(n)*sd(data$trend_BAU_signif)
  mu_ssp1_signif <- exp(mean(data$trend_SSP1_signif,na.rm=TRUE))
  var_ssp1_signif <- (sum(data$sd_SSP1_signif^2 + data$trend_SSP1_signif^2, na.rm = TRUE))/n - mean(data$trend_SSP1_signif,na.rm=TRUE)^2
  sd_ssp1_signif <- sqrt(mu_ssp1_signif^2*var_ssp1_signif)
  se_ssp1_signif <- mu_ssp1_signif/sqrt(n)*sd(data$trend_SSP1_signif)
  mu_ssp3_signif <- exp(mean(data$trend_SSP3_signif,na.rm=TRUE))
  var_ssp3_signif <- (sum(data$sd_SSP3_signif^2 + data$trend_SSP3_signif^2, na.rm = TRUE))/n - mean(data$trend_SSP3_signif,na.rm=TRUE)^2
  sd_ssp3_signif <- sqrt(mu_ssp3_signif^2*var_ssp3_signif)
  se_ssp3_signif <- mu_ssp3_signif/sqrt(n)*sd(data$trend_SSP3_signif)
  mu_nac_signif <- exp(mean(data$trend_nac_signif,na.rm=TRUE))
  var_nac_signif <- (sum(data$sd_nac_signif^2 + data$trend_nac_signif^2, na.rm = TRUE))/n - mean(data$trend_nac_signif,na.rm=TRUE)^2
  sd_nac_signif <- sqrt(mu_nac_signif^2*var_nac_signif)
  se_nac_signif <- mu_nac_signif/sqrt(n)*sd(data$trend_nac_signif)
  mu_nfn_signif <- exp(mean(data$trend_nfn_signif,na.rm=TRUE))
  var_nfn_signif <- (sum(data$sd_nfn_signif^2 + data$trend_nfn_signif^2, na.rm = TRUE))/n - mean(data$trend_nfn_signif,na.rm=TRUE)^2
  sd_nfn_signif <- sqrt(mu_nfn_signif^2*var_nfn_signif)
  se_nfn_signif <- mu_nfn_signif/sqrt(n)*sd(data$trend_nfn_signif)
  mu_nfs_signif <- exp(mean(data$trend_nfs_signif,na.rm=TRUE))
  var_nfs_signif <- (sum(data$sd_nfs_signif^2 + data$trend_nfs_signif^2, na.rm = TRUE))/n - mean(data$trend_nfs_signif,na.rm=TRUE)^2
  sd_nfs_signif <- sqrt(mu_nfs_signif^2*var_nfs_signif)
  se_nfs_signif <- mu_nfs_signif/sqrt(n)*sd(data$trend_nfs_signif)
  
  return(data.frame(mu_past,sd_past,se_past,mu_bau,sd_bau,se_bau,mu_ssp1,sd_ssp1,se_ssp1,mu_ssp3,sd_ssp3,se_ssp3,
                    mu_nac,sd_nac,se_nac,mu_nfn,sd_nfn,se_nfn,mu_nfs,sd_nfs,se_nfs,
                    mu_past_signif,sd_past_signif,se_past_signif,mu_bau_signif,sd_bau_signif,se_bau_signif,mu_ssp1_signif,sd_ssp1_signif,se_ssp1_signif,mu_ssp3_signif,sd_ssp3_signif,se_ssp3_signif,
                    mu_nac_signif,sd_nac_signif,se_nac_signif,mu_nfn_signif,sd_nfn_signif,se_nfn_signif,mu_nfs_signif,sd_nfs_signif,se_nfs_signif,n))
}

