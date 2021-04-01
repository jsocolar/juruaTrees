data_for_AGBmc2 <- readRDS("/Users/jacobsocolar/Dropbox/Work/juruaTrees/data_for_AGBmc2.RDS")
attach(data_for_AGBmc2)

wWD_3sp<-read.csv(file="/Users/jacobsocolar/Downloads/wWD_3sp.csv", 
                  header=T, sep=";", na=c('',' ','#DIV/0!', '#REF!', '#N/A'), 
                  row.names = NULL)
wWD_3gen<-read.csv(file="/Users/jacobsocolar/Downloads/wWD_3gen.csv", 
                   header=T, sep=";", na=c('',' ','#DIV/0!', '#REF!', '#N/A'), 
                   row.names = NULL)
wWD_3fam<-read.csv(file="/Users/jacobsocolar/Downloads/wWD_3fam.csv", 
                   header=T, sep=";", na=c('',' ','#DIV/0!', '#REF!', '#N/A'), 
                   row.names = NULL)

# Create datasets
stem_data_v <- stem_data[stem_data$Habitat == "vz", ]
stem_data_t <- stem_data[stem_data$Habitat == "tf", ]

#-----------------------------------------------#


sp_v <- wWD_3sp[wWD_3sp$Habitat == "vz", ]#92 obs
sp_v<- merge(sp_v, stem_data_v, by.x="sp", by.y = "sp", all.x = T) 
sp_v<-sp_v[,c(1,4:173)]#483 obs

sp_t <- wWD_3sp[wWD_3sp$Habitat == "tf", ]#32 obs
sp_t<- merge(sp_t, stem_data_t, by.x="sp", by.y = "sp", all.x = T) 
sp_t<-sp_t[,c(1,4:173)]#64 obs

gen_v <- wWD_3gen[wWD_3gen$Habitat == "vz", ]#75 obs
gen_v<- merge(gen_v, stem_data_v, by.x="gen", by.y = "gen", all.x = T) 
gen_v<-gen_v[,c(1,4:173)]#483 obs

gen_t <- wWD_3gen[wWD_3gen$Habitat == "tf", ]#27 obs
gen_t<- merge(gen_t, stem_data_t, by.x="gen", by.y = "gen", all.x = T) 
gen_t<-gen_t[,c(1,4:173)]#100 obs

fam_v <- wWD_3fam[wWD_3fam$Habitat == "vz", ]#32 obs
fam_v<- merge(fam_v, stem_data_v, by.x="fam", by.y = "fam", all.x = T) 
fam_v<-fam_v[,c(1,4:173)]#621 obs

fam_t <- wWD_3fam[wWD_3fam$Habitat == "tf", ]#17 obs
fam_t<- merge(fam_t, stem_data_t, by.x="fam", by.y = "fam", all.x = T) 
fam_t<-fam_t[,c(1,4:173)]#124 obs


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


sp_v<-completeFun(sp_v, c("wWD", "wWD_s"))#480 obs
sp_t<-completeFun(sp_t, c("wWD", "wWD_s"))#53 obs
gen_v<-completeFun(gen_v, c("wWD", "wWD_g"))#583 obs
gen_t<-completeFun(gen_t, c("wWD", "wWD_g"))#94 obs
fam_v<-completeFun(fam_v, c("wWD", "wWD_f"))#621 obs
fam_t<-completeFun(fam_t, c("wWD", "wWD_f"))#122 obs

HEsp_v<-completeFun(sp_v, c("H_E"))#430 obs
HEsp_t<-completeFun(sp_t, c("H_E"))#25 obs
HEgen_v<-completeFun(gen_v, c("H_E"))#517 obs
HEgen_t<-completeFun(gen_t, c("H_E"))#48 obs
HEfam_v<-completeFun(fam_v, c("H_E"))#551 obs
HEfam_t<-completeFun(fam_t, c("H_E"))#63 obs

HEstem_data_v <- completeFun(stem_data_v, c("H_E"))#568 obs
HEstem_data_t <- completeFun(stem_data_t, c("H_E"))#70 obs
HEworld_v <- completeFun(stem_data_v, c("H_E", "taxonomic_WD_matrix"))#568 obs
HEworld_t <- completeFun(stem_data_t, c("H_E", "taxonomic_WD_matrix"))#70 obs



dataset_list <- list(stem_data_v=stem_data_v, stem_data_t=stem_data_t, 
                     sp_v=sp_v, sp_t=sp_t, gen_v=gen_v, gen_t=gen_t, fam_v=fam_v, fam_t=fam_t, 
                     HEsp_v=HEsp_v, HEsp_t=HEsp_t, HEgen_v=HEgen_v, HEgen_t=HEgen_t, HEfam_v=HEfam_v, HEfam_t=HEfam_t,
                     HEstem_data_v=HEstem_data_v, HEstem_data_t=HEstem_data_t, HEworld_v=HEworld_v, HEworld_t=HEworld_t
                     )
for(i in 1:length(dataset_list)){
  dataset_list[[i]]$taxonomic_WD_matrix <- as.matrix(dataset_list[[i]]$taxonomic_WD_matrix)
}


saveRDS(dataset_list, "/Users/JacobSocolar/Desktop/dataset_list.RDS")


