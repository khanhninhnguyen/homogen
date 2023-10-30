
# draft -------------------------------------------------------------------

# support for segmentation
extract <- function(fichier, path_file, name.folder) {
  # In : fichier de donnees d'IWV journaliere
  # forme : yymmdd.HHMMSS  delta_IWV   IWV_ERAI    sigma_GPS      sigma_ERAI
  # Out : data.frame(name_station,date_num,deltaIWV)
  # renvoie NA si probleme
  
  # fichier='scatterIWV3_kit3.txt'
  
  donnees <- read.table(file = fichier, sep = "", skip = 0, header = F, na.strings = NaN, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  # donnees=read.table(file=fichier,sep="",skip=0,header=F,na.strings=NaN,colClasses=c("character",NA,NA,NA,NA))
  # NaN transformees en NA
  # Premiere colonne conservee en chaine de caractere (pour conserver les zeros des dates)
  
  colnames(donnees) <- c("date", "delta_IWV", "IWV_ERAI", "sigma_GPS", "sigma_ERAI")
  # head(donnees)
  # deltaIWV=ERAI - GPS// now it should be GPS-ERAI
  # On recupere le nom de la station
  # Le nom du fichier est sous la forme "scatterIWV3_alic.txt", les 4 dernières lettres composent le nom de la station
  name_station <- unlist(strsplit(fichier, "\\.")) # Le '.' separe le nom du fichier en deux chaines de caracteres
  if (name_station[length(name_station)] != "txt") {
    print("Le fichier n'est pas un .txt") # verification que c'est bien un fichier .txt
    return(NA)
  }
  name_station <- unlist(strsplit(fichier, "\\/")) # séparation du nom du fichier selon le '_'
  name_station <- substr(name_station[length(name_station)], 13, 16) # extract the name of station from scatter file
  
  # Les dates sont de la forme YY/MM/DD
  # Obtention d'un format YYYY/MM/DD pour la conversion en date numerique
  ind <- which(substr(donnees$date, 1, 1) == "9") # liste des indices des dates des annees 90
  if (length(ind) == 0) {
    donnees$date <- paste0("20", donnees$date)
  } else {
    donnees$date[ind] <- paste0("19", donnees$date[ind]) # on ajoute 19 pour les annees 90
    donnees$date[-ind] <- paste0("20", donnees$date[-ind]) # on ajoute 20 pour les annees 2000
  }
  
  date_num <- strptime(donnees$date, "%Y%m%d.%H", tz = "GMT") # conversion en date numerique
  class(date_num)
  unclass(head(date_num))
  
  # annee=substr(donnees$date,1,4) #les 4 premiers caractères de la colonne date
  # mois=substr(donnees$date,5,6) #les deux caracteres suivants
  # jour=substr(donnees$date,7,8)
  
  deltaIWV <- donnees$delta_IWV # deltaIWV = IWV_ERAI - IWV_GPS
  Y <- data.frame(name_station = name_station, date = date_num, signal = deltaIWV, ERAI = donnees$IWV_ERAI, GPS = donnees$IWV_ERAI - deltaIWV)
  Y$month <- as.factor(format(Y$date, format = "%m")) # extraction du mois de la date numerique et conversion en facteur
  Y$year <- as.factor(format(Y$date, format = "%Y")) # extraction de l'annee
  # Y$position=c(1:length(Y$month))
  # Y$annee=as.factor(Y$annee);Y$mois=as.factor(Y$mois);Y$jour=as.factor(Y$jour)  #Pour que l'estimateur robuste fonctionne
  
  file_name <- paste0(path_file, name.folder, "/", name_station, ".RData")
  save(Y, file = file_name)
  return(Y)
}
extract_ngl <- function(fichier, path_file, name.folder) {
  # NOTE THAT IT HAS BEEN MADE BY THE CHANGE IN ORDER OF COLUMN IN GPS AND ERA IN THE NGL RAW DATA
  # INWHICH ALL INFORMATION BETWEEN GPS AND ERA IS SWITCHED
  # In : fichier de donnees d'IWV journaliere
  # forme : yymmdd.HHMMSS  delta_IWV(GPS-ERA)   IWV_GPS    sigma_ERA5      sigma_GPS
  # Out : data.frame(name_station,date_num,deltaIWV)
  # renvoie NA si probleme
  
  # fichier='scatterIWV3_kit3.txt'
  
  donnees <- read.table(
    file = fichier,
    sep = "",
    skip = 1,
    header = FALSE,
    na.strings = NaN,
    colClasses = c("character", "numeric", "numeric", "numeric", "numeric")
  )
  # donnees=read.table(file=fichier,sep="",skip=0,header=F,na.strings=NaN,colClasses=c("character",NA,NA,NA,NA))
  # NaN transformees en NA
  # Premiere colonne conservee en chaine de caractere (pour conserver les zeros des dates)
  
  colnames(donnees) <- c("date", "delta_IWV", "IWV_GPS", "sigma_ERA5", "sigma_GPS")
  # head(donnees)
  # deltaIWV=ERAI - GPS// now it should be GPS-ERAI
  # On recupere le nom de la station
  # Le nom du fichier est sous la forme "scatterIWV3_alic.txt", les 4 dernières lettres composent le nom de la station
  name_station <- unlist(strsplit(fichier, "\\.")) # Le '.' separe le nom du fichier en deux chaines de caracteres
  if (name_station[length(name_station)] != "txt") {
    print("Le fichier n'est pas un .txt") # verification que c'est bien un fichier .txt
    return(NA)
  }
  name_station <- unlist(strsplit(fichier, "\\/")) # séparation du nom du fichier selon le '_'
  name_station <- substr(name_station[length(name_station)], 13, 16) # extract the name of station from scatter file
  
  # Les dates sont de la forme YY/MM/DD
  # Obtention d'un format YYYY/MM/DD pour la conversion en date numerique
  ind <- which(substr(donnees$date, 1, 1) == "9") # liste des indices des dates des annees 90
  if (length(ind) == 0) {
    donnees$date <- paste0("20", donnees$date)
  } else {
    donnees$date[ind] <- paste0("19", donnees$date[ind]) # on ajoute 19 pour les annees 90
    donnees$date[-ind] <- paste0("20", donnees$date[-ind]) # on ajoute 20 pour les annees 2000
  }
  
  date_num <- strptime(donnees$date, "%Y%m%d.%H", tz = "GMT") # conversion en date numerique
  class(date_num)
  unclass(head(date_num, 2))
  
  deltaIWV <- donnees$delta_IWV # deltaIWV = GPS-ERA5
  Y <- data.frame(
    name_station = name_station,
    date = date_num,
    signal = deltaIWV,
    ERAI = donnees$IWV_GPS - deltaIWV,
    GPS = donnees$IWV_GPS
  )
  Y$month <- as.factor(format(Y$date, format = "%m"))
  Y$year <- as.factor(format(Y$date, format = "%Y"))
  
  # file_name = paste0(path_file,name.folder,"/",name_station,".RData")
  # save(Y, file = file_name)
  return(Y)
}

a = read.table(file = paste0("/home/knguyen/Documents/PhD/Results/", "NGL6048.BM_BJ.date_mean.txt"), header = TRUE)
rpt_data$t_break = as.Date(rpt_data$t_break)
a <- infor_all %>% 
  full_join(rpt_data, 
            by = join_by(main == name_main, 
                         brp == t_break,
                         nearby == name_nearby))
