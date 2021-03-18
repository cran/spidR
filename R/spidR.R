#####spidR - Spider Biodiversity Tools
#####Version 0.1.0 (2021-03-17)
#####By Pedro Cardoso
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: None yet.
#####Changed from v0.0.0:
#####Everything new

#####required packages
library("graphics")
library("httr")
library("jsonlite")
library("rgbif")
library("rworldmap")
library("stats")
library("utils")
#' @import graphics
#' @import httr
#' @import jsonlite
#' @import rgbif
#' @import rworldmap
#' @import stats
#' @import utils

#auxiliary functions

wsc.api <- function(sp, api){
  results = c()

  #convert sp to lsid if a species name is given
  if(gregexpr(" ", sp)[[1]][1] > 0){
    wsc <- NULL
    wsc = wsc.data()
    sp = wsc[wsc$name == sp, ]$lsid
    sp = sub("urn:lsid:nmbe.ch:spidersp:", "", sp)
  }

  id = httr::GET(paste("https://wsc.nmbe.ch/api/lsid/urn:lsid:nmbe.ch:spidersp:", sp, "?apiKey=", api, sep = ""))
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Downloads WSC data.
#' @description Downloads the most recent data from the World Spider Catalogue.
#' @details The World Spider Catalog (2021) lists all currently valid species of spiders, from Clerck to date. Updated daily.
#' @return A matrix with all current species names and distribution. This should be used for other functions using wsc data.
#' @references World Spider Catalog (2021). World Spider Catalog. Version 22.0. Natural History Museum Bern, online at http://wsc.nmbe.ch. doi: 10.24436/2.
#' @examples \dontrun{
#' wsc.data()
#' }
#' @export
wsc.data <- function(){

  if(!exists("wsc") || ((Sys.time() - attributes(wsc)$lastUpdate) > 1440)){

    wsc <- NULL

    #fetch data from wsc (try both current and previous day)
    today = gsub("-", "", as.character(Sys.Date()))
    yesterday = gsub("-", "", as.character(Sys.Date()-1))
    wsc = tryCatch(read.csv2(paste("https://wsc.nmbe.ch/resources/species_export_", today, ".csv", sep = ""), sep = ","),
                   warning = function(x) x = read.csv2(paste("https://wsc.nmbe.ch/resources/species_export_", yesterday, ".csv", sep = ""), sep = ","))
    #clean data
    wsc[,1] = do.call(paste, wsc[, 4:5])
    colnames(wsc)[1:2] = c("name", "lsid")
    attr(wsc, 'lastUpdate') = Sys.time()

    #set wsc as global variable
    pos <- 1
    envir = as.environment(pos)
    assign("wsc", wsc, envir = envir)
  }
  return(wsc)
}

#' Check spider species names.
#' @description Check species names against the World Spider Catalogue.
#' @param spp A vector with species names.
#' @details This function will check if all species names in spp are updated according to the World Spider Catalogue (2021). If not, it returns a matrix with mismatches and closer names using fuzzy matching (Levenshtein edit distance).
#' @return If any mismatches, a matrix with species not found in WSC.
#' @references World Spider Catalog (2021). World Spider Catalog. Version 22.0. Natural History Museum Bern, online at http://wsc.nmbe.ch. doi: 10.24436/2.
#' @examples \dontrun{
#' wsc.names(spp = c("Iberesia machadonia", "Nemesia bacelari", "Amphiledorus ungoliantae"))
#' }
#' @export
wsc.names <- function(spp){

  wsc <- NULL
  wsc = wsc.data()

  mismatches = spp[!(spp %in% wsc[,1])]
  if(length(mismatches) == 0) {
    return("All species OK!")
  } else {
    mismatches = cbind(mismatches, rep(NA, length(mismatches)))
    colnames(mismatches) = c("Not Found", "Possible match")
    for(i in 1:nrow(mismatches)){
      d = adist(wsc$name, mismatches[i, 1])
      mismatches[i, 2] = wsc$name[which(d == min(d))[1]]
    }
    return(mismatches)
  }
}

#' Downloads WST data.
#' @description Downloads the most recent data from the World Spider Trait database.
#' @param spp A vector with species names.
#' @param trait A vector with required traits as abbreviations. Valid values can be found at: https://spidertraits.sci.muni.cz/traits
#' @details The World Spider Trait database (Pekar et al. 2021) has been designed to contain trait data in a broad sense, from morphological traits to ecological characteristics, ecophysiology, behavioural habits, and more (Lowe et al. 2020). This function will download everything available for the species given.
#' @return A matrix with trait data.
#' @references Lowe, E., Wolff, J.O., Aceves-Aparicio, A., Birkhofer, K., Branco, V.V., Cardoso, P., Chichorro, F., Fukushima, C.S., Goncalves-Souza, T., Haddad, C.R., Isaia, M., Krehenwinkel, H., Audisio, T.L., Macias-Hernandez, N., Malumbres-Olarte, J., Mammola, S., McLean, D.J., Michalko, R., Nentwig, W., Pekar, S., Petillon, J., Privet, K., Scott, C., Uhl, G., Urbano-Tenorio, F., Wong, B.H. & Herbestein, M.E. (2020). Towards establishment of a centralized spider traits database. Journal of Arachnology, 48: 103-109. https://doi.org/10.1636/0161-8202-48.2.103
#' @references Pekar, S., Cernecka, L., Wolff, J., Mammola, S., Cardoso, P., Lowe, E., Fukushima, C.S., Birkhofer, K. & Herberstein, M.E. (2021). The spider trait database. Masaryk University, Brno, URL: https://spidertraits.sci.muni.cz
#' @examples \dontrun{
#' wst.data("Atypus affinis")
#' wst.data("Zodarion costapratae", trait = "bole")
#' wst.data(c("Iberesia machadoi", "Zodarion costapratae"))
#' wst.data(c("Iberesia machadoi", "Zodarion costapratae"), trait = c("balo", "bole"))
#' }
#' @export
wst.data <- function(spp, trait = NULL){

  wsc <- NULL
  wsc = wsc.data()
  if(wsc.names(spp)[1] != "All species OK!")
    return("Species not found in WSC")

  results = c()
  for(s in spp){
    s = sub(" ", "+", s)
    id = httr::GET(paste("https://spidertraits.sci.muni.cz/backend/taxonomy?offset=0&limit=10&searchField=fullName&searchValue=", s, "&count=true", sep = ""))
    id = content(id)$items[[1]]$id
    traits = httr::GET(paste("https://spidertraits.sci.muni.cz/backend/data/family/*/genus/*/species/", id, "/original-name/*/trait-category/*/trait/*/method/*/location/*/country/*/dataset/*/authors/*/reference/*/row-link/*?offset=0&limit=10000", sep = ""))
    traits = as.data.frame(jsonlite::fromJSON(content(traits, "text"), flatten = TRUE)$items)
    if(nrow(traits) > 0 && !is.null(trait))
      traits = traits[traits$trait.abbrev %in% trait, ]
    if(nrow(traits) > 0)
      results = rbind(results, traits)
  }
  if(is.null(results))
    return("No results found!")
  colnames(results) = gsub("items.", "", colnames(results))
  return(results)
}

#' Downloads coordinate data.
#' @description Downloads coordinate data from GBIF and the World Spider Trait Database.
#' @param sp A single species name.
#' @details Outputs non-duplicate records with geographical (long, lat) coordinates. As always when using data from multiple sources the user should be careful and check if records "make sense" before using them.
#' @return A data.frame with longitude and latitude.
#' @references Pekar, S., Cernecka, L., Wolff, J., Mammola, S., Cardoso, P., Lowe, E., Fukushima, C.S., Birkhofer, K. & Herberstein, M.E. (2021). The spider trait database. Masaryk University, Brno, URL: https://spidertraits.sci.muni.cz
#' @examples \dontrun{
#' records("Pardosa hyperborea")
#' }
#' @export
records <- function(sp){

  wsc <- NULL
  wsc = wsc.data()
  if(wsc.names(sp)[1] != "All species OK!")
    return("Species not found in WSC")

  results = c()

  #get GBIF data
  gdata = occ_data(scientificName = sp)
  gdata = cbind(gdata$data$decimalLongitude, gdata$data$decimalLatitude)

  #get WST data
  sp = sub(" ", "+", sp)
  id = httr::GET(paste("https://spidertraits.sci.muni.cz/backend/taxonomy?offset=0&limit=10&searchField=fullName&searchValue=", sp, "&count=true", sep = ""))
  id = content(id)$items[[1]]$id
  traits = httr::GET(paste("https://spidertraits.sci.muni.cz/backend/data/family/*/genus/*/species/", id, "/original-name/*/trait-category/*/trait/*/method/*/location/*/country/*/dataset/*/authors/*/reference/*/row-link/*?offset=0&limit=10000", sep = ""))
  traits = as.data.frame(jsonlite::fromJSON(content(traits, "text"), flatten = TRUE)$items)
  tdata = cbind(traits$location.coords.lon, traits$location.coords.lat)

  #merge GBIF and WST data
  results = rbind(gdata, tdata)
  colnames(results) = c("long", "lat")

  #remove NA and duplicate values
  results = results[complete.cases(results), ]
  results = unique(results)

  return(results)

}

#' Map species range.
#' @description Maps species range according to the World Spider Catalogue and records according to GBIF and the World Spider Trait database.
#' @param sp A single species name.
#' @param countries Maps countries according to WSC.
#' @param records Maps records according to GBIF and WST.
#' @param window  Indicates if the map should open in a new window.
#' @details Countries based on the interpretation of the textual descriptions available at the World Spider Catalogue (2021). These might be only approximations to country level and should be taken with caution.
#' @return A world map with countries and records highlighted.
#' @references Pekar, S., Cernecka, L., Wolff, J., Mammola, S., Cardoso, P., Lowe, E., Fukushima, C.S., Birkhofer, K. & Herberstein, M.E. (2021). The spider trait database. Masaryk University, Brno, URL: https://spidertraits.sci.muni.cz
#' @references World Spider Catalog (2021). World Spider Catalog. Version 22.0. Natural History Museum Bern, online at http://wsc.nmbe.ch. doi: 10.24436/2.
#' @examples \dontrun{
#' map("Iberesia machadoi", window = TRUE)
#' map("Pardosa hyperborea", countries = FALSE)
#' }
#' @export
map <- function(sp, countries = TRUE, records = TRUE, window = FALSE){

  wsc <- NULL
  wsc = wsc.data()
  wscmap <- NULL
  data(wscmap, envir = environment())

  if(wsc.names(sp)[1] != "All species OK!")
    return("Species not found in WSC")

  if(countries){
    #get distribution
    distribution = wsc[wsc$name == sp, 10]
    distribution = tolower(strsplit(distribution, c("\\, | and | to |from |\\?|\\/|probably|possibly|introduced"))[[1]]) #split text in chunks
    distribution = sub("\\(.*", "", distribution) #remove everything after parentheses
    distribution = gsub("^\\s+|\\s+$", "", distribution) #remove trailing and white spaces

    #convert distribution to ISO codes
    iso = wscmap[wscmap[,1] %in% distribution, -1]
    iso = colSums(iso)
    iso = names(iso[iso > 0])

    #merge iso codes with map
    iso = unique(iso)
  } else {
    iso = c("STP")
  }
  iso = data.frame(code = iso, exists = rep(1, length(iso)))
  countryRegions <- joinCountryData2Map(iso, joinCode = "ISO3", nameJoinColumn = "code")

  #plot map
  if(window)
    mapDevice()
  mapCountryData(countryRegions, nameColumnToPlot = "exists", mapTitle = sp, catMethod = "categorical", addLegend = FALSE)
  if(records)
    points(records(sp))

}

#' Distribution of spiders at the World Spider Catalogue.
#'
#'A dataset that links species distribution descriptions with the map using the ISO3 code
#'
#' @docType data
#' @keywords datasets
#' @name wscmap
#' @usage data(wscmap)
#' @format A matrix with regions and corresponding ISO3 codes.
NULL
