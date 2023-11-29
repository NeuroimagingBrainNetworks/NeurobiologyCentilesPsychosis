
################################################################################
# SCRIPT TO CALCULATE THE AVERAGE OF THE DIFFERENT BRAIN VOLUME REGIONS BETWEEN 
# BOTH HEMISPHERES, LEFT AND RIGHT (LH & RH), FROM A .CSV FILE.
################################################################################

#  Copyright (C) 2023 University of Seville
# 
#  Written by Natalia García San Martín (ngarcia1@us.es)
# 
#  This file is part of Neurobiology Centiles Psychosis toolkit.
# 
#  Neurobiological Centiles Mapping toolkit is free software: 
#  you can redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
# 
#  Neurobiological Centiles Mapping toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Neurobiological Centiles Mapping toolkit. If not, see 
#  <https://www.gnu.org/licenses/>.

rm(list=ls()) # Previous data cleaning in memory

location <-'D:/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code'
location <-'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code'
# Change to the file path of the file to be used
setwd(paste0(location, "/Data/Centiles/PE"))
# setwd(paste0(location, "/Data/Centiles/FEP"))
NEWDATA <- read.csv("GM_volumes_ALSPAC.csv") 
# NEWDATA <- read.csv("GM_volumes_PAFIP.csv") 


# Set the list of names of the regions of interest (columns)
names <- c('bankssts',
             'caudalanteriorcingulate',
             'caudalmiddlefrontal',
             'cuneus',
             'entorhinal',
             'frontalpole',
             'fusiform',
             'inferiorparietal',
             'inferiortemporal',
             'insula',
             'isthmuscingulate',
             'lateraloccipital',
             'lateralorbitofrontal',
             'lingual',
             'medialorbitofrontal',
             'middletemporal',
             'paracentral',
             'parahippocampal',
             'parsopercularis',
             'parsorbitalis',
             'parstriangularis',
             'pericalcarine',
             'postcentral',
             'posteriorcingulate',
             'precentral',
             'precuneus',
             'rostralanteriorcingulate',
             'rostralmiddlefrontal',
             'superiorfrontal',
             'superiorparietal',
             'superiortemporal',
             'supramarginal',
             'temporalpole',
             'transversetemporal')

d<-c()

for(n in names){
  both_hemispheres <- data.frame(NEWDATA[paste0('lh_',n,'_volume')],NEWDATA[paste0('rh_',n,'_volume')]) # Merge the hemispheres of each region.
  media <- rowMeans(both_hemispheres) # Take the average of each hemisphere
  d <- cbind(d,media) # Merge the averages of all regions
  
}
colnames(d) <- c(names) # Assign names to the columns


METADATA <- read.csv("metadata_rafa.csv") # Open the file containing the age, sex and study information

for (i in 1:length(METADATA$age_at_scan)){
  METADATA$age_months=METADATA$age_at_scan*12 # Convert the age to months
  METADATA$age_days=METADATA$age_at_scan*365.245 + 280 # Convert the age to days
  # METADATA$age_days=METADATA$age_months*30
  if  (METADATA$Sex..1.male.[i]==1){ # Convert the sex (1 or 2) to "Male" or "Female" 
    METADATA$sex[i]="Male"
  } else if(METADATA$Sex..1.male.[i]==2){
    METADATA$sex[i]="Female"
  }
  METADATA
  if  (METADATA$FJPL163_clinical_num[i]==1){ # Convert the diagnosis
    METADATA$dx[i]="CN"
  } else if(METADATA$FJPL163_clinical_num[i]==2){
    METADATA$dx[i]="FEP_1"
  } else if(METADATA$FJPL163_clinical_num[i]==3){
    METADATA$dx[i]="FEP_2"
  }else if(METADATA$FJPL163_clinical_num[i]==4){
    METADATA$dx[i]="FEP_3"
  }

}


# Add "METADATA" to the volume dataframe with the SAME COLUMN NAMES
METADATA<-data.frame(METADATA)

datos <- cbind(age_months=METADATA$age_months,
                age_days=METADATA$age_days,
                sex=METADATA$sex,
                Ventricles=rep(NA,nrow(d)),
                WMV=rep("",nrow(d)),
                sGMV=rep("",nrow(d)),
                GMV=rep("",nrow(d)),
                etiv=rep("",nrow(d)),
                study=METADATA$study,
                fs_version=rep('FS53',nrow(d)), 
                run=rep(1,nrow(d)),
                session=rep(1,nrow(d)),
                dx=METADATA$dx,
                country=rep(NA,nrow(d)),
                participant=METADATA$ID,
                CT=rep("",nrow(d)),
                SA=rep("",nrow(d)),d)
 
datos <- data.frame(datos)

numeric_variables<-c(1,2,5:8,11,12,16:ncol(datos)) 
for (i in numeric_variables){
  datos[[i]] <- as.numeric(datos[[i]]) # convert numeric variables to double
}

# Save the file in .csv format
write.csv(datos,'mean_volumes.csv',quote=FALSE, row.names = FALSE)  

