################################################################################
# Scrpit to plot regional brain maps of: (1) predicted centiles/effect sizes, 
# (2) empirical centiles/effect sizes, and (3) molecular and micro-architectural 
# features from .csv files
################################################################################

#  Copyright (C) 2023 University of Seville
# 
#  Written by Natalia García San Martín (ngarcia1@us.es)
# 
#  This file is part of Neurobiology Centiles Psychosis toolkit.
# 
#  Neurobiology Centiles Mapping toolkit is free software: 
#  you can redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
# 
#  Neurobiology Centiles Mapping toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Neurobiology Centiles Mapping toolkit. If not, see 
#  <https://www.gnu.org/licenses/>.

rm(list=ls()) # Previous data cleaning in memory

# Change file path as desired
location = 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/'
# location = 'D:/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/'
setwd(location)


# Libraries
library(dplyr)
library(ggseg)
library(ggplot2)

############### PREDICTED MAPS ###############
# Load predicted centiles and effect sizes of CCA-PCA diagnoses
Y_pred_dx <-read.csv(paste0(location,"Data/PCA_CCA/Y_pred_dx.csv"), row.names = 1)
Y_pred_dx <- subset(Y_pred_dx, !(rownames(Y_pred_dx) == "FEP_1"))
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "SCZ_Relative"] <- "SCZ-relatives"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "Schizoaffective_Disorder_Relative"] <- "SAD-relatives"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "Schizoaffective_Disorder"] <- "SAD"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "FEP_1"] <- "FEP"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "PE_1"] <- "PE-suspected"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "PE_2"] <- "PE-definite"
rownames(Y_pred_dx)[rownames(Y_pred_dx) == "PE_3"] <- "PE-clinical"

# Load predicted centiles and effect sizes of CCA-PCA groups
Y_pred_groups <-read.csv(paste0(location,"Data/PCA_CCA/Y_pred_groups.csv"), row.names = 1) 
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G0"] <- "SCZ and SAD-relatives"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G1"] <- "PE"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G1_5"] <- "FEP"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2"] <- "SCZ and SAD-chronic"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_without_ASRB"] <- "SCZ and SAD-chronic_without_ASRB"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_without_BSNIP"] <- "SCZ and SAD-chronic_without_BSNIP"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_without_LA5c"] <- "SCZ and SAD-chronic_without_LA5c"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_without_MCIC"] <- "SCZ and SAD-chronic_without_MCIC"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_vs_G1_5"] <- "SCZ and SAD-chronic vs FEP"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_vs_G1"] <- "SCZ and SAD-chronic vs PE"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G2_vs_G0"] <- "SCZ and SAD-chronic vs SCZ and SAD-relatives"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G1_5_vs_G1"] <- "FEP vs PE"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G1_5_vs_G0"] <- "FEP vs SCZ and SAD-relatives"
rownames(Y_pred_groups)[rownames(Y_pred_groups) == "G1_vs_G0"] <- "PE vs SCZ and SAD-relatives"

Y_pred <- rbind(Y_pred_dx,Y_pred_groups) 


############### EMPIRICAL MAPS ###############

dx_total <- list('SCZ-relatives','SAD-relatives','PE-suspected','PE-definite','PE-clinical',
           'FEP','SCZ','SCZ_without_ASRB','SCZ_without_BSNIP','SCZ_without_LA5c','SCZ_without_MCIC',
           'SAD','SCZ_and_SAD-relatives','PE','FEP','SCZ_and_SAD-chronic','SCZ_and_SAD-chronic_without_ASRB',
           'SCZ_and_SAD-chronic_without_BSNIP','SCZ_and_SAD-chronic_without_LA5c','SCZ_and_SAD-chronic_without_MCIC',
           'SCZ_and_SAD-chronic_vs_SCZ_and_SAD-relatives','PE_vs_SCZ_and_SAD-relatives',
           'SCZ_and_SAD-chronic_vs_PE','SCZ_and_SAD-chronic_vs_FEP','FEP_vs_PE','FEP_vs_SCZ_and_SAD-relatives')

for (i_tot in 1:length(dx_total)){
  dx <- dx_total[i_tot]
  print(dx)


if (dx == 'PE' || grepl('and',dx) || grepl('vs',dx)){
  folder_dx <- 'groups'
}else if (grepl('PE',dx)){
  folder_dx <- 'PE'
}else if (grepl('SCZ',dx) || grepl('SAD',dx)){
  folder_dx <- 'SCZ_SAD'
}else if (grepl('FEP',dx)){
  folder_dx <- 'FEP'
}

type <- 'mean' # centiles
route <- paste0(location, "Data/Centiles/",folder_dx,"/",type,"_",dx,".csv")
type <- 'centiles'
if (grepl('vs',dx)){
  type <- 'effsizes'
  route <- paste0(location, "Data/Effect sizes/",type,"_",dx,".csv")
}


# We load empirical centiles/effect sizes
NEWDATA <- read.csv(route)
# We load significant_values
route_significant_values <- paste0(location, "Data/significant_values.csv")
significant_values <- read.csv(route_significant_values)
column_to_plot_significant_values = paste0(type,'_',dx)

file_name <- basename(route)
file_name <- gsub(".csv","",file_name)
names(NEWDATA) <-c('regions',file_name)
column_to_plot = file_name


# We set the list of region of interest names (Desikan-Killiany atlas). The order is crucial, 
# especially for later data representation.
zones <- c( 'bankssts',
            'caudal anterior cingulate',
            'caudal middle frontal',
            'cuneus',
            'entorhinal',
            'frontal pole',
            'fusiform',
            'inferior parietal',
            'inferior temporal',
            'insula',
            'isthmus cingulate',
            'lateral occipital',
            'lateral orbitofrontal',
            'lingual',
            'medial orbitofrontal',
            'middle temporal',
            'paracentral',
            'parahippocampal',
            'pars opercularis',
            'pars orbitalis',
            'pars triangularis',
            'pericalcarine',
            'postcentral',
            'posterior cingulate',
            'precentral',
            'precuneus',
            'rostral anterior cingulate',
            'rostral middle frontal',
            'superior frontal',
            'superior parietal',
            'superior temporal',
            'supramarginal',
            'temporal pole',
            'transverse temporal')


# Formatting data to construct a table with the format required by the library that will represent it.
# The table consists of three columns: region (standardized name of the region), values (numeric values for each region: centiles or eff sizes),
# and groups (subdivisions of the regions. In this case, all 34 regions will always have the same group label.

someData <- tibble(
  region = zones, 
  values = NEWDATA[[(column_to_plot)]],
  groups = c(rep(sub('_',' ',sub(paste(".*", dx, sep = ""), dx, column_to_plot)), 34))          
)

someData_significant_values <- tibble(
  region = zones, 
  !!ifelse(grepl('vs', column_to_plot_significant_values), "values", "values") := significant_values[[gsub('-','.',column_to_plot_significant_values)]],
  groups = c(rep(sub('_',' ',sub(paste(".*", dx, sep = ""), dx, column_to_plot_significant_values)), 34))          
)

# Obtain region coordinates to plot regional names on the map
data(dk)
indices <- which(dk$data$region %in% zones)
total_zones <- dk$data$region[indices]
coord_x <- numeric(length(indices))
coord_y <- numeric(length(indices))
j = 1
for (i in indices){
  total_coord <- unlist(dk$data$geometry[i])
  coord_x[j] <- mean(total_coord[1:length(total_coord)/2])
  coord_y[j] <- mean(total_coord[(length(total_coord)/2 + 1):length(total_coord)])
  j = j +1
}
coord = data.frame(total_zones,coord_x,coord_y)

labels <- na.omit(dk$data$label)
labels <- labels[-c(23,45)] # delete corpus callosum
coord <- cbind(labels,coord)
coord_lh <- coord[1:(nrow(coord)/2),]
coord_lh <- coord_lh[-c(27),] # delete the wrong areas

# Adjust position
coord_lh$coord_x <-c(206,80,270,250,200,275,51,175,73,30,42,150,105,10,70,216,137,187,73,145,95,577,393,537,460,475,310,408,660,520,507,397,335,525,105,426,600,587,311,353,585)
coord_lh$coord_y <-c(86,160,45,117,31,75,53,46,105,63,85,155,145,120,195,185,60,137,15,87,72,118,103,17,25,90,305,53,40,34,176,77,324,140,298,127,62,160,240,345,13)

# Configuration of the desired limits of representation (depends on whether it is for centiles or effect sizes)
limits_matrix <- cbind(significant_values,t(Y_pred))
if (type == 'centiles') {
  limits <- c(min(unlist(select(limits_matrix, -1, -contains('vs'), -contains('FEP1'))), na.rm = TRUE),
              max(unlist(select(limits_matrix, -1, -contains('vs'), -contains('FEP1'))), na.rm = TRUE))
} else if (type == 'effsizes') {
  limits <- c(min(unlist(select(limits_matrix, -1, -!contains('vs'))), na.rm = TRUE),
              max(unlist(select(limits_matrix, -1, -!contains('vs'))), na.rm = TRUE))
}

# Map the data
someData %>%
  group_by(groups) %>%
  ggplot()+ scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = ifelse(type == 'effsizes', 0, ifelse(type == 'centiles', 0.5, NA)),
    labels = function(x) sprintf("%.2f", x),
    #space = "Lab",
    limits = limits
    # limits = c(min(0.46), max(0.53))
    #na.value = "grey50",
    #aesthetics = "fill"
  ) +
  geom_brain(atlas = dk,
             position = position_brain(hemi ~ side),
             show.legend = TRUE,
             aes(fill = values)) +
  guides(fill = guide_colorbar(title = NULL)) + 
  
  # # Uncomment for plotting regional labels
  # Region labels
  # geom_text(data = coord_lh[coord_lh$coord_y<200,],
  #           x = coord_lh[coord_lh$coord_y<200,]$coord_x,
  #           y = coord_lh[coord_lh$coord_y<200,]$coord_y,
  #           label = coord_lh[coord_lh$coord_y<200,]$total_zones,
  #           size = 3,
  #           vjust = 1.5,
  #           color = "black") +
  
  # facet_wrap(~groups)
  labs(title = dx) +
  # theme(plot.title = element_text(hjust = 0.5),
  #       panel.background = element_rect(fill = '#B1C8F5'),
  #       panel.grid = element_blank(),
  #       panel.border = element_blank(),
  #       legend.background = element_rect(fill = "#B1C8F5"))
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())


file_name <- basename(route)
# folder <- basename(dirname(route))
if (grepl('Effect sizes',route)){
  folder <- 'effsizes'
}else{
  folder <- 'centiles'
}

file_name <- gsub(".csv",".png",file_name)

ggsave(file.path(paste0(location,"Data/maps/empirical_maps/",folder,'/'),file_name), width = 8, height = 5, dpi = 300)


someData_significant_values %>%
  group_by(groups) %>%
  ggplot()+ scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = ifelse(type == 'effsizes', 0, ifelse(type == 'centiles', 0.5, NA)),
    labels = function(x) sprintf("%.2f", x),
    #space = "Lab",
    limits = limits
    # limits = c(0, 1)
    #na.value = "grey50",
    #aesthetics = "fill"
  ) +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),show.legend = TRUE,
             aes(fill = !!sym(colnames(someData_significant_values)[2]))) +
  guides(fill = guide_colorbar(title = NULL)) +
  
  # Region labels
  # geom_text(data = coord_lh[coord_lh$coord_y<200,],
  #           x = coord_lh[coord_lh$coord_y<200,]$coord_x,
  #           y = coord_lh[coord_lh$coord_y<200,]$coord_y,
  #           label = coord_lh[coord_lh$coord_y<200,]$total_zones,
  #           size = 3,
  #           vjust = 1.5,
  #           color = "black") +
  
  # facet_wrap(~groups)
  labs(title = dx) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())

ggsave(file.path(paste0(location,"Data/maps/empirical_maps/significant_regions/",folder,'/'),paste0(column_to_plot_significant_values,'.png')), width = 8, height = 5, dpi = 300)


}

############### MOLECULAR and MICRO-ARCHITECTURAL MAPS ###############

molecular_maps <- read.csv(paste0(location,"Data/Molecular_maps/molecular_map_hemi_ordered.csv"), row.names = 1)
colnames(molecular_maps)[7] <- 'α4β2'
molecular_maps_names <- colnames(molecular_maps)
molecular_maps_names[1:6] <- sub('X5','5-',colnames(molecular_maps)[1:6])
molecular_maps_names[12] <- 'GABAa/bz'
molecular_maps_names[23:24] <- list('Neuro-Ex','Neuro-In')
molecular_maps_names[27] <- 'Layer I'
molecular_maps_names[28] <- 'Layer II'
molecular_maps_names[29] <- 'Layer III'
molecular_maps_names[30] <- 'Layer IV'
molecular_maps_names[31] <- 'Layer V'
molecular_maps_names[32] <- 'Layer VI'
molecular_maps_names[33] <- 'gene expression PC1'
molecular_maps_names[34] <- 'T1w:T2w (myelin)'
molecular_maps_names[35] <- 'cortical thickness'
molecular_maps_names[36] <- 'developmental expansion'
molecular_maps_names[37] <- 'evolutionary expansion'
molecular_maps_names[38] <- 'cerebral blood flow'
molecular_maps_names[39] <- 'cerebral blood volume'
molecular_maps_names[40] <- 'cerebral metabolic rate of oxygen'
molecular_maps_names[42] <- 'allometric scaling (NIH)'
molecular_maps_names[43] <- 'allometric scaling (PNC)'
molecular_maps_names[44] <- 'neurotrasmitter PC1'
molecular_maps_names[45] <- 'synapse density'
molecular_maps_names[46] <- 'glycolytic index'
  

for (i in 1:ncol(molecular_maps)){
  
  someData_molecular_maps <- tibble(
    region = zones, 
    values = molecular_maps[,i],
    groups =  c(rep(colnames(molecular_maps)[i], 34))           
  )
  
  someData_molecular_maps %>%
    group_by(groups) %>%
    ggplot()+ scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      # midpoint = 0.5,
      labels = function(x) sprintf("%.2f", x),
      #space = "Lab",
      limits = c(min(molecular_maps),max(molecular_maps))
      # limits = c(0, 1)
      #na.value = "grey50",
      #aesthetics = "fill"
    ) +
    geom_brain(atlas = dk, 
               position = position_brain(hemi ~ side),show.legend = TRUE,
               aes(fill = !!sym(colnames(someData_molecular_maps)[2]))) +
    guides(fill = guide_colorbar(title = NULL)) +
    
    # facet_wrap(~groups) +
    labs(title = molecular_maps_names[i]) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank())
    
  ggsave(file.path(paste0(location,"Data/maps/molecular_maps"),paste0(gsub('\\.','_',colnames(molecular_maps)[i]),'.png')), width = 8, height = 5, dpi = 300)
  
}

############### PREDICTED MAPS ###############

for (i in row.names(Y_pred)){
  print(i)
  someData_predicted <- tibble(
    region = zones, 
    values = t(Y_pred[rownames(Y_pred) == i,]),
    groups = c(rep(i, 34))          
  )
  
  if (!grepl('vs',i)) {
    limits <- c(min(Y_pred[setdiff(1:nrow(Y_pred),grep("vs", rownames(Y_pred))), ]),
                max(Y_pred[setdiff(1:nrow(Y_pred),grep("vs", rownames(Y_pred))), ]))
  } else if (grepl('vs',i)) {
    limits <- c(min(Y_pred[grep("vs", rownames(Y_pred)), ]),
                max(Y_pred[grep("vs", rownames(Y_pred)), ]))
  }
  
  someData_predicted %>%
    group_by(groups) %>%
    ggplot()+ scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = ifelse(grepl('vs',i), 0, ifelse(!grepl('vs',i), 0.5, NA)),
      labels = function(x) sprintf("%.2f", x),
      #space = "Lab",
      limits = limits
      # limits = c(min(-0.81), max(0.4))
      #na.value = "grey50",
      #aesthetics = "fill"
    ) +
    geom_brain(atlas = dk,
               position = position_brain(hemi ~ side),
               show.legend = TRUE,
               aes(fill = values)) +
    guides(fill = guide_colorbar(title = NULL)) + 
    
    # facet_wrap(~groups)
  labs(title = i) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())
  
  ggsave(file.path(paste0(location,"Data/maps/predicted_maps"),paste0(gsub(' ','_',i),'.png')), width = 8, height = 5, dpi = 300)
  
}

'Done!'


