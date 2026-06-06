library(matrixStats)
library(clue)

# Prep coord files and try running rotate.parcellation function
source('rotate_parcellation.r')

LH_coords_Schaefer100 <- "LH_Schaefer100_Coords.txt"


rotate.parcellation(coord.l=LH_coords_Schaefer100,
                    nrot=1,
                    method='hungarian')
