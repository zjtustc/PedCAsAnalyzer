Args <- commandArgs(TRUE);

Chr <- Args[1]
gbandFile <- Args[2]
parent_info_het_file <- Args[3]
parent_info_hom_file <- Args[4]
parent_info_ref_file <- Args[5]
parent_fig_out <- Args[6]

plot.chromosome <- function(Chr, gbandFile, parent_info_het_file, parent_info_hom_file, parent_info_ref_file, parent_fig_out){ # Chr == X
  #read in data1
  gband <- read.table(file = gbandFile, head = F, sep = '\t')
  Het <- read.table(file=parent_info_het_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  Hom <- read.table(file=parent_info_hom_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  Ref <- read.table(file=parent_info_ref_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  # statistic
  parent_info_marker_num <- nrow(Het) + nrow(Hom) + nrow(Ref)
  parent_info_het_per <- nrow(Het)/parent_info_marker_num
  parent_info_het_per_new <- paste('n=', nrow(Het), '; percent:', signif(parent_info_het_per,4)*100, '%', spe='')
  parent_info_hom_per <- nrow(Hom)/parent_info_marker_num
  parent_info_hom_per_new <- paste('n=', nrow(Hom), '; percent:', signif(parent_info_hom_per,4)*100, '%', spe='')
  parent_info_ref_per <- nrow(Ref)/parent_info_marker_num
  parent_info_ref_per_new <- paste('n=', nrow(Ref), '; percent:', signif(parent_info_ref_per,4)*100, '%', spe='')

  #read in data2 of cytoband
  df.peak = gband
  chrom_start <- 0
  chrom_end <- max(df.peak[df.peak$V1 == Chr, ]$V3)
  
  bed_table <- df.peak[df.peak[, 1] == Chr, ]

  track_start <- bed_table[, 2]
  track_end <- bed_table[, 3]
  track_value <- bed_table[, 5]
  track_info <- bed_table[, 6]
  
  ## use 255 gray pattern
  alpha_vector <- floor(track_value / 600 * 255) 
  alpha_vector[alpha_vector > 255] = 255
  color_vector <- rgb(t(col2rgb("black")), alpha = alpha_vector, maxColorValue = 255)
  
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  track_y1 <- rep(0.025, nrow(bed_table))
  track_y2 <- rep(0.055, nrow(bed_table))
  
  track_y1[track_info == "acen"] <- 0.035
  track_y2[track_info == "acen"] <- 0.045

  #plot the image
  ## set the imag parameters
  png(parent_fig_out, 
      width = 4, height = 4, units = "in", 
      res=800, pointsize = 6)
  par(
    mar = c(5,5,2,2), 
    xaxs = "i", yaxs = "i",
    cex.axis = 1.7, cex.lab = 2 
  )

  plot(x=c(0,chrom_end/1e6), y = c(-0.1,-0.1), xlab = paste("Chr", Chr, ", Location (Mb)", spe=''), ylab = "Mutation Type", ylim = c(0,0.7), yaxt = "n")
  axis(2, at=c(0, 0.2, 0.4, 0.6), labels=c('', 'Het', 'Hom', 'Ref'), cex.axis = 1.7)
  points(x = Het$Start/1e6,y = rep(c(0.2), each = nrow(Het)),pch =3,cex=1.2, col = rgb(205, 149, 12, 180, maxColorValue = 255))
  text(5, 0.23, labels=parent_info_het_per_new, cex=1.5, adj = 0)
  points(x = Hom$Start/1e6, y = rep(c(0.4), each = nrow(Hom)), pch=3, cex=1.2, col = rgb(0, 0, 255, 180, maxColorValue = 255))
  text(5, 0.43, labels=parent_info_hom_per_new, cex=1.5, adj = 0)
  points(x = Ref$Start/1e6, y = rep(c(0.6), each = nrow(Ref)), pch=3, cex=1.2, col = rgb(0, 139, 139, 180, maxColorValue = 255))
  text(5, 0.63, labels=parent_info_ref_per_new, cex=1.5, adj = 0)

  rect(track_x1/1e6,track_y1,track_x2/1e6,track_y2,col = color_vector,border = T,lwd = 1)
  dev.off()
}

plot.chromosome(Chr, gbandFile, parent_info_het_file, parent_info_hom_file, parent_info_ref_file, parent_fig_out)
