#!/usr/bin/env R
#options(echo=FALSE)
#args = commandArgs(trailingOnly=TRUE)

#file<-args[1]
library("vroom")
input_file<-vroom("all_60.mafs")


# Expected input file format:
# chromo  position        major   minor   unknownEM       nInd
# QVIC01000001.1  2969    T       A       0.000002        1
# QVIC01000001.1  2970    C       A       0.000002        1
# QVIC01000001.1  2971    C       A       0.000002        1
# QVIC01000001.1  2972    C       A       0.000002        1
# QVIC01000001.1  2973    A       C       0.000002        1
# QVIC01000001.1  2974    G       A       0.000002        1
# QVIC01000001.1  2975    G       A       0.000002        1
# QVIC01000001.1  2976    G       A       0.000002        1
# QVIC01000001.1  2977    C       A       0.000002        1


####################
#Description:
#Parse file for occurrences of number of individuals, nInd, over a user-defined threshold
#Once this threshold is met, look at the following positions, in order, to define a region
#of the chromosome where the threshold is met. If a region dips under the threshold, it
#is not included in the output. Additionally, the region must have a length greater than
#or equal to the default length, def_length and less than or equal to the upper length,
#end_length
#
#Examples:
#If def_length is 134 and end_length is 300, then a continious region where nInd is above
#the threshold nInd_cutoff will only be reported if the region's length is >= 134 and <= 300.
#Therefore if the continious region above the threshold is 301 bases long, it is not reported.
#Also, if the region 267 bases long and meeets the threshold on all bases but #134, neither
#the first half (133 bases long) nor the second half (133 bases long) will meet the criteria
#and that region will not be reported.


####################
# Variables to tune

def_length <- 134
end_length <- 300
nInd_cutoff <- 70 #threshold


# Start and stop will be based on position of chromosome
start <- -1     #Set default start to -1
stop <- -1      #Set default stop to -1

# Scaff1 and scaff 2 will be based on name of chromosome
scaff1 <- ""
scaff2 <- ""

length <- 0;    #length of frame (positions in a row that have a min # of individuals)
ave <- 0;         #average number of individuals per frame

row_index <- 1
total_rows <- nrow(input_file)- 1
####################

while (row_index <= total_rows){ #Iterate through all the rows in input_file
  
  row <- input_file[row_index,] #Store entire row
  
  # row format:
  # chromo  position        major   minor   unknownEM       nInd(number of individuals)
  nInd <- row[[6]]
  pos <- row[[2]]
  chromo <- row[[1]]
  
  if (nInd >= nInd_cutoff & start == -1) { # Start of frame: If this specific chromosome position has a min # of individuals and previous position didn't (or was set to default)
    start <- pos;
    stop <- pos;                        # Once in this if statement, start is no longer default -1
    
    scaff1 <- chromo
    scaff2 <- chromo
    
    length <- 1; #It is the start of the frame, therefore length is 1
    ave <- nInd;
    
  } else if (nInd >= nInd_cutoff & start != -1) { #Middle of frame: If this specific chromosome position has a min # of individuals AND is not the first position in a row to have a min # of individuals
    
    stop <- pos
    scaff2 <- chromo
    
    length <- length + 1                        # Increase the length of the frame (positions in a row that have a min # of individuals)
    ave <- ave + nInd;          # Sum the nInd together (from frame) for future average calculation
    
  } else if (nInd <= nInd_cutoff & start != -1) { #If this position is breaks the frame (Used to have continious positions in a row that have a min # of individuals but this position does not meet the threshold)
    ave <- ave/length   #calculate the average
    
    if ((length >= def_length) & (length <= end_length) & (scaff1 == scaff2) & ((stop - start + 1) >= def_length)) { # if length of frame is the default/wanted length, verified manually with stop-start +1, and on the same chromosome
      cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", scaff1, start, scaff2, stop, length, ave))         # then print the: chromosome_name      start_position          chromosome_name     stop_position   length_of_frame (should be default length)    average_number_of_individuals_sequenced
      
    }
    
    # Reset start and stop positions, chromosome scaffs, and frame stats BUT continue from where you left off/the next line. Therefore you cant have overlapping frames
    start <- -1
    stop <- -1
    scaff1 <- ""
    scaff2 <- ""
    length <- 0
    ave <- 0
  }
  row_index <- row_index + 1 #Continue to next row regardless if this current position meets the # of individuals threshold
}

#Print last region if it is valid (would not be printed in the above loop if the last position did not dip under the # of individuals threshold)
if (start != -1) {
  ave <- ave/length     #calculate the average
  
  if ((length >= def_length) & (length <= end_length) & (scaff1 == scaff2) & ((stop - start + 1) >= def_length)) { # if length of frame is the default/wanted length, verified manually with stop-start +1, and on the same chromosome
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", scaff1, start, scaff2, stop, length, ave))   # then print the: chromosome_name      start_position          chromosome_name     stop_position   length_of_frame (should be default length)    average_number_of_individuals_sequenced
  }
}
