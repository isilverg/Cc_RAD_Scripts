#!/usr/bin/env R
#args = (commandArgs(TRUE))
#print(args[1])
#file = eval(parse(text=args[1]))

#work this out after
library("vroom")
setwd("/gpfs/research/fuenteslab/cm_project/parse/")
tabs<-vroom("all_50.mafs")


####################
#Description:
#Parse file for occurrences of number of individuals, nInd, over a user-defined threshold
#Once this threshold is met, look at the following positions, in order, to define a region
#of the chromosome where the threshold is met. If a region dips under the threshold, it 
#is not included in the output. Additionally, the region must have a length equal
#to the default length, def_length

#send output to file
sink("cm_134_50.txt", split=FALSE, append=FALSE)

####################
# Variables to tune

def_length <- 134
nInd_cutoff <- 50 #threshold


# Start and stop will be based on position of chromosome
start <- -1	#Set default start to -1
stop <- -1	#Set default stop to -1 

# Scaff1 and scaff 2 will be based on name of chromosome	
scaff1 <- ""
scaff2 <- ""

length <- 0;	#length of frame (positions in a row that have a min # of individuals)
ave <- 0;	  #average number of individuals per frame

row_index <- 1
total_rows <- nrow(tabs)- 1
####################

while (row_index <= total_rows){
  #print("aa")
  row <- tabs[row_index,] #Store line in row (a list): index 1-6: chromo, position, major, minor, unknownEM, nInd (number of individuals)
  
  nInd <- row[[6]]
  pos <- row[[2]]
  chromo <- row[[1]]
  
  if (nInd >= nInd_cutoff & start == -1) { # Start of frame: If this specific chromosome position has a min # of individuals and previous position didn't (or was set to default)
    start <- pos;
    stop <- pos;			# Once in this if statement, start is no longer default -1
    
    scaff1 <- chromo
    scaff2 <- chromo
    
    length <- 1; #It is the start of the frame, therefore length is 1
    ave <- nInd; 
    
  } else if (nInd >= nInd_cutoff & start != -1) { #Middle of frame: If this specific chromosome position has a min # of individuals AND is not the first position in a row to have a min # of individuals
    
    stop <- pos
    scaff2 <- chromo
    
    length <- length + 1			# Increase the length of the frame (positions in a row that have a min # of individuals)
    ave <- ave + nInd;		# Sum the nInd together (from frame) for future average calculation 
    
  } else if (nInd <= nInd_cutoff & start != -1) { #If this position is breaks the frame (positions in a row that have a min # of individuals)
    ave <- ave/length	#calculate the average
    
    if ((length == def_length) & (scaff1 == scaff2) & (stop - start + 1) == def_length ) { # if length of frame is the default/wanted length, verified manually with stop-start +1, and on the same chromosome
      cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", scaff1, start, scaff2, stop, length, ave))	 # then print the: chromosome_name	start_position  	chromosome_name	    stop_position   length_of_frame (should be default length)    average_number_of_individuals_sequenced
      
    }
    
    # Reset start and stop positions, chromosome scaffs, and frame stats BUT continue from where you left off/the next line. Therefore you cant have overlapping frames
    start <- -1	
    stop <- -1
    scaff1 <- ""
    scaff2 <- ""
    length <- 0
    ave <- 0
  }
  row_index <- row_index + 1
}

#print last frame if it's still >= threshold (since above it is only printed when the position breaks the frame)

if (start != -1) {
  ave <- ave/length	#calculate the average
  
  if ((length == def_length) & (scaff1 == scaff2) & (stop - start + 1) == def_length ) { # if length of frame is the default/wanted length, verified manually with stop-start +1, and on the same chromosome
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", scaff1, start, scaff2, stop, length, ave))	 # then print the: chromosome_name	start_position  	chromosome_name	    stop_position   length_of_frame (should be default length)    average_number_of_individuals_sequenced
  } 
}
