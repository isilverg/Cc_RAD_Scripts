#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

# Variables to tune
$def_length = 134	
$nInd_cutoff = 50; #min # of individuals

# Start and stop will be based on position of chromosome
$start = -1;	#Set default start to -1
$stop = -1;	#Set default stop to -1 

# Scaff1 and scaff 2 will be based on name of chromosome	
$scaff1 = "";
$scaff2 = "";

$length = 0;	#length of frame (positions in a row that have a min # of individuals)
$ave = 0;	#average

 
while (<FILE>) {	#Read each line

	$line = $_; 			# Store line from loop in variable line
	chomp($line);			# Trim line
	@tabs = split(/\t/,$line);	# Separate line, by tabs, into a list
					# tabs index 0-5: chromo, position, major, minor, unknownEM, nInd (number of individuals)
					
	$nInd= $tabs[5];		
	$pos = $tabs[1];
	$chromo = $tabs[0];

	# Default variables start, stop, scaff1, scaff2 don't change until the first position where the min # of individuals is achieved
	

	if ($nInd >= $nInd_cutoff && $start == -1) { # Start of frame: If this specific chromosome position has a min # of individuals and previous position didn't (or was set to default)
		$start = $pos;
		$stop = $pos;			# Once in this, start is no longer default -1

		$scaff1 = $chromo;
		$scaff2 = $chromo;

		$length = 1;
		$ave =nInd; 

	} elsif ($nInd >= $nInd_cutoff && $start != -1) { #Middle of frame: If this specific chromosome position has a min # of individuals AND is not the first position in a row to have a min # of individuals
		
		$stop = $pos;
		$scaff2 = $chromo;

		$length++;			# Increase the length of the frame (positions in a row that have a min # of individuals)
		$ave = $ave + $nInd;		# Sum the nInd together (from frame) for future average calculation 

	} elsif ($nInd <= $nInd_cutoff && $start != -1) { #If this position is breaks the frame (positions in a row that have a min # of individuals)
		$ave = $ave/$length;	#calculate the average

		if ($length == def_length && $scaff1 eq $scaff2 && ($stop - $start + 1) == def_length) { # if length of frame is the default/wanted length, verified manually with stop-start +1, and on the same chromosome
			print "$scaff1\t$start\t$scaff2\t$stop\t$length\t$ave\n";			 # then print the: chromosome_name	start_position  	chromosome_name	    stop_position   length_of_frame (should be default length)    average_number_of_individuals_sequenced
		}

		$start = -1;	#Reset start and stop positions, chromosome scaffs, and frame stats BUT continue from where you left off/the next line. Therefore you cant have overlapping frames
		$stop = -1;
		$scaff1 = "";
		$scaff2 = "";
		$length = 0;
		$ave = 0;
	}
}
close FILE;


