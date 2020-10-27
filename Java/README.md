# FileProcess

  Java file processor to take the output of Hybond and collate it into a single .csv file for interpretation.

__Usage:__ 
+ Argument 1 = filename without the frame number in it
+ Argument 2 = Number of files that need to be processed (number of frames)
+ Argument 3 = Line number within the .csv files being processed that holds the value wanted
+ Argument 4 = Name of the output file
+ Argument 5 = the type of file being read (options are)  
  * m - minus interaction (attractive)
  * p - positive interaction (repulsive)               
  * t - total interaction (sum of previous files)
