# Hybond
---------
__Hydrogen Bond Strength Analysis Tool__

This tool has been developed to provide an insight into a hydrogen bond interaction over a period of time.
Using a trajectory file for a MD simulation this program can extract the frames and characterise a hydrogen bond
over every frame. 

__USAGE:__

_Hybond only:_
To run the program use ./Hybond. This will start the prompts for the Hybond program. The prompts are as follows:
+ Skip the file process - only skip if the frames and the distance file of the target atoms has already been made
+ Name of .dcd file - file path to the .dcd file that is being analysed
+ First atom of analysis - the pdb index of the first atom being looked at
+ Second atom of analysis - the pdb index of the second atom being looked at
+ Name of the xyz file - file path to the .xyz file that holds the atom types (if pdb format will need to use OpenBabel to convert)

__After the frames are finished the program will send another prompt__
+ Name of the distance file - output file name for the distances between target atoms

_Bonder only:_
To run Bonder additional arguments are required when the program is called.

+ Call the program normally with ./Hybond but dont start running yet
+ The program now needs the Bonder arguments from the command line
+ The first argument is the type of interaction - line, trig, quad, all _(only the line method is currently implemented)_
+ The next arguments use the typical Bonder flags which are summarised below
         
_Both Programs:_
To run both it is as simple as adding the arguments to the ./Hybond command line call for Bonder.
Then following all the prompts from the Hybond program to process the .dcd and .xyz file.
+ A typical call to the program would be: ./Hybond line -1 2343 -2 3444 -o unknownBond -c 0.3 -r 0.05
Note: values of -c and -r might need to be evaluated and changed for best results
         
__Bonder Flags/Arguments:__ 
+      -i is the input file to be used by the program _(Hybond doesnt need this argument)_
+      -r is the resolution (sixe of a voxel in the intergration grid)
+      -c is the RDG cut-off value and specifies the maximum value for the RDG
+      -o is the output file (Hybond makes use of this argument)
+      -q is the cubesize which determines the resoltion of the cubefiles produced 
+      -1 is the first atom for use in the line, trig or quad modes
+      -2 is the second atom for use in the line, trig or quad modes
+      -3 is the thrid atom for use in the trig or quad modes _(Hybond does not have functionality for this)_
+      -4 is the fourth atom for use in quad mode only _(Hybond does not have functionality for this)_

Visit the Bonder GitHub for a more detailed overview of the program and installation requirements.
[Bonder GitHub repository](https://github.com/danielaschipper/Bonder-promolecular "Bonder Promolecular GitHub page")
