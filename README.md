# mrRevBayes
translator from MrBayes to RevBayes

### Warnings
Do double check the algorithm. Bayes blocks, as I have seen, are written so differently from one another. Most of the translations I have made using my code has been correct, but a few had small mistakes. Just check the config files made by ```mrRevBayes()``` before moving onto ```mrRevBayes2()```.

### Translating MrBayes to RevBayes
This R library translates Bayes blocks from MrBayes into RevBayes scripts. This is, of course, important if you're interested in rerunning analyses that you've previously done in MrBayes in RevBayes.

### Progress
I have two versions of ```mrRevBayes()``` which is the function that takes a MrBayes Bayes block and turns it into a configuration file for RevBayes. The first version works quite robustly. I have tested it on about 40 datasets. The second version does not work as well, but has the structure to work better eventually (once the bugs are exterminated). For now, I am still developing off of the first version of ```mrRevBayes()``` and will fix the second version (1.1) once I have updated ```mrRevBayes2()``` which creates the RevBayes script from the configuration file.

## Two Functions
```mrRevBayes()```
  - Input: Nexus file with a Bayes Block
  - Ouput: A configuration file for RevBayes
 
 ```mrRevBayes2()```
  - Input: Configuration file (.csv)
  - Output: RevBayes script
  
## Bugs in mrRevBayes() -- First Version
This function has a bug in the loops below ```readLsets()```. For most studies, it works well, but there is a problem when we go through this process:

  1. Unlink all
  2. Link some back together
  
Not all the numbers are unique, and therefore, some partitions that are not supposed to be linked are specified as linked.

## Current Development
I am working on ```mrRevBayes2()``` to write scripts for datasets that are not saturated-partitioned.
