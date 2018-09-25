# FID-A
Toolbox for simulation and processing of in-vivo magnetic resonance spectroscopy data

Getting Started (requires a working internet connection):
1. EITHER a) Go to https://github.com/CIC-methods/FID-A and click on the "Clone or Download" -> "Download ZIP" button to download the package.  Unzip the downloaded "FID-A-master.zip" file and move the resulting folder "FID-A-master" into the directory where you would like to store FID-A and its contents (e.g. /Users/JohnDoe/MatlabStuff/).
   
   OR     b) open a terminal, navigate to the directory where you would like to store FID-A and its contents (e.g. /Users/JohnDoe/MatlabStuff/), and type "git clone https://github.com/CIC-methods/FID-A.git".  This option allows you to easily merge software updates or improvements at a later date using the "git pull" command.  

2. Open MATLAB and add the FID-A directory to your path using: "addpath(genpath('/Users/JohnDoe/MatlabStuff/FID-A-master'));" (if you used option a above) or "addpath(genpath('/Users/JohnDoe/MatlabStuff/FID-A'));" (if you used option b above).

3. Now you are ready to start using the toolbox!

IMPORTANT NOTE:  FID-A contains a folder of example MRS data that is useful for testing purposes.  However, the ONLY way to get this example data is to use option 1b above (The example data will not get downloaded using option 1a).  Also:  In order for option 1b to work properly, you need to have git installed, AND you also need to download and install the 'large file storage' (lfs) extension for git (https://git-lfs.github.com).  Follow the online instructions for installing git-lfs, and once you have successfully installed it, option 1b will properly download the example data (as well as the rest of the FID-A repository).  
