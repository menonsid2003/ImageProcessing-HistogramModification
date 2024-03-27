# ImageProcessing-HistogramModification
This project can take a greyscale and color image and perform a series of modifications. This primarily includes histogram stretching in order to improve contrast. The user can select to only perform this on a specific color channel, all 3 (for ppms), or on simply a grey scale image (pgm). 

Do not delete any of the folders, they are required as output paths.
The format of the parameters.txt file is as follows:

input.pgm output.pgm(for task 1 from last hw, does not matter here, but it is still required)
commandName(eg: allstretch slidingwindow colorstretch)sub roiX roiY roiSize A B
allstretch roiX roiY roiSize A B
colorstretch filename.ppm colorChannel(R,G,B,A) roiX roiY roiSize A B

lines beginning with '// ' will be read as a comment and subsequently ignored by the program.

example:
men.pgm task1_out.pgm
// allstretch 50 20 200 0 255
allstretch 100 200 200 0 250
slidingwindow 0 0 500 30 150 5
// colorstretch valtr.ppm R 0 0 2000 0 255
// colorstretch valtr.ppm B 0 0 2000 0 255
colorstretch cat.ppm A 50 100 600 40 220

The outputs will be saved into folders for organization.

compile using:
g++ -std=c++11 main.cpp
and run with:
./a.out or ./a depending on if macOS or Windows respectively.
do not need to recompile if you want to change parameters, just run ./a again.

Test images I used are included: men.pgm, baboon.pgm, valtr.ppm, cat.ppm
