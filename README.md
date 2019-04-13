# Texture Synthesis

This an implementation in C++ of the paper : [Parallel Controllable Texture Synthesis](http://hhoppe.com/paratexsyn.pdf).
The goal is to design an algorithm that given a small exemplar image produce a big image with the same structure as the small image.
The created image has to be different from the input image and the output should also not repeat.

The code is parallelized with OpenMP. To install it, use:
```
sudo apt install libomp-dev
```
We also use the library boost to manage directories and file. To install it, use:
```
sudo apt install libboost-filesystem-dev
```
Finally we use gtkmm for the GUI. To install it, use:
```
sudo apt install libgtkmm-3.0-dev 
```

You can compile the code with the Makefile by typing
```
make
```

Then to execute the algorithm you have to type:
```
./main <filename> [-c] [-t]
```
The first parameter is the name of the image given in input.
The argument `-c` is optional and if it is given then some images called coherence maps are computed. It may take a while. The argument `-t` is optional and if it is given then the image in input is torified.

Here are two examples of images created by the algorithm :  

![example1](ims/2.png) ![result1](ims/2/res.png)

![example1](ims/1.png) ![result1](ims/1/res.png)

## GUI

You can also use a GUI. This GUI allow you to change parameters and see the result in the window in real time. There are some BUGs when the sample is big and generated files are big. To compile this GUI you have to run:
```
make g
```
Then to execute the GUI you have to run:
```
./gui [filename]
```
Where filename is the name of the file you want to load. This parameter is optional and the file can be loaded in the GUI.

## More Images

You can obtain more images with their coherence and sometimes with the result by running:
```
cd ims
./get_more_images.sh
```