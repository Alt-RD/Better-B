# Better-B project: thermal modeling of beehives

This repository contains the hygro-thermal modeling tools of beehives carried out by [Alt-RD research group](https://www.alt-rd.com/) in the [European project Better-B](https://www.better-b.eu/).  

This work was supported by the Better-B project, which has received funding from the European Union, the Swiss State Secretariat for Education, Research and Innovation (SERI) and UK Research and Innovation (UKRI) under the UK government's Horizon Europe funding guarantee (grant number 10068544). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union, European Research Executive Agency (REA), SERI or UKRI. Neither the European Union nor the granting authorities can be held responsible for them.

The GitHub repository currently contains about 90 source files written in [Octave programming language](https://www.octave.org/). These files are named and grouped into directories based on the operations they perform. There are two main directories at the root of the repository:
1. The “HiveModel” directory contains the files related to the hive hygro-thermal models (about 20 files) implemented based on the HiveTemp library.
2. The “HiveTemp” directory contains all files related to the HiveTemp library (about 70 files).
The HiveTemp library is a set of functions which perform basic operations. They must be combined to define the geometric shape of the hive and to solve the heat and mass transfer equations. To keep the HiveTemp library general and independent of any specific application, it was chosen to separate the hive hygro-thermal models from of the HiveTemp library.  

Since both the hive models and the HiveTemp library are in full development, no comprehensive documentation has been written yet. Instead, it was preferred to add a brief documentation inside the source files, below the licence header. This approach allows easy update of the documentation and prevent mismatches with the current source code version that would be misleading for the user.


## HiveModel
This directory contains the hygro-thermal models of three types of hive:
1. Dadant hive
2. Log hive
3. Ecological hive

![HiveImage](https://github.com/user-attachments/assets/d42927de-d909-455f-8b64-32b35f0b57d4)

## HiveTemp library
The HiveTemp library contains functions to implement hygro-thermal models of beehive. Although this library has been developed by Alt-RD for thermal modelling of hives, it is not restricted to this field. 
The general principle of this library is to define complex shapes (like beehives) by defining and connecting elementary shapes like rectangular cuboids. A space discretization using the Finite Volume Method (FVM) [^1] is then applied to each elementary volume to get a matrix formulation of the heat and mass transfer equations. The matrices of all elementary volumes are finally merged into larger matrices to get one single matrix equation which is solved numerically by applying a time discretization. 

![HiveTempPrinciple](https://github.com/user-attachments/assets/827dc8c5-7591-4b3e-a41d-af6264992334)

[^1]: Eymard R., Gallouet T., Herbin R., The Finite Volume Method, Handbook of Numerical Analysis, North Holland, 2000, P. Ciarlet et J.L. Lions, 713-1020
