# dec. 22,2017. Now uses Matlab GUI to organize/modify/display simulations.
# mcxyz add Ryx, use ../sims/ folder for storing input/output files.
All input files (*_H.mci, *_T.bin) and output files (*_F.bin, *_Ryx.bin, *_Rd.dat)
are stored in folder ../simsLibrary/*, where * denotes the name of the run, eg., skinvessel,
which creates skinvessel_H.mci, skinvessel_T.bin, etc.

mcxyz is a 3D multi-voxel Monte Carlo simulation of light transport within biological tissue.
The program is written in ANSI standard C.
It reads input files:
  yourname.mci, which controls the Monte Carlo simulation,
  yourname_T.bin, which assigns tissue types to the voxels.
It outputs:
  yourname_F.bin, which is the output of fluence [J/cm^2/J.delivered],
  yourname_Ryx.bin, which is the escaping flux across the tissue's top surface (zsurf),
  yourname_Rd.dat, which is a single value for the total diffuse reflectance.
While any program can create the input files, this .git uses a MATLAB GUI
to accomplish this task.
The program mcxyz.c is described at http://omlc.org/software/mc/mcxyz/index.html .
