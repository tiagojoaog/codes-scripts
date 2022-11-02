To make a GIF, first load imageIO module by typing 'module load Python/Imageio/3/Imageio-2.16.2' on the terminal

Then type on the terminal:

ase_make_gifs name=relax.traj angle=-90 direction=x savename=movie.gif interval=1

Name - The name of the animated trajectory. Could be relax.traj, vib.0.traj XDATCAR etc...
Angle - The rotation angle of the movie on a given direction compared to the standard view.
Direction - The direction where the rotation is taking place
Savename - The name of the output GIF
interval - Provides a GIF with frames=(total images)/interval.        
