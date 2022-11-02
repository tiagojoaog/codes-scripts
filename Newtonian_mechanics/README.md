
- Integrators class has Euler, second order Runge-Kutta and fourth order Runge-kutta methods. It also defines all relevant body properties (such as Sun). 
- Motion class uses integrators to calculate Kutta coefficients, adaptive step, energies and updates acceleration,velocities and positions.

## Install
* Pre-requirement 
    * vpython,
    * numpy,
    * keyboard,
    * imageio,
    * matplotlib,
    * jupyter-packaging,
    * hyperlink,
    * cryptography,
    * txaio

* run `python setup.py develop` to get the code installed in develop mode.
* run `python setup.py install` to get the code installed in standard mode.

## Usage

For usage see the provided examples.
