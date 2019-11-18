
# PURPOSES OF THE PROJECT:



In this Project our purpose is that of performing a Monte Carlo simulation of the motion of a given particle on a double well potential under the action of white noise.

We have a single particle living in a potential V(x) of the form:

V(x) = (x^2 - ax) (x^2 + bx)

And want to perform a Monte Carlo simulation integrating the motion with a first order Symplectic Integrator of equations:

p(i+1) = p(i) - gamma*p(i)*dt - V(i)''*dt + eps*csi(i)*sqrt(dt)
x(i+1) = x(i) + p(i+1)*dt

In order to compare numerical results (such as spatial and impulse distributions and variations of probabilities in time) with the theoretical ones. 

All variables explanations are contained in **params.txt**










## HOW TO INSTALL Kramers:



Simply run on your Terminal, in the Directory you want:
'>git clone https://github.com/Giomu/Kramers.git   


Please, before running Kramers, check that you have already installed the following libraries. For any information on these libraries we refer to the official python's documentations as follows:



> random       : https://docs.python.org/3.6/library/random.html

> numpy        : https://numpy.org/doc/1.17/user/setting-up.html

> logging      : https://docs.python.org/3/library/logging.html

> configparser : https://docs.python.org/3/library/configparser.html

> sys          : https://docs.python.org/3/library/sys.html

> scipy        : https://www.scipy.org/install.html










## HOW TO USE Kramers:



1. Open **params.txt** and fill it with the input values required. Please respect the indications given. All parameters explanations will be contained in this file.

2. Run 'Kramers.py'. At the end of the simulation five new txt files will be saved in the Directory used:

	- **Posizioni.txt**       : txt file containing the positions of all the particles for all times
	- **Impulsi.txt**          : txt file containing the values of the impulses for all the particles at any time
	- **LeftFraction.txt**   : txt file containing the numbers of particles contained in the Left hole of the Potential
	- **RightFraction.txt** : txt file containing the numbers of particles contained in the Right hole of the Potential
	- **Time.txt**              : txt file containing all of the time instants 

And a print message will advise you that the simulation is finished asking to run **Graphic.py** .

3. Run **Graphic.py** in order to obtain a visualization of the numerical results previously obtained. Five plots will be showed. The first two will report the spatial and the impulse distributions, the second two will represent the variations of the numerical probabilities compared with the theoretical ones, while the last histogram will compare these probabilities at the end of the simulation.


The file **Kram_Functions.py** is a file containing all the functions created for **Kramers.py**, while **Test_Kram_Functions.py** is a file containing a few examples of tests performed on **Kram_Functions.py** in order to test the robustness and the correctness of the functions created 







