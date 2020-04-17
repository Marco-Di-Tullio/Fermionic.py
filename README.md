# Fermionic operators

This programs numerically defines fermionic operators through a sparse matrix representation. The only input the program needs is the dimension, and the corresponding one body fermionic operators will be automatically defined. The you can define states in the base of posible states of the specified dimension, and calculate some of its properties.


It's done in Python, Mathematica and Julia. In the examples folder there is a comparison in the efficency between the 3 languages. Julia is more than 100 times faster than Python, and more than 10 times faster than Mathematica, being hence the clear winner for this task.

In Python, the class State() gives access to several properties of a given state, such as the one body matrix and the one body entropy.


![](/images/quantuminfo.png)
