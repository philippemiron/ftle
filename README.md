ftle

Finite time Lyapunov exponent (FTLE) from analytical velocity field

-  velocity field expression parsing (string in config file)
	- using the beautiful tool exprtk https://code.google.com/p/exprtk/

- the method is explained in multiple papers from  G. Haller, F. Lekien. and S.C. Shadden.
	- http://en.wikipedia.org/wiki/Lyapunov_exponent
	- http://mmae.iit.edu/shadden/LCS-tutorial/FTLE-derivation.html

Usage:
1. `make with the Makefile or with CLion or other IDE`
2. `./main/ftle.bin config.txt  (or other configs in /run/)`
3. `Output file (ftle.dat) will be created in running directory`
4. `Open in Tecplot to visualize FTLE field`
