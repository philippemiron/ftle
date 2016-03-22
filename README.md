ftle

Copyright (C) 2014  Philippe Miron

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses


Finite time Lyapunov exponent (FTLE) from analytical velocity field

-  velocity field expression parsing (string in config file)
	- done using the beautiful tool exprtk, more info here https://code.google.com/p/exprtk/

- the method is explain in multiple article from  G. Haller, F. Lekien. and S.C. Shadden.
	- http://en.wikipedia.org/wiki/Lyapunov_exponent
	- http://mmae.iit.edu/shadden/LCS-tutorial/FTLE-derivation.html

Usage:
1. `make with the Makefile or with CLion or other IDE`
2. `./main/ftle.bin config.txt  (or other configs in /run/)`
3. `Output file (ftle.dat) will be created in running directory`
4. `Open in Tecplot to visualize FTLE field`
