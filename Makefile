platform=$(shell uname)
arch=$(shell uname -m)

# define the right os
# essentially for mkl path
# other library problems..
ifeq ($(platform),Linux)
makefile=Makefile_linux
else
makefile=Makefile_macos
endif

defrule:
	make -f $(makefile)

clean:
	make -f $(makefile) clean
