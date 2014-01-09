#  CVS information:
#  $Revision: 1.10 $
#  $Date: 2002/06/18 04:33:40 $
#  $Author: rohl $

# NNMAKE Makefile
# modified from JJG's ROSETTA Makefile  4/10/01         
# 7/3/01 added comments, new flags as in ROSETTA makefile

# requires makefile, make.system, and make.deps.sh

# Included in the package: 
#
# --Automatic dependency file generation, provided by make.deps.sh. 
#   If you change a header file (e.g. misc.h), all object files that
#   depend on this header will be re-compiled
#
# --Automatic compiler choice, provided by make.system depending on the 
#   machine you are on, or you can set the COMPILER variable yourself
#   before invoking make (see functions below)
#
# --Automatic compile flags defined in make.system, for various 
#   types of compilation.
#   Available options:
#
#		make optim        (default)
#		make debug
#		make debugoptim   (only available with gnu f77)
#		make profile      (not yet tested)
#
#   Be sure to 'make clean' when changing from one to another.
#
# --Object files are named with their compiler type, so you can work
#   with more than one compiler at once.  Also, object files are now
#   hidden 'dot' files, e.g. .gnu.main_rosetta.o.
#


# These function definitions in your .bashrc are convenient to use on mary:
#
#    function mkab  () { setenv COMPILER absoft;make $@ }
#    function mkpgi () { setenv COMPILER pgi;make $@ }




#specify COMPILER=gnu, lnx, absoft, alpha, or suse; or let make.system guess:
#COMPILER=absoft

include make.system   # compiler choice and flag settings

# define base name for executable:
BASE_NAME=pNNMAKESHORT

##################### List all Rosetta source files #####################

SRC = \
make_ss_nn.f            all_ss.f                move_frag.f\
dme_native.f            score_fragment.f        dme_score.f\
make_dat.f              read_phd.f              write_names.f\
homologs.f              input_pdb_nn.f          find_nn.f\
read_vall.f             convert1_num.f          \
convert3_num.f          makeoutput.f            main.f\
lower2upper.f           sizemer_dme.f           homologs_nr.f   \
read_rdb.f              input_checkpoint.f      read_jufo.f\
read_jones.f            read_dsc.f              setup.f     \
file_manager.f          zscore_dme.f            loops_nn.f\
constraints_nn.f        dipolar_nn.f            dipolar.f
 

# add the compiler name onto each object file name
OBJS=$(patsubst %.f, .$(COMPILER).%.o, $(SRC))
# creates a list of dependancy files to create/include
DEPS=$(patsubst %.f, .%.d, $(SRC))


# targets that add special case flags
optim: FFLAGS+=$(FOPTIMFLAGS)
optim: compile

debug: FFLAGS+=$(FDEBUGFLAGS)
debug: compile

# warning: with any non-gnu compiler, the debugoptim executable will fail
debugoptim: FFLAGS+=$(FOPTIMFLAGS) $(FDEBUGFLAGS) 
debugoptim: compile

profile: FFLAGS+=$(FPROFILEFLAGS)  # not yet supported
profile: compile

regular: compile


##################### compilation rules #####################

# rule to compile executable
compile: print $(BASE_NAME).$(COMPILER)
$(BASE_NAME).$(COMPILER) :  ${OBJS}
	$(F77) $(FFLAGS) -o $@ $(OBJS) $(LINKFLAGS)

# rule to compile object files:
.$(COMPILER).%.o: %.f
	$(F77) $(FFLAGS) -c -o $@ $<

# rule to make depency files:
.%.d: %.f
	@./make.deps.sh $< $@ $(<:%.f=%.o) 



##################### miscellaneous rules #####################

# all these targets are fake (they don't correspond to a filename)
.PHONY: clean superclean print tags compile optim debug debugoptim profile benchmarks

# set up the TAGS file for jumping to subroutines in Emacs 
# (see [f1],[f2],[f3] in ~jeff/.emacs)
tags:
	etags --language=fortran *h *f


# clean-up commands

clean:
	rm -f .*.d .$(COMPILER).*.o $(BASE_NAME).$(COMPILER) core .ftnchek*

superclean: clean
	rm -f .*.*.o $(BASE_NAME).* 

almostsuperclean:
	rm -f .*.*.o .*.d core .ftnchek*



print:
	$(system-print)

# ftnchek is a utility to check the syntax of fortran files.  it might help
# you find function calls with incorrect arguments, undeclared variables, etc.
ftnchek:
	/users/jeff/bin/ftnchek  -noextern -declare -library -usage=303 -array=2 -columns=132 -pretty=no-long-line -nonovice -quiet $(SRC)



##################### include dependency makefiles #####################

include $(DEPS)

