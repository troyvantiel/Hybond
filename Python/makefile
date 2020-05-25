g++ -o Hybond Hybond.cpp -std=c++11 -O3 -I/usr/include/python3.6 -L/usr/include/python3.6
-lpython3.6m

#
# 'make depend' uses makedepend to automatically generate dependencies
#				(dependencies are added to end of Makefile)
# 'make'		build executable file 'mycc'
# 'make clean'	removes all .o and executable files
#

# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g

# define and directories containing header files other than /usr/include
# python needs to be included as it is not a standard header file
#
INCLUDES = -I/usr/include/python3.6

# define library paths in addition to /usr/lib
# if I wanted to include libraries not in /usr/lib I'd specify
# their path using -Lpath, something like:
LFLAGS = -L/usr/include/python3.6

# define any libraries to link into executable:
# If I want to link in libraries (libx.so or libx.a) I use the -llibname
# option, something like (this will link in libmylib.so and libm.so)
LIBS = -lpython3.6m

# define the C source files
SRCS = Hybond.cpp

# define the C object files
#
# This uses Suffix Replacement within a macro:
#	$(name:string1=string2)
#		for each word in 'name' replace 'string1' and 'string2'
# Below we are replacing the suffix .cpp of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

#define the executable file
MAIN = Hybond

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by 
# deleting dependencies appended to the file from 'make append'
#

.PHONY: depend clean

all:	$(MAIN)
		@echo Simple compiler named Hybond has been complied

$(MAIN): $(OBJS)
		$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .cpp's
# it uses automatic variables $<: the name of the prerequisite of 
# the rule(a .c file) and $@: the name of the target of the rule(a .o file)
# (see the gnu make manual section about automatic variables)
.cpp.o:
		$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	





	



