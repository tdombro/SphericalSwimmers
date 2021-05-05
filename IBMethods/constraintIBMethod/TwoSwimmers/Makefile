######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
#IBAMR_SRC_DIR   = /nas/longleaf/apps-dogwood/ibamr/2018-03/sfw/ibamr/IBAMR
#IBAMR_BUILD_DIR = /nas/longleaf/apps-dogwood/ibamr/2018-03/sfw/ibamr/ibamr-objs-opt
IBAMR_SRC_DIR   = $(HOME)/longleaf/sfw2/ibamr/IBAMR
IBAMR_BUILD_DIR = $(HOME)/longleaf/sfw2/ibamr/ibamr-objs-opt

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

## Needed for Xcode to capture compiler errors and warnings.
ifdef XCODE_VERSION_ACTUAL
CXXFLAGS += -fno-color-diagnostics
endif

######################################################################
## Build the application.
##
## NOTE: The following assumes that all .cpp files in the present
##       directory are used to build the executable.
SRC = $(wildcard *.cpp)  
CPPFLAGS += -MD -MP
PDIM = 2
OBJS = $(SRC:%.cpp=%.o) $(IBAMR_LIB_2D) $(IBTK_LIB_2D)
main2d: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o $@

clean:
	$(RM) *.o *.lo *.objs *.ii *.int.c *.d
	$(RM) -r .libs
	#$(RM) main2d

-include $(SRC:%.cpp=%.d)

