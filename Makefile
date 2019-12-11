#
# General Makefile for the OfflineUser package
#
#
# If the environment variable AUGEROFFLINEROOT is not set
# AND the executable 'auger-offline-config' is not in your PATH
# the definition of the following variable is required
#
#AUGEROFFLINEROOT = /data/userdata9/holt-e/offline_141017/install_RdObserver
#
# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
USER_SRCS = $(wildcard *.cc)
#
# Give your executable a name
#
EXE = $(patsubst %.cc,%, $(wildcard *.cc))
#
#############################################################

## You should not need to change anything below this line ###

# Authors: T. Paul, S. Argiro, L. Nellen, D. Veberic
# $Id: Makefile.in 18270 2011-01-10 15:29:32Z huege $
# Send bug reports to http://www.auger.unam.mx/bugzilla/

.PHONY: all depend clean

ifdef AUGEROFFLINEROOT
  AUGEROFFLINECONFIG = $(AUGEROFFLINEROOT)/bin/auger-offline-config
else
  AUGEROFFLINECONFIG = auger-offline-config
endif

OBJS = $(USER_SRCS:.cc=.o)

CPPFLAGS    = $(shell $(AUGEROFFLINECONFIG) --cppflags)
CXXFLAGS    = $(shell $(AUGEROFFLINECONFIG) --cxxflags)
LDFLAGS     = $(shell $(AUGEROFFLINECONFIG) --ldflags)
MAIN        = $(shell $(AUGEROFFLINECONFIG) --main)
CONFIGFILES = $(shell $(AUGEROFFLINECONFIG) --config)
XMLSCHEMALOCATION = $(shell $(AUGEROFFLINECONFIG) --schema-location)

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^  $(CXXFLAGS) $(LDFLAGS) -lMinuit

%: %.in
	@echo -n "Generating $@ file..." 
	@sed -e 's!@''CONFIGDIR@!$(CONFIGFILES)!g;s!@''SCHEMALOCATION@!$(XMLSCHEMALOCATION)!g' $< >$@
	@echo "done"

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(CPPFLAGS) -MM $^ > $@

clean:
	- rm -f *.o  *.so *.ps core $(USER_XMLS) $(USER_XSDS) Make-depend



-include Make-depend
