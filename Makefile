CC = g++
CFLAGS = -c -g -fPIC -O2 -Wall -Wextra -pedantic -Wno-long-long -Wshadow -Werror -ansi -std=c++11
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

LIBS = -L$(PANDORA_DIR)/lib -lPandoraSDK $(shell root-config --glibs --evelibs)
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

PROJECT_INCLUDE_DIR = $(PANDORA_DIR)/PandoraMonitoring/include/
PROJECT_LIBRARY = $(PANDORA_DIR)/lib/libPandoraMonitoring.so

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
INCLUDES += -I$(shell root-config --incdir)

DEFINES  = -DMONITORING=1
DEFINES += -DROOT_EVE=1

SOURCES = $(wildcard $(PANDORA_DIR)/PandoraMonitoring/src/*.cc)
OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

all: library

library: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(PROJECT_LIBRARY)

-include $(DEPENDS)

%.o:%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(PROJECT_LIBRARY)

install:
ifdef INCLUDE_TARGET
	rsync -r --exclude=.svn $(PROJECT_INCLUDE_DIR) ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PROJECT_LIBRARY) ${LIB_TARGET}
endif
