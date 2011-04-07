#Path to project directory
PROJECT_DIR = YOUR_PATH_HERE

#Path to project dependency
PANDORAPFANEW_DIR = YOUR_PATH_HERE

DEFINES = -DROOT_EVE=1

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)/include/
PROJECT_SOURCE_DIR  = $(PROJECT_DIR)/src/
PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib/

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
INCLUDES += -I$(PANDORAPFANEW_DIR)/Framework/include/

CC = g++
CFLAGS = -c -Wall -g -w -fPIC -O2
CFLAGS += $(INCLUDES)
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES = $(wildcard $(PROJECT_SOURCE_DIR)*.cc)

OBJECTS = $(SOURCES:.cc=.o)

LIBS = $(shell $(ROOTSYS)/bin/root-config --glibs) -lEve
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS  = $(shell $(ROOTSYS)/bin/root-config --auxcflags)
LDFLAGS += $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libPandoraMonitoring.so

all: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $(DEFINES) $< -o $@

clean:
	rm -f ../src/*.o
	rm -f $(LIBRARY)
