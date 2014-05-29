DEFINES  = -DMONITORING=1
DEFINES += -DROOT_EVE=1

INCLUDES  = -I$(PANDORA_DIR)/PandoraMonitoring/include
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
INCLUDES += -I$(shell root-config --incdir)

CC = g++
CFLAGS = -c -g -fPIC -O2 -Wall -Wextra -pedantic -Wno-long-long -Wshadow -Werror -ansi
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES = $(wildcard $(PANDORA_DIR)/PandoraMonitoring/src/*.cc)

OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

LIBS = -L$(PANDORA_DIR)/lib -lPandoraSDK $(shell root-config --glibs --evelibs)
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS  = $(shell root-config --auxcflags)
LDFLAGS += $(LIBS) -Wl,-rpath

LIBRARY = $(PANDORA_DIR)/lib/libPandoraMonitoring.so

all: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

-include $(DEPENDS)

%.o:%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

install:
ifdef INCLUDE_TARGET
	rsync -r --exclude=.svn $(PANDORA_DIR)/PandoraMonitoring/include/ ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PANDORA_DIR)/lib/libPandoraMonitoring.so ${LIB_TARGET}
endif

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(LIBRARY)
