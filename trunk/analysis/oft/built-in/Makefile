TARGET = forward

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)
CFLAGS = -I./ -I./include/

GXX = g++

all: $(TARGET)

$(TARGET) : $(TARGET).c
	$(GXX) -o $@ $(TARGET).c $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	@rm -f $(TARGET) *.o
