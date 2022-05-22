#TARGET1 = CalcMapPoints
TARGET1 = CalcMapPointsGrid
#TARGET2 = MakeData
TARGET2 = MakeDataG4
#TARGET3 = AnaData
#TARGET4 = MakeMap
TARGET4 = MakeMapGrid

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lMinuit
CFLAGS = -I./ 

GXX = g++

all: $(TARGET1) $(TARGET2) $(TARGET4)

$(TARGET1) : $(TARGET1).C
	$(GXX) -o $@ $(TARGET1).C $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

$(TARGET2) : $(TARGET2).C
	$(GXX) -o $@ $(TARGET2).C $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

$(TARGET4) : $(TARGET4).C
	$(GXX) -o $@ $(TARGET4).C $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	@rm -f $(TARGET1) $(TARGET2) $(TARGET4) *.o
