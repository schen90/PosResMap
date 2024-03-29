# PosResMap

###folder
py:      python3 code
macros:  C++ code
npy:     npy files for pulse shape db and position
pulsedb: rootfiles for pulse shape db and position, converted from npy files
Map:     output position resolution Map

###steps
1. compile

   > cd macros
   > make
   > cd ..

2. convert npy to rootfile

   > py py/npy2root.py [dbname]
   > cd pulsedb/macros
   > make
   > cd ..
   > ./macros/MakeSignalBasis [dbname]  <==  add Xtalk and PreAmp response
   

3. calculate Map points position from a large Geant4 simulation data

   > ./macros/CalcMapPointsGrid  <==  generate map points at grid position


4. get simulated pulse shape in one crystal (uniform distribution, linear combination of pulse shape db from step 3)

   > Agata -Path ./trunk/ -b G4mac/sim1gamma.mac -run 0
   
   > ./macros/MakeDataG4 [inputfile] [outputfile]


5. analysis pulse shape data from step 4 to get analyzed position

   > ./macros/GridSearch/build/AnaDataGridSearch [inputfile] [outputfile] [detid] [dbfile]


6. update Map for type x detector using (analyzed position - simulated position) from step 5

   #./macros/MakeMapGrid [dettype] [PSAfile] [Mapfile]
   > ./macros/MakeMapGrid 0 rootfiles/typeA/AnaData_run Map/MapGrid.dat  <==  map at grid points
