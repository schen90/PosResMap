#ifndef SIGMAMAP_H
#define SIGMAMAP_H

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <TVector3.h>
#include <TMath.h>
#include <TMatrixD.h>

// position dependent positon resolution map
bool kposmapgrid = kFALSE;

// rotation matrix
const int NCrystals = 180;
TMatrixD Rt[NCrystals]; // 3x3 rot matrix world frame -> detector frame
TMatrixD Rt2[NCrystals]; // 3x3 rot matrix detector frame -> world frame
TMatrixD Tr[NCrystals]; // 3x1 trans matrix detector frame -> world frame

// segmentation
const int NCryType = 3;
const int Nsectors = 6;
const int Nslices  = 6;


// map at (x,y,z) grid point
const int MaxSteps = 50;
double GridDist = 2; // 2mm grid
double XYZrange[NCryType][3][2]; // range of XYZ
double MapG[NCryType][MaxSteps][MaxSteps][MaxSteps][3]; // position res map (sigma_phi,sigma_r,sigma_z) at grid point
int    SegG[NCryType][MaxSteps][MaxSteps][MaxSteps]; // segment at grid point
bool   FlgG[NCryType][MaxSteps][MaxSteps][MaxSteps];

// rotation matrix
void LoadMatrix(string LoopUpTable){
  ifstream fin(LoopUpTable.c_str());
  if(!fin){cerr<<"Cannot find CrystalPositionLookUpTable"<<endl; return;}
  cout<<"\e[1;32m find CrystalPositionLookUpTable... \e[0m"<<endl;

  int ir;
  int dummy_i;
  for(int i=0; i<NCrystals; i++){
    ir = -1;
    fin >> ir >> dummy_i;
    if(ir<0 || ir>=NCrystals) break;

    Tr[ir].ResizeTo(3,1);  Tr[ir].Zero();
    Rt[ir].ResizeTo(3,3);  Rt[ir].Zero();
    Rt2[ir].ResizeTo(3,3); Rt2[ir].Zero();
    for(int ix=0; ix<3; ix++) fin >> Tr[ir](ix,0);
    for(int ix=0; ix<3; ix++){
      fin >> dummy_i;
      for(int ix2=0; ix2<3; ix2++) fin >> Rt[ir](ix,ix2);
      for(int ix2=0; ix2<3; ix2++) Rt2[ir](ix,ix2) = Rt[ir](ix,ix2);
    }
    Rt[ir].Invert();    // change to rot from world frame -> detector frame
  }
  fin.close();
  return;
}


// read map at grid point
void LoadMapGrid(string mapfilename){
  const int kMaxBufLen = 500;
  ifstream fin(mapfilename.c_str());
  if(!fin){ cerr<<"cannot open map file: "<<mapfilename<<" !!!"<<endl; return;}

  cout<<"\e[1;32m read position resolution map at grid point \e[0m"<<endl;
  char *buffer = new char[kMaxBufLen];

  // init
  for(int itype=0; itype<NCryType; itype++){
    for(int iaxis=0; iaxis<3; iaxis++){
      XYZrange[itype][iaxis][0] = 1;
      XYZrange[itype][iaxis][1] = 0;
    }
    for(int ix=0; ix<MaxSteps; ix++)
      for(int iy=0; iy<MaxSteps; iy++)
        for(int iz=0; iz<MaxSteps; iz++){
          MapG[itype][ix][iy][iz][0] = -1;
          MapG[itype][ix][iy][iz][1] = -1;
          MapG[itype][ix][iy][iz][2] = -1;
	  SegG[itype][ix][iy][iz] = -1;
	  FlgG[itype][ix][iy][iz] = false;
        }
  }

  // read Map
  while(!fin.eof()){
    fin.getline(buffer,kMaxBufLen);

    if(strncmp(buffer,"#range",6)==0){
      cout<<"reading range for each crystal type"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype;
      while(1){
        fin >> itype;
        if(itype==-1) break;
        int tmp;
        for(int iaxis=0; iaxis<3; iaxis++){
          fin >> XYZrange[itype][iaxis][0] >> XYZrange[itype][iaxis][1] >> tmp;
	  XYZrange[itype][iaxis][0] -= 2*GridDist;
	  XYZrange[itype][iaxis][1] += 2*GridDist;
	  tmp += 4; // extend the range a bit
          if(tmp>MaxSteps) {cerr<<"change MaxSteps to >= "<<tmp<<endl; return;}
        }
      }
    }

    if(strncmp(buffer,"#Map",4)==0){
      cout<<"reading Map"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype, iseg, ipos[3];
      double pos[3],res[3];
      while(1){
        fin >> itype >> iseg >> pos[0] >> pos[1] >> pos[2];
        if(itype==-1) break;

        for(int iaxis=0; iaxis<3; iaxis++){
          if(pos[iaxis]-XYZrange[itype][iaxis][0]<-0.1 ||
             pos[iaxis]-XYZrange[itype][iaxis][1]>0.1){
            cout<<Form("axis%d: %.3f outside range %.3f ~ %.3f",iaxis,pos[iaxis],XYZrange[itype][iaxis][0],XYZrange[itype][iaxis][1])<<endl;
            ipos[iaxis] = -1;
          }else{
            ipos[iaxis] = (int)((pos[iaxis] - XYZrange[itype][iaxis][0]) / GridDist + 0.5);
            if(fabs(pos[iaxis]-(XYZrange[itype][iaxis][0]+GridDist*ipos[iaxis]))>0.1){
              cout<<Form("axis%d: %.3f not a grid point",iaxis,pos[iaxis])<<endl;
              ipos[iaxis] = -1;
            }
          }
        }

        for(int i=0; i<3; i++) fin>>res[i];
        if(ipos[0]<0 || ipos[1]<0 || ipos[2]<0) continue;
        for(int i=0; i<3; i++){
          MapG[itype][ipos[0]][ipos[1]][ipos[2]][i] = res[i];
        }
	SegG[itype][ipos[0]][ipos[1]][ipos[2]] = iseg;
	FlgG[itype][ipos[0]][ipos[1]][ipos[2]] = true;
      }
    }

  }

  fin.close();

  // extend the range a bit
  for(int itype=0; itype<NCryType; itype++)
    for(int ix=1; ix<MaxSteps-1; ix++)
      for(int iy=1; iy<MaxSteps-1; iy++)
	for(int iz=1; iz<MaxSteps-1; iz++)
	  if(SegG[itype][ix][iy][iz]<0){

	    for(int it=0; it<18; it++){
	      int tmpx = ix; int tmpy = iy; int tmpz = iz;
	      if     (it==0) tmpx--; else if(it==1) tmpx++;
	      else if(it==2) tmpy--; else if(it==3) tmpy++;
	      else if(it==4) tmpz--; else if(it==5) tmpz++;

	      else if(it==6){ tmpx--; tmpy--;}else if(it==7){ tmpx--; tmpy++;}
	      else if(it==8){ tmpx++; tmpy--;}else if(it==9){ tmpx++; tmpy++;}
	      else if(it==10){ tmpx--; tmpz--;}else if(it==11){ tmpx--; tmpz++;}
	      else if(it==12){ tmpx++; tmpz--;}else if(it==13){ tmpx++; tmpz++;}
	      else if(it==14){ tmpy--; tmpz--;}else if(it==15){ tmpy--; tmpz++;}
	      else if(it==16){ tmpy++; tmpz--;}else if(it==17){ tmpy++; tmpz++;}


	      if(FlgG[itype][tmpx][tmpy][tmpz]){
		SegG[itype][ix][iy][iz] = SegG[itype][tmpx][tmpy][tmpz];
		for(int iaxis=0; iaxis<3; iaxis++)
		  MapG[itype][ix][iy][iz][iaxis] = MapG[itype][tmpx][tmpy][tmpz][iaxis];
		break;
	      }
	    }

	  }


  kposmapgrid = kTRUE;
  return;
}


// position dependent position resolution [map at grid point]
int possigmagrid(int itype, double posxrel, double posyrel, double poszrel, double uncertainty3[]){
  if(!kposmapgrid) LoadMapGrid("include/Map/MapGrid.dat");

  uncertainty3[0] = 6.0; // average phi resolution 2 mm 
  uncertainty3[1] = 1.0; // average r resolution 1 mm
  uncertainty3[2] = 1.5; // average z resolution 1.5 mm

  double pos[3] = {posxrel,posyrel,poszrel};
  int ipos[3];
  // find nearest grid
  for(int ix=0; ix<3; ix++){
    ipos[ix] = (int)((pos[ix] - XYZrange[itype][ix][0]) / GridDist + 0.5);
  }

  // if cannot find map point
  for(int ix=0; ix<3; ix++) if(ipos[ix]<0 || ipos[ix]>MaxSteps-1) return -1;
  if(SegG[itype][ipos[0]][ipos[1]][ipos[2]]<0) return -1;

  // if find available map point
  for(int ix=0; ix<3; ix++) uncertainty3[ix] = MapG[itype][ipos[0]][ipos[1]][ipos[2]][ix];
  return SegG[itype][ipos[0]][ipos[1]][ipos[2]];
}


int FindNearestGrid(int itype, double &xrel, double &yrel, double &zrel){
  if(!kposmapgrid) return -1;
  
  double pos[3] = {xrel, yrel, zrel}; // mm
  int ipos[3];
  // find nearest grid
  for(int ix=0; ix<3; ix++){
    ipos[ix] = (int)((pos[ix] - XYZrange[itype][ix][0]) / GridDist + 0.5);
  }

  // if cannot find map point
  for(int ix=0; ix<3; ix++) if(ipos[ix]<0 || ipos[ix]>MaxSteps-1) return -1;
  if(SegG[itype][ipos[0]][ipos[1]][ipos[2]]<0) return -1;

  // if find available map point
  for(int ix=0; ix<3; ix++){
    pos[ix] = ipos[ix]*GridDist + XYZrange[itype][ix][0];
  }
  xrel = pos[0];  yrel = pos[1];  zrel = pos[2];
  return SegG[itype][ipos[0]][ipos[1]][ipos[2]];
}



void smearpointsmap(int number, double energy[],
		    double posx[], double posy[], double posz[],
		    double posxrel[], double posyrel[], double poszrel[],
		    int nseg[], int ndet[],
		    double errorfc[],double step, int nbsteps){

  int i,j;
  double tru,err,ener_u,ener_p;

  for (i=0;i<number;i++) {

    int itype = ndet[i]%3;
    
    ener_u = sqrt(1 + energy[i]*3.7)/2.35; // fwhm/2.35 in keV
    //ener_p = uncertainty*( (2.7+6.2*sqrt(0.1/energy[i]))/4.66 ); //from NIMA638
    double escale = (2.7+6.2*sqrt(0.1/energy[i]))/4.66; // energy dependence

    double uncertainty[3];
    int seg = possigmagrid(itype,posxrel[i]*10,posyrel[i]*10,poszrel[i]*10,uncertainty);
    if(seg<0) {
      //cout<<Form("cannot find resolution for det %d seg %d: %.2f %.2f %.2f",ndet[i],seg,posxrel[i]*10,posyrel[i]*10,poszrel[i]*10)<<endl;
      return;
    }
    
    TVector3 tmpposrel(posxrel[i]*10,posyrel[i]*10,0);
    double tmpx, tmpy;
    double tmpphi = tmpposrel.Phi();
    double tmpr   = tmpposrel.Mag();
    double tmpz   = poszrel[i]*10;
    uncertainty[0] = escale*uncertainty[0]/tmpr; // change to rad
    uncertainty[1] = escale*uncertainty[1]; // mm
    uncertainty[2] = escale*uncertainty[2]; // mm

    int tmpseg = -1;
    int nrepeat = 0;
    
    //while(tmpseg!=nseg[i]){
    while(tmpseg!=seg){

      nrepeat++;
      if(nrepeat>100){
	//cout<<Form("cannot smear point at det %d seg %d: %.2f %.2f %.2f",ndet[i],seg,posxrel[i]*10,posyrel[i]*10,poszrel[i]*10)<<endl;

	tmpx = posxrel[i]*10;
	tmpy = posyrel[i]*10;
	tmpz = poszrel[i]*10;
	tmpseg = FindNearestGrid(itype,tmpx,tmpy,tmpz);
	
	continue;
      }
      
      // phi
      tru = drand48(); if (tru < errorfc[0]) err = step;
      for(j=0;j<=nbsteps;j++){
	if(errorfc[j] > tru ){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48(); if(tru<=0.5) err=-err;      
      tmpphi = tmpposrel.Phi() + uncertainty[0]*err;

      // r
      tru = drand48(); if (tru < errorfc[0]) err = step;
      for(j=0;j<=nbsteps;j++){
	if(errorfc[j] > tru){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48(); if(tru<=0.5) err=-err;      
      tmpr = tmpposrel.Mag()+uncertainty[1]*err;

      // z
      tru=drand48(); if (tru < errorfc[0]) err = step;
      for(j=0;j<nbsteps;j++){
	if(errorfc[j] > tru){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48(); if(tru<=0.5) err=-err;      
      tmpz = poszrel[i]*10 + uncertainty[2]*err;

      TVector3 tmpposrel2;
      tmpposrel2.SetMagThetaPhi(tmpr,TMath::Pi()/2,tmpphi);
      tmpx = tmpposrel2.x();
      tmpy = tmpposrel2.y();
      
      tmpseg = FindNearestGrid(itype,tmpx,tmpy,tmpz);

    }
    posxrel[i] = tmpx/10;  posyrel[i] = tmpy/10;  poszrel[i] = tmpz/10;
    TMatrixD PosRel(3,1);
    PosRel(0,0)=tmpx;  PosRel(1,0)=tmpy;  PosRel(2,0)=tmpz;
    TMatrixD Pos(3,1);
    Pos = Rt2[ndet[i]] * PosRel + Tr[ndet[i]];
    posx[i] = Pos(0,0)/10;  posy[i] = Pos(1,0)/10;  posz[i] = Pos(2,0)/10;
    


    // energy
    tru=drand48(); if (tru < errorfc[0]) err = step;
    for(j=0;j<=nbsteps;j++){
      if(errorfc[j] > tru ){
	err = j*step;
	break;
      }
    }      
    tru=drand48(); if(tru<=0.5) err=-err;

    energy[i] = energy[i]+(err*ener_u)/1000; // add error in MeV
    if(energy[i]<0)energy[i] = fabs((err*ener_u)/1000);

  }

}


#endif
