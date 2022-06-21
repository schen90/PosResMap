#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include "/usr/include/x86_64-linux-gnu/sys/types.h"
#include "/usr/include/x86_64-linux-gnu/sys/stat.h"
//#include "/usr/include/sys/types.h"
//#include "/usr/include/sys/stat.h"
#include <unistd.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <iostream>

using namespace std;

// compile with g++ forward.c -I$ROOTSYS/include `root-config --libs` `root-config --glibs` -bind_at_load -lm -o forward

#define POSRESMAP
#ifdef POSRESMAP
#include "include/SigmaMap.h"
#endif


#define MAXYSPN 600
#define YSPNT 600
#define XSPN 8192

#define PI 3.1415927
#define rho_ge 5.32 /* g/cm3 */
#define Z_ge 32 
#define A_ge 74
#define N_av 6.022e23
#define r_0 2.8179e-15
#define mec2 0.511 /* electron mass MeV */
#define SQ(x) ((x)*(x))
#define CB(x) ((x)*(x)*(x))
#define alpha 0.0072993 /* fine structure constant 1/137 */




#define betafactor sqrt(1-SQ(voverc))

int spn_data[XSPN][MAXYSPN];


int main()
{
  
  char fn[250],buffer[250],filename[250],rwfile[250],rootfile[250];
  FILE *fp,*fp0,*fp1;
  int counter,particle_code;
  int i,j,n,ir,iff,nn;
	int k,kmax;
  int l,m,o,p,q,qq;
    int itot,nbrtotgam;
	double core[180],coren[180];
	  int nintu,nfiringdet;
  double u2,ddu,ee1u,facu,erfcu[1000];
  double etot,probtot[800],sigma_thet,power;
  double angn[800],alfamax,minprobtrack,minprobsing,yfl;
  int sn[800][800],numn[800],order[800];
  double angth,angph,etotalen,etotale,escatter,escattern,deltaesc,deltaescn,etest=0,total;
  int flagu[800],flagpair[800];
  int itotcurr,totnumbofgammas,nseg[800],ndet[800];
	long int det_count[180];
	double rmin,r[800][800],r_ge[800][800],costheta,gammacore[800];
  double costest,angtheta[800],angphi[800],alfa;
  double Ax,Ay,Az,Bx,By,Bz,vectproduct,r_center,r_vacuum,ri_orig,r_vacuumi;
    double cosangle=0,CI;
  double mincosdif,costest1,probtest,probtest1;
    double probtest3,probtest2,probtest4,photodistmax,test_sigma_thet,probtestpair;
  double e[800],x[800],y[800],z[800],et[800],xrel[800],en[800],yrel[800],zrel[800],probcomp;
  double cross1,cross2,gamma[100],gamma1[100];
    double Px[800],Py[800],Pz[800],firstx[800],firsty[800],firstz[800],dist[800];
    double totenergy,totnenergy;
    double secondx[800],secondy[800],secondz[800],dist2[800];
  double d_res,ri,rj,lambda2;
  double cross, lambda;
  int nb_int,numinter[800],mult,readflag,nshoot,iseg[800],interactionpair[800][20],nhit[180];
    int nintprev,nprev,nopair,flagseg,flagerr,flagpack;
  int number_of_unflagged;
  double sigma_pos,energythresh,lambda1;
  int nodirect,temp1[800],interaction[800][800],nmult[800];
  double coef1,coef2,ercos,eres,voverc,vrecoilx,vrecoily,vrecoilz;
  double deltaalfa,radius,radius_out,diffenergy;
  /* inner and outer radius now read from file */
    int  dummy_i,equal,flagdet[180];;
  double  xsource, ysource, zsource, dx, dy, dz, dd, dd_ge,cosTheta;
   
 char first_char='$';
int zero = 0L,l1,l2;
 int ntot=0;
	double esum,ecalore;
double rotzx[180],rotzy[180],rotzz[180],TrZ[180],rotxx[180],rotxy[180],rotxz[180],TrX[180],rotyx[180],rotyy[180],rotyz[180],TrY[180];
	double rxx[180],rxy[180],rxz[180],ryx[180],ryy[180],ryz[180],rzx[180],rzy[180],rzz[180],radiusrel;
	double detR[180];
    int numevt,corsphere;

    
    
  double sig_compt(double E);
  /* total compton cross section, cm2/atom , E en MeV */
   double sig_abs(double E); 
  /* total photoelectric cross section , cm2/atom, E en MeV */
  double range_process(double sig); 
  /* interaction length in cm = 1/SIGMA_macro */
  double proba(double range, double distance);
  /* probability of interaction at a given distance given a certain range */
  void swap(double v[],int m, int l);
  /* pair production cross section */
 void swapi(int v[],int m, int l);
  /* pair production cross section */
  double sig_pair(double E); 
 /* error in cosine */
  double err_cos(double x, double xp, double xpp, double y, double yp, double ypp, double z, double zp, double zpp, double r, double rp, double err);
    int packpoints(int number, double posx[], double posy[], double posz[], double posxrel[], double posyrel[], double poszrel[], double transx[], double transy[], double transz[],double energy[], int intnb[], int intorder[],int nseg[], int ndet[], double resolution);

  void smearpoints(int number, double energy[], double posx[], double posy[], double posz[],double errorfc[],double step, int nbsteps, double uncertainty); 
	void  wsec( char file_l[], int matrix_l[][MAXYSPN], int yspn_l);


  /*********************************************************************/
  /**** energies in MeV and positions in cm ****************************/
  /*********************************************************************/ 
 
    sprintf(fn,"OFTinput");
    fp1=fopen(fn,"r");
    if(!fp1){
        printf("cannot find OFTinput \n");
       
    }
    else{
    
    fgets(buffer,150,fp1);
    /* required emitted event multiplicity - or concatenation required in case of exp. data */
    /* number of emitted gammas to treat or detected events to treat in case of exp. data */
    sscanf(buffer,"%d %d ",&mult,&totnumbofgammas);
    memset(buffer,zero,sizeof(buffer));
    // psa hits file to be read
    fgets(buffer,150,fp1);
    sscanf(buffer,"%s ",filename);
    memset(buffer,zero,sizeof(buffer));
    // output root file
    fgets(buffer,150,fp1);
    sscanf(buffer,"%s ",rootfile);
    memset(buffer,zero,sizeof(buffer));
    // output radware file
    fgets(buffer,150,fp1);
    sscanf(buffer,"%s ",rwfile);
    /* position uncertainty in mm for costheta */
    /* threshold for track acceptance */
    /* threshold for single int acceptance */
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp1);
    sscanf(buffer,"%lf %lf %lf ",&sigma_thet,&minprobtrack,&minprobsing);
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp1);
    sscanf(buffer,"%d %d %d ",&nodirect, &nopair, &corsphere);/* if set to 1 do not check sing. interaction or pair production events */
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp1);
    sscanf(buffer,"%lf %d %d %d ",&energythresh, &flagseg, &flagerr, &flagpack);/* energy threshold in MeV, if set to 1 pack all interactions in a segment, smear x,y,z, with gaussian of sigma_pos at 100 keV & smear energies sqrt(1 + e*3.7)/2.35 , pack points within d_res of each other */
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp1);
    /* source position */
    sscanf(buffer,"%lf %lf %lf ",&xsource, &ysource, &zsource);
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp1);
    /* recoil velocity/c and direction */
    sscanf(buffer,"%lf %lf %lf %lf",&voverc,&vrecoilx,&vrecoily,&vrecoilz);

    fclose(fp1);
 
    
  /********************************************************************/
  /**************** parameter initialisation **************************/
  /********************************************************************/
    TFile *file = new TFile(rootfile,"recreate");
    
    TH1F *spectrumin_raw = new TH1F("core_energy_raw","energy",20000,0,20000);
    TH1F *spectrumin = new TH1F("core_energy","energy",20000,0,20000);
    TH1F *spectrumin_segtocore = new TH1F("core_energy_segtocore","energy",20000,0,20000);
    TH2F *core_vs_id = new TH2F("core_energy_vs_core_id","energy",20000,0,20000,180,0,180);
    TH2F *gammagammacore = new TH2F("core_gammagamma","energy",1000,0,2000,1000,0,2000);

    TH1F *spectrumin_calorraw = new TH1F("calor_energy_raw","energy",20000,0,20000);
    TH1F *spectrumin_calor = new TH1F("calor_energy","energy",20000,0,20000);
    
    TH1F *spectrumint = new TH1F("interaction_point_energy","energy",20000,0,20000);
    
    
    TH1F *spectrumout = new TH1F("TrackedDop","energy",20000,0,20000);
    TH1F *spectrumout_nodop = new TH1F("TrackednoDop","energy",20000,0,20000);
    TH1F *spectrumout_sing = new TH1F("TrackedDop_sing","energy",20000,0,20000);
    TH1F *spectrumout_sing_nodop = new TH1F("TrackednoDop_sing","energy",20000,0,20000);
    
    TH2F *coinc = new TH2F("tracked_gammagamma","energy",1000,0,2000,1000,0,2000);
    
   
    TH2F *specvsprob = new TH2F("cluster_energy_vs_prob","energy",1000,0,1,10000,0,20000);
    TH2F *specvsprobsing = new TH2F("sing_int_cluster_energy_vs_prob","energy",1000,0,1,10000,0,20000);
    TH2F *specvsrange = new TH2F("sing_cluster_energy_vs_range","energy",1000,0,9,10000,0,20000);

    TH2F *specvsprobpair = new TH2F("pair_prod_cluster_energy_vs_prob","energy",1000,0,1,10000,0,20000);
   
    
    TH2F *XY = new TH2F("XYrel","pos",100,-49.5,50.5,100,-49.5,50.5);
    TH2F *XZ = new TH2F("XZrel","pos",100,-49.5,50.5,100,-0.5,99.5);
    TH2F *YZ = new TH2F("YZrel","pos",100,-49.5,50.5,100,-0.5,99.5);

    TH2F *XYs = new TH2F("XYrels","smear pos",100,-49.5,50.5,100,-49.5,50.5);
    TH2F *XZs = new TH2F("XZrels","smear pos",100,-49.5,50.5,100,-0.5,99.5);
    TH2F *YZs = new TH2F("YZrels","smear pos",100,-49.5,50.5,100,-0.5,99.5);
    
    
    TH2F *XYabs = new TH2F("XYabs","pos",500,-500,500,500,-500,500);
    TH2F *XZabs = new TH2F("XZabs","pos",500,-500,500,500,-500,500);
    TH2F *YZabs = new TH2F("YZabs","pos",500,-500,500,500,-500,500);
    

    
    TH1F *stat_int= new TH1F("Nbdetectedhits","fold",100,-0.5,99.5);
    TH1F *stat_det = new TH1F("Nbfiringdets","fold",180,-0.5,179.5);
    TH1F *stat_det_id = new TH1F("Ndets_dist","detectornbdistribution",180,0,180);

    
 


   nprev = 0;
   nintprev = 0;
   readflag = 0;
 
 
    for (j=0;j<100;j++){
    
			gamma[j]=0;
			gamma1[j]=0;
		
     }
	
  
  k=0;
  n=0;

  itot = 0; /* number of reconstructed gamma-rays */

  d_res=0.5; /* min distance between points below which they are smeared*/

    
  nbrtotgam = 0; /* number of incident gamma-rays */
  nshoot = 0; /* number of events with given multiplicity */
  
  nintu=1000; /* initialisation for random number vs uncertainty procedure */
  ddu=0.005; /* step in random number vs uncertainty procedure */
  ee1u=0; /* initial value for random number vs uncertainty procedure */
   
   photodistmax = 4. ; /*min distance in cm between 2 points to process single int. point */
     
   sigma_pos=0.2 ; /* width in mm of gaussian position smearing at 100 keV  sigma = sigma_pos * sqrt(0.1/e) */
	
   eres = 3e-3; /* energy resolution */
  
    alfa = 0.05; // min angular separation in phi and theta in rad
    deltaalfa = 0.05; // step in alpha
  
    kmax = 10; /* max number of interactions for a given incident gamma-ray */


  /**********************************************************/
  /* mapping random number vs uncertainty in units of sigma */
  /**********************************************************/
  
   
  u2=ddu;
  
  for(j=1;j<=nintu;j++) {   
    ee1u=ee1u+exp(-(u2*u2/2.));   
    erfcu[j]=ee1u;
    u2=u2+ddu;
  }
  
  facu=ee1u;
  for(j=1;j<=nintu;j++) {
    erfcu[j]=erfcu[j]/facu; /* random => ddu * j  in units of sigma_pos */  
  }
  
  /*********************************************************************/
  /* reading in source position and interaction positions and energies */
  /* and opening files for spectra- performances - detailed output *****/
  /*********************************************************************/
  
	

   	
	
	sprintf(fn,"CrystalPositionLookUpTable");
	fp0=fopen(fn,"r");
    if(!fp0){
        printf("Cannot find CrystalPositionLookUpTable \n");
      
    }
    else{
	j=0;
	for(i=0;i<180;i++){
        flagdet[i]=0;
		fgets(buffer,150,fp0);
			//	  printf("%s \n",buffer);
		sscanf(buffer,"%d %d %lf %lf %lf ",&ir,&dummy_i,&TrX[j],&TrY[j],&TrZ[j]);
		memset(buffer,zero,sizeof(buffer));
		fgets(buffer,150,fp0);
	//	printf("%s \n",buffer);
		sscanf(buffer,"%d %lf %lf %lf  ",&dummy_i,&rotxx[j],&rotxy[j],&rotxz[j]);
		//		  printf("ir = %d \n",ir);
		memset(buffer,zero,sizeof(buffer));
		fgets(buffer,150,fp0);
				 //   printf("%s \n",buffer);
		sscanf(buffer,"%d %lf %lf %lf  ",&dummy_i,&rotyx[j],&rotyy[j],&rotyz[j]);
		memset(buffer,zero,sizeof(buffer));
		fgets(buffer,150,fp0);
			 //   printf("%s \n",buffer);
		sscanf(buffer,"%d %lf %lf %lf  ",&dummy_i,&rotzx[j],&rotzy[j],&rotzz[j]);
	
					
		
				
		
		detR[j]=rotxx[j]*(rotyy[j]*rotzz[j]-rotzy[j]*rotyz[j]) - rotxy[j]*(rotyx[j]*rotzz[j]-rotzx[j]*rotyz[j])+rotxz[j]*(rotyx[j]*rotzy[j]-rotzx[j]*rotyy[j]);
		rxx[j]=(rotyy[j]*rotzz[j]-rotzy[j]*rotyz[j])/detR[j];
		rxy[j]=-(rotxy[j]*rotzz[j]-rotxz[j]*rotzy[j])/detR[j];
		rxz[j]=(rotxy[j]*rotyz[j]-rotyy[j]*rotxz[j])/detR[j];
		ryx[j]=-(rotyx[j]*rotzz[j]-rotzx[j]*rotyz[j])/detR[j];
		ryy[j]=(rotxx[j]*rotzz[j]-rotxz[j]*rotzx[j])/detR[j];
		ryz[j]=-(rotxx[j]*rotyz[j]-rotyx[j]*rotxz[j])/detR[j];
		rzx[j]=(rotyx[j]*rotzy[j]-rotzx[j]*rotyy[j])/detR[j];
		rzy[j]=-(rotxx[j]*rotzy[j]-rotzx[j]*rotxy[j])/detR[j];
		rzz[j]=(rotxx[j]*rotyy[j]-rotxy[j]*rotyx[j])/detR[j];
		
	 		
		
		j++;
	}
	
	fclose(fp0);

#ifdef POSRESMAP
	LoadMatrix("CrystalPositionLookUpTable");
#endif

  fp = fopen(filename,"r");
  if(!fp)
    {
      printf("could not open file %s\n",filename);
      
    }
  else
    {
      printf("opened file %s\n",filename);
    
  
  
  counter=0;
  
  
  /******************************************************************/
  /*********** read buffer at the beginning of file *****************/
  /******************************************************************/
  
  while(1) {
    memset(buffer,zero,sizeof(buffer));
    fgets(buffer,150,fp);
 
  
    if (strncmp(buffer,"SUMMARY",7)==0) {
      
      sscanf(buffer,"%s %lf %lf %d %d %d %d %d %d %d %d ",fn,&radius,&radius_out,&dummy_i,&dummy_i,&dummy_i,&dummy_i,&dummy_i,&dummy_i,&dummy_i,&dummy_i);
       radius=radius/10; // mm to cm
      radius_out=radius_out/10;
      printf("radius inner = %f outer = %f \n",radius,radius_out);
    }

  

    if (strncmp(buffer,&first_char,1)==0) break;
  }
  

 readfile: ; // beginning of an event


  

  
  /************************************************************************/
  /**************** read interaction points and energies  ****************/
  /***********************************************************************/
  
 
  
	for(i=0;i<180;i++){
		core[i]=0; // raw core energy
        coren[i]=0; // core energy with threshold
        nhit[i]=0 ;// number of energy depositions per detector with threshold
		
	}
	
  
 
  ir=0;
  nb_int=0;
  alfa =0.05;
  i=0;

  
  while(1) { 
    
    
  readread: // next interaction point 
    

    if(nintprev > totnumbofgammas) { // maximum nb of gamma-rays treated
      readflag =1;
 
		break;
    }
      if(fgets(buffer,sizeof(buffer),fp)== NULL){
		  printf("end of file  \n");
		  readflag=1;
		  break;
	  }
	  
	  sscanf(buffer,"%d",&particle_code);
	  
	
	  if(particle_code==-101)goto readread;
	  
    
    if(counter > mult)break;
    
    sscanf(buffer,"%d %lf %lf %lf %lf %d",&iff,&e[ir],&x[ir],&y[ir],&z[ir],&nprev);
  
	  
	
     
    
    if(readflag ==1) break;
    
    if(iff == -1) {  
        
       
        
        numevt=nprev;

		total = e[ir];
		etest=0;
      nintprev++;
		if(counter==mult)break;
      counter++;
      i=0;
      goto readread;      
    }
   
	
    numinter[ir]=nintprev;
    order[ir] = i;
    i++;
    nseg[ir]=nprev; // !!!!  0->35 with 0->5 down the detector
    ndet[ir]=iff;
      flagdet[ndet[ir]]=1;
	  det_count[ndet[ir]]++;
	   etest = etest+e[ir];
    e[ir]=1e-3*e[ir];
  
	// coordinates in individual detector refrence frame
      
	  xrel[ir]=rxx[ndet[ir]]*(x[ir]-TrX[ndet[ir]])+rxy[ndet[ir]]*(y[ir]-TrY[ndet[ir]])+rxz[ndet[ir]]*(z[ir]-TrZ[ndet[ir]]);
	  yrel[ir]=ryx[ndet[ir]]*(x[ir]-TrX[ndet[ir]])+ryy[ndet[ir]]*(y[ir]-TrY[ndet[ir]])+ryz[ndet[ir]]*(z[ir]-TrZ[ndet[ir]]);

	  zrel[ir]=rzx[ndet[ir]]*(x[ir]-TrX[ndet[ir]])+rzy[ndet[ir]]*(y[ir]-TrY[ndet[ir]])+rzz[ndet[ir]]*(z[ir]-TrZ[ndet[ir]]);

      // position of front face of detector for each interaction point in abs. coordinates
      
      Px[ir] = TrX[ndet[ir]]/10;
      Py[ir] = TrY[ndet[ir]]/10;
      Pz[ir] = TrZ[ndet[ir]]/10;
      
      radiusrel=sqrt(Px[ir]*Px[ir]+Py[ir]*Py[ir]+Pz[ir]*Pz[ir]);
   

	  
	  XZ->Fill(xrel[ir],zrel[ir]);
	  XY->Fill(xrel[ir],yrel[ir]);
	  YZ->Fill(yrel[ir],zrel[ir]);
	 
	  XYabs->Fill(x[ir],y[ir]);
	  XZabs->Fill(x[ir],z[ir]);
	  YZabs->Fill(y[ir],z[ir]);
	  
	  x[ir]=x[ir]/10.;
	  y[ir]=y[ir]/10.;
	  z[ir]=z[ir]/10.;
      
      xrel[ir]=xrel[ir]/10.;
      yrel[ir]=yrel[ir]/10.;
      zrel[ir]=zrel[ir]/10.;

	  
		

    nb_int++;
    ir++;
      
      if(ir>800){
          printf("too many interaction points in this event....there must be a problem \n");
          readflag=1;
          break;
      }
    
  }
    
    
    
	
	 	
  if(nb_int != 0)nshoot++;
 
  nbrtotgam = nbrtotgam+mult;
  itotcurr=itot;
  
  if((nintprev)%1000 ==0){
    //printf("%d events \n",nintprev);
    cout<<"\r"<<nintprev<<" events "<<flush;
   }
  
 	
	
	/***********************************************************************************/
	/**************** manipulate data and construct input spectra ***********************/
	/**********************************************************************************/
	
 
  	
  
  if(flagpack==1)nb_int=packpoints(nb_int, x, y, z,xrel,yrel,zrel,Px,Py,Pz,e,numinter,order,nseg,ndet,d_res);

  if(flagerr==1){
#ifdef POSRESMAP
    smearpointsmap(nb_int, e, x, y, z, xrel, yrel, zrel, nseg, ndet, erfcu, ddu, nintu);
#else
    smearpoints(nb_int, e, x, y, z, erfcu, ddu, nintu,sigma_pos);
#endif
  }

  for(i=0; i<nb_int; i++){
    XZs->Fill(xrel[i]*10,zrel[i]*10);
    XYs->Fill(xrel[i]*10,yrel[i]*10);
    YZs->Fill(yrel[i]*10,zrel[i]*10);
  }
  
   if(flagseg==1){
        
        for(i=0;i<nb_int;i++) {
            
            for(j=i+1;j<nb_int;j++){
                if(e[i]!=0 && ndet[i]==ndet[j] && nseg[i]==nseg[j]){
                    
                    en[i]=e[i]+e[j];
                    x[i]=((e[i]*x[i])+(e[j]*x[j]))/en[i];
                    y[i]=((e[i]*y[i])+(e[j]*y[j]))/en[i];
                    z[i]=((e[i]*z[i])+(e[j]*z[j]))/en[i];
                    
                    
                    e[j]=0;
                    
                    e[i]=en[i];
                    
                }
                
            }
            
        }
        
        
        j=0;
        for(i=0;i<nb_int;i++){
            if(e[i]!=0){
                e[j]=e[i];
                x[j]=x[i];
                y[j]=y[i];
                z[j]=z[i];
                numinter[j]=numinter[i];
                nseg[j]=nseg[i];
                ndet[j]=ndet[i];
                order[j]=order[i];
                xrel[j]=xrel[i];
                yrel[j]=yrel[i];
                zrel[j]=zrel[i];
                Px[j]=Px[i];
                Py[j]=Px[i];
                Pz[j]=Pz[i];
                j++;
            }
            
        }
        nb_int = j;
    }
	
    
    /********** setting thresholds on detected hit energy ***************/
    
    
   
    totenergy=0;
    totnenergy=0;
    
    for(i=0;i<nb_int;i++){
       core[ndet[i]]=core[ndet[i]]+e[i];
        totenergy=totenergy+e[i];

        if(e[i]<energythresh)e[i]=0;
            
        
        
    }
  

    if(totenergy>0 && totenergy*1e3<20000)spectrumin_calorraw->Fill(totenergy*1e3);
   
    
	if(nb_int<100)stat_int->Fill(nb_int);
  


	
	for(i=0;i<nb_int;i++){
        iseg[i]=1;
        }
	
    l=0;
	
	nfiringdet=0;
    ecalore=0;
    
	for(i=0;i<nb_int;i++){
                ecalore=ecalore+e[i];
		if(e[i]<20000)spectrumint->Fill(e[i]*1e3);
		
		if(iseg[i]==1){
            if(core[ndet[i]]>0 && core[ndet[i]]*1e3<20000)spectrumin_raw->Fill(core[ndet[i]]*1e3);

			nfiringdet++;
			if(ndet[i]<180)stat_det_id->Fill(ndet[i]);
			esum=e[i];
			for(j=i+1;j<nb_int;j++){
				if(ndet[j]==ndet[i]){
					esum=esum+e[j];
					iseg[j]=0;
				}
			}
			
            if(esum<20000){
                
               
            }

			
		
			coren[ndet[i]]=esum;
            if(coren[ndet[i]]>0 && coren[ndet[i]]*1e3<20000){
                spectrumin->Fill(coren[ndet[i]]*1e3);
          
                                                                          
                                                            }
            gammacore[l]=esum*1e3;
          
          
            l++;
		}
		
		
	}
    
    if(ecalore>0 && ecalore*1e3<20000){
        
        yfl=drand48();
        j=(int) (ecalore*1e3 + yfl);
        
      
            spectrumin_calor->Fill(ecalore*1e3);
            if(j>0 && j<8192)spn_data[j][3]++;
            
        

    }
    
    for(i=0;i<l;i++){
     
            for(j=i+1;j<l;j++){
                if(gammacore[i]<2000 && gammacore[j]<2000){
                    
                    
                    gammagammacore->Fill(gammacore[i],gammacore[j]);
                    gammagammacore->Fill(gammacore[j],gammacore[i]);
                    
                   

                    
                    
                }
                
            }
        
    }
	
	if(nfiringdet>0 && nfiringdet<=180)stat_det->Fill(nfiringdet);
	
	
	
	  
 /**** redistributing the missing energy on enery depositions *****/
    
    for(i=0;i<nb_int;i++){
        iseg[i]=1;
        if(e[i]>0)nhit[ndet[i]]++;
    }

    
        
    for (i=0;i<nb_int;i++){
        if(iseg[i]==1){
            
            diffenergy=core[ndet[i]]-coren[ndet[i]];
                if(e[i]>0)e[i]=e[i]+diffenergy/nhit[ndet[i]];
                esum=e[i];
                for(j=i+1;j<nb_int;j++){
                    if(ndet[j]==ndet[i]){
                        if(e[j]>0)e[j]=e[j]+diffenergy/nhit[ndet[i]];
                        esum=esum+e[j];
                        iseg[j]=0;
                    }
                }
            if(esum>0 && esum*1e3<20000){
                spectrumin_segtocore->Fill(esum*1e3);
                core_vs_id->Fill(esum*1e3,ndet[i]);
                yfl=drand48();
                k=(int) (1e3*esum+yfl);
                if(k<8192)spn_data[k][0]++;
   
            }
            
            }
        }
        
    

  /*****************************************************/
  /******** compute distance source-int and ************/
  /*************theta and phi for each int *************/
  /*****************************************************/
  

  
  for (i=0; i< nb_int; i++) {
  
 
    r[i][i] = sqrt(SQ(x[i]-xsource)+SQ(y[i]-ysource)+SQ(z[i]-zsource));
    angtheta[i] = acos((z[i]-zsource)/r[i][i]);
    angphi[i] = atan2((y[i]-ysource),(x[i]-xsource));
    if(angphi[i] <0)angphi[i]=2*PI+angphi[i];

  }
  
  /*****************************************************/
  /********* sort according to increasing theta*********/
  /*****************************************************/
  
  
  
  for (i=0; i< nb_int;i++)
    {
      for (j=i+1;j<nb_int; j++)
	{	  	  	  
	  if(angtheta[j] < angtheta[i])
	    {
	      swap(e, i, j);
	      swap(x, i, j);
	      swap(y, i, j);
	      swap(z, i, j);
			swap(xrel, i, j);
			swap(yrel, i, j);
			swap(zrel, i, j);
            swap(Px,i,j);
            swap(Py,i,j);
            swap(Pz,i,j);
			swapi(nseg,i,j);
			swapi(ndet,i,j);
	      swapi(numinter,i,j);
	      swapi(order,i,j);
	      swap(angtheta, i, j);
	      swap(angphi, i, j);
	    }	  
	}  
		
    }
  
  
 
  
  /*****************************************************************/
  /******** computing all distances between interactions   *********/
  /*****************************************************************/
  
  
  
  for (i=0;i<nb_int;i++) {

    r[i][i] = sqrt(SQ(x[i]-xsource)+SQ(y[i]-ysource)+SQ(z[i]-zsource));

    ri_orig=sqrt(SQ(x[i])+SQ(y[i])+SQ(z[i])); //distance from geometrical point to i
    
    Ax = -x[i]; 
    Ay = -y[i];
    Az = -z[i];
      
    
    
    Bx=xsource-x[i];
    By=ysource-y[i];
    Bz=zsource-z[i];
    
    vectproduct = sqrt(SQ((Ay*Bz)-(Az*By))+SQ((Az*Bx)-(Ax*Bz))+SQ((Ax*By)-(Ay*Bx)));
    
      
      
    /* distance from center to line joining points : rcenter = sintheta * distance center-i */
    
    r_center = vectproduct/r[i][i]; 
    
    r_ge[i][i]=r[i][i];
    
    if(r_center <= radius) {
      r_vacuum=sqrt(SQ(radius)-SQ(r_center)); // distance from perpendicular intersection to i-source
      // from origin to inner surface of shell
      
      r_vacuumi=sqrt(SQ(ri_orig)-SQ(r_center)); // distance from perpendicular intersction to i-source
      //from origin to i
      
      r_ge[i][i]=r_vacuumi-r_vacuum;
      
        
    }
      
      
      // correcting for proper Ge distances into Ge - this adds efficiency below 50 keV
      
      if(corsphere==1){ //
          
         
        cosangle = (-Bx*Px[i] - By*Py[i] - Bz*Pz[i])/(r[i][i] * radiusrel); // all geometries

          CI = zrel[i]/cosangle; // distance from front face to interaction point i
          
      
          r_ge[i][i]=CI;
       
          
      }
    
    if(r_ge[i][i] < 0)r_ge[i][i]=0.03;
	

      
    

    
    for (j=i+1;j<nb_int;j++) {
      
      r[i][j]= sqrt(SQ(x[i]-x[j])+SQ(y[i]-y[j])+SQ(z[i]-z[j]));
      r[j][i]=r[i][j];
      
     
      ri = ri_orig;
      rj = sqrt(SQ(x[j]) + SQ(y[j]) + SQ(z[j]));;
      if (ri > rj) {
		k = j;
		l = i;
      } else {
		k = i;
		l = j;
      }
      
      /********************************************************/
      /***** determine effective length in Ge *****************/
      /********************************************************/	  
      
      
      
      /* vector (interaction point i - geometrical center) */
	Ax = -x[k]; 
	Ay = -y[k];
	Az = -z[k];
	
	
	
	/* vector (interaction point i - interaction point j  */ 	  
	Bx= x[l]-x[k];
	By= y[l]-y[k];
	Bz= z[l]-z[k];
	
	costest = Ax * Bx + Ay * By + Az * Bz;
	costest1 = sqrt(SQ(Ax) + SQ(Ay) + SQ(Az));
	costest = costest / (costest1 * r[i][j]);
	
	
	/* norm of vector product = sintheta * distance center-i * distance i-j */
	vectproduct = sqrt(SQ((Ay*Bz)-(Az*By))+SQ((Az*Bx)-(Ax*Bz))+SQ((Ax*By)-(Ay*Bx)));
	
	/* distance from center to line joining points : rcenter = sintheta * distance center-i */	 
	r_center = vectproduct/r[i][j]; 
	
	r_ge[i][j] = r[i][j];
	
	
	if((r_center < radius)  && (costest > 0)) {
	  /* if < radius, the photon does not go through */
	  /* Ge all the way from j to i */
	  /* right triangle:  r_vacuum^2 + rcenter^2 = radius^2 */
	  /* 2*r_vacuum = distance to take away from rmin */
	  r_vacuum = sqrt(SQ(radius)-SQ(r_center));
	  r_ge[i][j]= r[i][j] - 2*r_vacuum;
	  
	  
	}
	
	if(r_ge[i][j] < 0)r_ge[i][j]=0.02;
	r_ge[j][i]=r_ge[i][j];
		


    }
  }        
  
  
  
  
  
  /****************************************************/
  /*************** find cluster ***********************/
  /****************************************************/
  
  n=0;

/* calculating maximum opening angle from the number of interaction points */ 
 
  power = pow(((nb_int+2)/3.),0.9);
  alfamax=acos(1-2/power);
	

 findcluster: ;
 
 
  
  for(i=0; i<nb_int;i++) flagu[i]=0;
	  for(i=0;i<nb_int;i++) {

    rmin=0;
    et[n]=0;
    k=0;
    if(flagu[i] ==0) {
      sn[n][k] = i;
      et[n] = et[n]+e[i];
      angth=angtheta[sn[n][k]];
      
      angph = angphi[sn[n][k]];
      k++;
      flagu[i]=1;
      
      for(j=0;j<nb_int;j++) {
	if(j!=i && flagu[j]==0) {

	  costest1 = (cos(angth)*cos(angtheta[j]));
	  costest1=costest1+((sin(angth)*sin(angtheta[j]))*cos(angphi[j]-angph));
	  costest1=acos(costest1);
	  
	  if((fabs(costest1) <= alfa)) {

	    et[n] = et[n]+e[j];
	    sn[n][k] = j;
	    k++;
	    flagu[j]=1;
	    
	    angth=angtheta[j];  
	    angph=angphi[j];
	   
	  }
	  
	  
	}
	if(k == kmax) break;
	
      }
		
		equal=0;
      numn[n]=k;
      angn[n]=alfa;
		
		// accept new cluster only if unique
		for(nn=0;nn<n;nn++){
			
			if(numn[nn]==numn[n] && fabs(et[nn]-et[n])<1e-6)
			{
				equal=1;
				break;	
			}
			   }
		
     if(equal==0)n++;
     
    }

    
  }
  ntot=ntot+n;
	


  /**************************************************/ 
  /*********** loop on alpha values *****************/
  /**************************************************/ 
  
  if(alfa < alfamax) {
    alfa =alfa + deltaalfa;
    goto findcluster;
  }
  
    
   
  
  /*******************************************************/
  /******** compute figure of merit of clusters **********/
  /*******************************************************/
  
  for(i=0;i<n;i++) {
    etotale = et[i];
    flagu[i]=0;
    flagpair[i]=0;
    
    if(numn[i]==1)goto suite; /* keep single int points for later */
    mincosdif=0;
    
      
      if(nopair==0){
          for(j=0;j<numn[i];j++){
              
              if((etotale-1.022)>=-sqrt(numn[i]*SQ(eres)) && fabs(e[sn[i][j]]-et[i]+1.022)<=sqrt((numn[i]+1)*SQ(eres))){
                  
                  flagpair[i]=1;
                  interactionpair[i][0]=sn[i][j];
                  cross= sig_pair(etotale);
                  cross1 = sig_compt(etotale);
                  cross2 = sig_abs(etotale);
                  coef1 = cross/(cross+cross1+cross2);
                  lambda1=range_process(cross);
                  probtestpair=coef1*proba(lambda1,r_ge[sn[i][j]][sn[i][j]]);
                  specvsprobpair->Fill(probtestpair,etotale*1e3);
              }
              
          }
      }
    for(j=0;j<numn[i];j++) {
           
   
      for(l=0;l<numn[i];l++) {

	if(l!=j) {
	 
	  
	  costheta = acos(z[sn[i][j]]/r[sn[i][j]][sn[i][j]]);
	  costheta=(x[sn[i][j]]-xsource)*(x[sn[i][l]]-x[sn[i][j]])+(y[sn[i][j]]-ysource)*(y[sn[i][l]]-y[sn[i][j]])+(z[sn[i][j]]-zsource)*(z[sn[i][l]]-z[sn[i][j]]);
	  costheta= costheta/(r[sn[i][j]][sn[i][j]]*r[sn[i][j]][sn[i][l]]);
		test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][j]])/2.35;
	  ercos= err_cos(xsource,x[sn[i][j]],x[sn[i][l]],ysource,y[sn[i][j]],y[sn[i][l]],zsource,z[sn[i][j]],z[sn[i][l]],r[sn[i][j]][sn[i][j]],r[sn[i][j]][sn[i][l]],sigma_thet);
	  
	  escattern = etotale/(1+(etotale/mec2)*(1-costheta));
	  escatter = etotale - e[sn[i][j]];
	  
	  deltaescn = SQ(escattern)*ercos/mec2;
	  deltaesc = sqrt((numn[i]+1)*SQ(eres));
	  
	 
	  probtest=exp(-2*SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc))); 
	
 
	     if(probtest < 0.005) goto nextsecond;
	  cross1 = sig_compt(etotale);
	  cross= sig_pair(etotale);
	  coef1 = cross1/(cross1+sig_abs(etotale)+cross);
	  cross2 = sig_compt(escatter);
	  cross=sig_pair(escatter);
	  coef2 = cross2/(cross2+cross+sig_abs(escatter));
	  if(numn[i] == 2) {
	    cross2=sig_abs(escatter);	   
	    coef2= cross2/(cross2+cross+sig_compt(escatter));
	  }
	  lambda1=range_process(cross1);
	  lambda2=range_process(cross2);
	  probtest=probtest*SQ(coef1*proba(lambda1,r_ge[sn[i][j]][sn[i][j]]))*coef2*proba(lambda2,r_ge[sn[i][j]][sn[i][l]]);
	  
	
	  
	  for(m=0;m<numn[i];m++) {
	    if(m!=j && m!=l) {
	      
	      costheta=(x[sn[i][l]]-x[sn[i][j]])*(x[sn[i][m]]-x[sn[i][l]])+(y[sn[i][l]]-y[sn[i][j]])*(y[sn[i][m]]-y[sn[i][l]])+(z[sn[i][l]]-z[sn[i][j]])*(z[sn[i][m]]-z[sn[i][l]]);
	      costheta= costheta/(r[sn[i][j]][sn[i][l]]*r[sn[i][l]][sn[i][m]]);
			test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][l]])/2.35;
	      ercos= err_cos(x[sn[i][j]],x[sn[i][l]],x[sn[i][m]],y[sn[i][j]],y[sn[i][l]],y[sn[i][m]],z[sn[i][j]],z[sn[i][l]],z[sn[i][m]],r[sn[i][j]][sn[i][l]],r[sn[i][l]][sn[i][m]],sigma_thet);
	      etotalen=etotale - e[sn[i][j]];
	      escattern = etotalen/(1+(etotalen/mec2)*(1-costheta));
	      escatter = etotalen - e[sn[i][l]];
	      
	      deltaescn = SQ(escattern)*ercos/mec2;
	      deltaesc = sqrt((numn[i]+2)*SQ(eres));
	      
	      probcomp=exp(-SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));
 
	          if(probcomp < 0.005) goto nextthird;
	      cross1 = sig_compt(etotalen);
	      cross=sig_pair(etotalen);
	      coef1 = cross1/(cross1+cross+sig_abs(etotalen));
	      cross2 = sig_compt(escatter);
	      cross=sig_pair(escatter);
	      coef2 = cross2/(cross2+cross+sig_abs(escatter));
	      if(numn[i] == 3) {
		cross2=sig_abs(escatter);
		coef2= cross2/(cross2+cross+sig_compt(escatter));
		
	      }
	      lambda1=range_process(cross1);
	      lambda2=range_process(cross2);
	      
	      probcomp = probtest*probcomp*coef2*proba(lambda2,r_ge[sn[i][l]][sn[i][m]]);
	      
	      
	      for(o=0;o<numn[i];o++) {
		if(o!=j && o!=l && o!=m) {
		  costheta=(x[sn[i][m]]-x[sn[i][l]])*(x[sn[i][o]]-x[sn[i][m]])+(y[sn[i][m]]-y[sn[i][l]])*(y[sn[i][o]]-y[sn[i][m]])+(z[sn[i][m]]-z[sn[i][l]])*(z[sn[i][o]]-z[sn[i][m]]);
		  costheta= costheta/(r[sn[i][l]][sn[i][m]]*r[sn[i][m]][sn[i][o]]);
			test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][m]])/2.35;
		  ercos= err_cos(x[sn[i][l]],x[sn[i][m]],x[sn[i][o]],y[sn[i][l]],y[sn[i][m]],y[sn[i][o]],z[sn[i][l]],z[sn[i][m]],z[sn[i][o]],r[sn[i][l]][sn[i][m]],r[sn[i][m]][sn[i][o]],sigma_thet);
		  etotalen=etotale - e[sn[i][j]] - e[sn[i][l]];
		  escattern = etotalen/(1+(etotalen/mec2)*(1-costheta));
		  escatter = etotalen - e[sn[i][m]];
		  
		  deltaescn = SQ(escattern)*ercos/mec2;
		  deltaesc = sqrt((numn[i]+3)*SQ(eres));
		  
		  probtest1=exp(-SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));
		
		  if(probtest1 < 0.005) goto nextfourth;
		  cross1 = sig_compt(etotalen);
		  cross=sig_pair(etotalen);
		  coef1 = cross1/(cross1+cross+sig_abs(etotalen));
		  cross2 = sig_compt(escatter);
		  cross=sig_pair(escatter);
		  coef2 = cross2/(cross2+cross+sig_abs(escatter));
		  if(numn[i] == 4) {
		    cross2=sig_abs(escatter);
		    coef2= cross2/(cross2+cross+sig_compt(escatter));
		    
		  }
		  lambda1=range_process(cross1);
		  lambda2=range_process(cross2);
		  probtest1=probcomp*probtest1*coef2*proba(lambda2,r_ge[sn[i][m]][sn[i][o]]);
		  
		  for(p=0;p<numn[i];p++) {
		    
		    if(p!=o && p!=m && p!=l && p!=j) {
		      
		      costheta=(x[sn[i][o]]-x[sn[i][m]])*(x[sn[i][p]]-x[sn[i][o]])+(y[sn[i][o]]-y[sn[i][m]])*(y[sn[i][p]]-y[sn[i][o]])+(z[sn[i][o]]-z[sn[i][m]])*(z[sn[i][p]]-z[sn[i][o]]);
		      costheta= costheta/(r[sn[i][m]][sn[i][o]]*r[sn[i][o]][sn[i][p]]);
				test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][o]])/2.35;
		      ercos= err_cos(x[sn[i][m]],x[sn[i][o]],x[sn[i][p]],y[sn[i][m]],y[sn[i][o]],y[sn[i][p]],z[sn[i][m]],z[sn[i][o]],z[sn[i][p]],r[sn[i][m]][sn[i][o]],r[sn[i][o]][sn[i][p]],sigma_thet);
		      etotalen=etotale - e[sn[i][j]] - e[sn[i][l]] - e[sn[i][m]];
		      escattern = etotalen/(1+(etotalen/mec2)*(1-costheta));
		      escatter = etotalen - e[sn[i][o]];
		      
		      deltaescn = SQ(escattern)*ercos/mec2;
		      
		      deltaesc = sqrt((numn[i]+4)*SQ(eres));
		      
		      probtest2 = exp(-SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));
		       if(probtest2 < 0.0005) goto nextfifth;
		      cross1 = sig_compt(etotalen);
		      cross=sig_pair(etotalen);
		      coef1 = cross1/(cross1+cross+sig_abs(etotalen));
		      cross2 = sig_compt(escatter);
		      cross=sig_pair(escatter);
		      coef2 = cross2/(cross2+cross+sig_abs(escatter));
		      if(numn[i] == 5) {
			cross2=sig_abs(escatter);
			coef2= cross2/(cross2+cross+sig_compt(escatter));
			
		      }
		      lambda1=range_process(cross1);
		      lambda2=range_process(cross2);
		      probtest2=probtest1*probtest2*coef2*proba(lambda2,r_ge[sn[i][o]][sn[i][p]]);
		      
		      
		      for(q=0;q<numn[i];q++) {
			if(q!=j && q!=l && q!=m && q!=o && q!=p ) {
			  
			  costheta=(x[sn[i][p]]-x[sn[i][o]])*(x[sn[i][q]]-x[sn[i][p]])+(y[sn[i][p]]-y[sn[i][o]])*(y[sn[i][q]]-y[sn[i][p]])+(z[sn[i][p]]-z[sn[i][o]])*(z[sn[i][q]]-z[sn[i][p]]);
			  costheta= costheta/(r[sn[i][o]][sn[i][p]]*r[sn[i][p]][sn[i][q]]);
				test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][p]])/2.35;
			  ercos= err_cos(x[sn[i][o]],x[sn[i][p]],x[sn[i][q]],y[sn[i][o]],y[sn[i][p]],y[sn[i][q]],z[sn[i][o]],z[sn[i][p]],z[sn[i][q]],r[sn[i][o]][sn[i][p]],r[sn[i][p]][sn[i][q]],sigma_thet);
			  etotalen=etotale - e[sn[i][j]] - e[sn[i][l]] - e[sn[i][m]] - e[sn[i][o]];
			  escattern = etotalen/(1+(etotalen/mec2)*(1-costheta));
			  escatter = etotalen - e[sn[i][p]];
			  
			  deltaescn = SQ(escattern)*ercos/mec2;
			  
			  deltaesc = sqrt((numn[i]+5)*SQ(eres));
			  
			  probtest3=exp(-SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));
			  
			  cross1 = sig_compt(etotalen);
			  cross=sig_pair(etotalen);
			  coef1 = cross1/(cross1+cross+sig_abs(etotalen));
			  cross2 = sig_compt(escatter);
			  cross=sig_pair(escatter);
			  coef2 = cross2/(cross2+cross+sig_abs(escatter));
			  if(numn[i] == 6) {
			    cross2=sig_abs(escatter);
			    coef2= cross2/(cross2+cross+sig_compt(escatter));
			    
			  }
			  lambda1=range_process(cross1);
			  lambda2=range_process(cross2);
			  probtest3=probtest2*probtest3*coef2*proba(lambda2,r_ge[sn[i][p]][sn[i][q]]);
			 
			  
			  for(qq=0;qq<numn[i];qq++) {
			    if(qq!=j && qq!=l && qq!=m && qq!=o && qq!=p && qq != q ) {
			      
			      costheta=(x[sn[i][q]]-x[sn[i][p]])*(x[sn[i][qq]]-x[sn[i][q]])+(y[sn[i][q]]-y[sn[i][p]])*(y[sn[i][qq]]-y[sn[i][q]])+(z[sn[i][q]]-z[sn[i][p]])*(z[sn[i][qq]]-z[sn[i][q]]);
			      costheta= costheta/(r[sn[i][p]][sn[i][q]]*r[sn[i][q]][sn[i][qq]]);
					test_sigma_thet=sigma_thet*sqrt(0.1/e[sn[i][q]])/2.35;
			      ercos= err_cos(x[sn[i][p]],x[sn[i][q]],x[sn[i][qq]],y[sn[i][p]],y[sn[i][q]],y[sn[i][qq]],z[sn[i][p]],z[sn[i][q]],z[sn[i][qq]],r[sn[i][p]][sn[i][q]],r[sn[i][q]][sn[i][qq]],sigma_thet);
			      etotalen=etotale - e[sn[i][j]] - e[sn[i][l]] - e[sn[i][m]] - e[sn[i][o]] -e[sn[i][p]];
			      escattern = etotalen/(1+(etotalen/mec2)*(1-costheta));
			      escatter = etotalen - e[sn[i][q]];
			      
			      deltaescn = SQ(escattern)*ercos/mec2;
			      
			      deltaesc = sqrt((numn[i]+6)*SQ(eres));
			      
			      probtest4=exp(-SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));
			      
			      cross1 = sig_compt(etotalen);
			      cross=sig_pair(etotalen);
			      coef1 = cross1/(cross1+cross+sig_abs(etotalen));
			      cross2 = sig_compt(escatter);
			      cross=sig_pair(escatter);
			      coef2 = cross2/(cross2+cross+sig_abs(escatter));
			      if(numn[i] == 6) {
				cross2=sig_abs(escatter);
				coef2= cross2/(cross2+cross+sig_compt(escatter));
				
			      }
			      lambda1=range_process(cross1);
			      lambda2=range_process(cross2);
			      probtest4=probtest3*probtest4*coef2*proba(lambda2,r_ge[sn[i][q]][sn[i][qq]]);
			      probtest4=pow(probtest4,1/13.);
			      
			      if(probtest4 > mincosdif) {
				mincosdif=probtest4;
				interaction[i][0]=sn[i][j];
				interaction[i][1]=sn[i][l];
				interaction[i][2]=sn[i][m];
				interaction[i][3]=sn[i][o];
				interaction[i][4]=sn[i][p];
				interaction[i][5]=sn[i][q];	
				interaction[i][6]=sn[i][qq];
				
			      }
			      
			    }
			  }
			  
			  
			  if(numn[i]==6) {
			    probtest3=pow(probtest3,1/11.);
			    if(probtest3 > mincosdif) {
			      mincosdif=probtest3;
			      interaction[i][0]=sn[i][j];
			      interaction[i][1]=sn[i][l];
			      interaction[i][2]=sn[i][m];
			      interaction[i][3]=sn[i][o];
			      interaction[i][4]=sn[i][p];
			      interaction[i][5]=sn[i][q];		    
			      
			    }
			    
			  }
			}
			
		      }
		      
		      if(numn[i]==5) {
			
			probtest2=pow(probtest2,1/9.);
			
			if(probtest2 > mincosdif) {
			  mincosdif=probtest2;
			  
			  interaction[i][0]=sn[i][j];
			  interaction[i][1]=sn[i][l];
			  interaction[i][2]=sn[i][m];
			  interaction[i][3]=sn[i][o];
			  interaction[i][4]=sn[i][p];
			  
			  
			}
		      }
		      
		    }
		  nextfifth: ;
		  }
		  
		  if(numn[i]==4){
		    probtest1=pow(probtest1,1/7.);
		    if(probtest1 > mincosdif) {
		      mincosdif = probtest1;
		      
		      interaction[i][0]=sn[i][j];
		      interaction[i][1]=sn[i][l];
		      interaction[i][2]=sn[i][m];
		      interaction[i][3]=sn[i][o];
   		      
		    }
		  }
		}
	      nextfourth:;		
		
	      }
	      if(numn[i]==3) {
		probcomp=pow(probcomp,1/5.);
		if(probcomp > mincosdif) {
		  mincosdif = probcomp;
		  interaction[i][0]=sn[i][j];
		  interaction[i][1]=sn[i][l];
		  interaction[i][2]=sn[i][m];
 		}
		
	      }
	    }
	  nextthird:;
	  }
	  if(numn[i]==2) {
	    probtest=pow(probtest,1/3.);
	    if(probtest > mincosdif) {
	      mincosdif=probtest;
	      interaction[i][0]=sn[i][j];
	      interaction[i][1]=sn[i][l];
	    }
	  }
	}
   nextsecond:;  
      }
   
 
    }
    
  suite:if(numn[i]==1){
    mincosdif=minprobtrack; 
    interaction[i][0]=sn[i][0];
  }
	
    probtot[i]=mincosdif;
      if(flagpair[i]==1 && probtestpair> probtot[i]){ // pair production event and not good compton
          interaction[i][0]=interactionpair[i][0];
          probtot[i]=probtestpair;
          
      }

  
  }
  
  
  /***********************************************************************/
  /************ sort clusters according to figure of merit ****************/ 
  /***********************************************************************/
 
 
  // single interactions are awarded for the time being the minimum figure of merit 
  
  for(i=0;i<n;i++) {
    for(j=i+1;j<n;j++) {
      if(probtot[i] < probtot[j]) {
	swap(probtot, i, j);
	swap(angn,i,j);
	swap(et, i, j);
	swapi(flagu, i, j);
	swapi(flagpair, i, j);
	for (l=0;l<numn[i];l++)
	  temp1[l] = interaction[i][l];
	for(l=0;l<numn[j];l++)
	  interaction[i][l]=interaction[j][l];
	for(l=0;l<numn[i];l++)
	  interaction[j][l]=temp1[l];
	for(l=0;l<numn[i];l++)
	  temp1[l] = sn[i][l];
	for(l=0;l<numn[j];l++)
	  sn[i][l]=sn[j][l];
	for(l=0;l<numn[i];l++)
	  sn[j][l]=temp1[l];
        swapi(numn, i, j);      
	
      }
    }
  }
    
    
    
    

   
  for(i=0;i<n;i++) {
 
    if(flagu[i]==0) {
      for(k=0;k<numn[i];k++) {
	for(l=i+1;l<n;l++) {
	  
	  for(m=0;m<numn[l];m++) {
	    if(sn[i][k] == sn[l][m]){
	      
	      flagu[l]=1;
	      break;
	    }
	  }
	  
	}
	
      }   
      
    }
    
  }

  number_of_unflagged=0;
   for(i=0;i<n;i++){
    if(flagu[i]==0)number_of_unflagged++;

  }
  

	l1=0;

    
    
    for(j=0;j<n;j++){
        for(k=0;k<numn[j];k++){
            iseg[sn[j][k]]=1;
        }
    }
    
    
  for(j=0;j<n;j++){
  
      
	  if(numn[j]>1 && et[j]*1e3<20000 && probtot[j]<1)specvsprob->Fill(probtot[j],et[j]*1e3);
		  
    if(flagu[j]==0 && numn[j]>1) {
     
      dx=x[interaction[j][0]]-xsource;
      dy=y[interaction[j][0]]-ysource;
      dz=z[interaction[j][0]]-zsource;
      dd = SQ(dx)+SQ(dy)+SQ(dz);
      if(dd==0)etot=et[j];
      if(dd!=0){
	cosTheta = (dx*vrecoilx + dy*vrecoily + dz*vrecoilz)/sqrt(dd);
	etot=et[j]*(1-voverc*cosTheta)/betafactor;
      }
     
	
		
			
      if(probtot[j] >=minprobtrack){
          

          

		  gamma[l1]=etot*1e3;
          dist[l1]=dd;
          nmult[l1]=numn[j];
          firstx[l1]=x[interaction[j][0]]-xsource;
          firsty[l1]=y[interaction[j][0]]-ysource;
          firstz[l1]=z[interaction[j][0]]-zsource;
          
          secondx[l1]=x[interaction[j][1]]-x[interaction[j][0]];
          secondy[l1]=y[interaction[j][1]]-y[interaction[j][0]];
          secondz[l1]=z[interaction[j][1]]-z[interaction[j][0]];
          
          dist2[l1]=SQ(secondx[l1])+SQ(secondy[l1])+SQ(secondz[l1]);
          
		  l1++;
		  	if(etot*1e3<20000)spectrumout->Fill(etot*1e3);
		  if(et[j]*1e3<20000){
			  spectrumout_nodop->Fill(et[j]*1e3);
			
		  }
		  yfl=drand48();
		  k=(int) (1e3*etot+yfl);
		  if(k<8192)spn_data[k][1]++;
        if(k<8192)spn_data[k][2]++; // no sing

		
		  
		
    
      itot++;
	       }
		

 
    }
  }
  
    
  i=0;
  
  for(j=0;j<n;j++) {
    if(numn[j]==1 && flagu[j]==0)i++;
    
  }
  
 
  
  /**********************************/
  /* test of direct photopic events */
  /**********************************/
  
	l2=0;
  
	if(nodirect==0) { //  treat single interactions
  j=i;
  
  for(l=0;l<n;l++) {
	  
    if(numn[l]==1 && flagu[l]==0) {
      rmin=50;
	
      for(m=0;m<nb_int;m++) {
	if(m!=sn[l][0]){
	  if(r[sn[l][0]][m] < rmin)
	    rmin = r[sn[l][0]][m];
	}
	
      }
     
	cross1=sig_abs(et[l]);
	
	cross = sig_pair(et[l]);
	cross2=cross1+cross+sig_compt(et[l]);
	lambda=range_process(cross2);
			  probcomp = proba(lambda,r_ge[sn[l][0]][sn[l][0]])*cross1/cross2;//*(1-log(cross1/cross2));

        dd_ge=r_ge[sn[l][0]][sn[l][0]];
		

      	
	
		if(rmin > photodistmax) {
			if(sqrt(probcomp)<1 && et[l]*1e3<20000)specvsprobsing->Fill(sqrt(probcomp),et[l]*1e3);
            specvsrange->Fill(r_ge[sn[l][0]][sn[l][0]],et[l]*1e3);
            
 
	  dx=x[interaction[l][0]]-xsource;
	  dy=y[interaction[l][0]]-ysource;
	  dz=z[interaction[l][0]]-zsource;
	  dd = SQ(dx)+SQ(dy)+SQ(dz);
			if(dd==0){
				dd=d_res*d_res;
				cosTheta = (dx*vrecoilx + dy*vrecoily + dz*vrecoilz)/sqrt(dd);
				etot=et[l]*(1-voverc*cosTheta)/betafactor;
				
			}
	  if(dd!=0){
	    cosTheta = (dx*vrecoilx + dy*vrecoily + dz*vrecoilz)/sqrt(dd);
	    etot=et[l]*(1-voverc*cosTheta)/betafactor;
	  }

			
			
            
		  if(sqrt(probcomp) > minprobsing) {
			  
          
            	
			  gamma[l1]=etot*1e3;
              dist[l1]=dd;
            
              firstx[l1]=x[interaction[l][0]]-xsource;
              firsty[l1]=y[interaction[l][0]]-ysource;
              firstz[l1]=z[interaction[l][0]]-zsource;
              nmult[l1]=1;
			  l1++;
              
			  gamma1[l2]=etot*1e3;
			  l2++;
			  if(etot*1e3<20000){
				  spectrumout->Fill(etot*1e3);
				   spectrumout_sing->Fill(etot*1e3);
			  }
			    
			  if(et[l]*1e3<20000){
				  spectrumout_nodop->Fill(et[l]*1e3);
				  spectrumout_sing_nodop->Fill(et[l]*1e3);
			  }
			 
			  
			  
			  yfl=drand48();
			  k=(int) (1e3*etot+yfl);
			  if(k<8192)spn_data[k][1]++;
             
			
	  
	
	  itot++;

	  j--;
	 		}
		  
		      }
		


    }
    
    
  }
	}
    
    
  // analysis of l1 tracked gammas
	
	for(i=0;i<l1;i++){
        
		
		for(k=i+1;k<l1;k++){
            
       
            
			if(gamma[i]<2000 && gamma[k]<2000){

                
                coinc->Fill(gamma[i],gamma[k]);
                coinc->Fill(gamma[k],gamma[i]);
			}
		}
		
	}
	
   
	
	
	
    
    
 
  
  if(readflag == 0) {
    counter =1;
	  i=0;
      ecalore=0;
      esum=0;
    goto readfile ; /* go to read next event */
  }
  
  /**********************************************************************/
  /************* end of tracking procedure *****************************/
  /*********************************************************************/
  
  cout<<"\r"<<nintprev<<" events "<<endl;
   

 
   printf("total numb of gamma-rays shot = %d \n",nbrtotgam);
  printf("total numb of gamma-rays with interaction = %d \n",nshoot);
  printf("total numb of gamma-rays detected = %d \n",itot);
	
	  

 

		

    wsec(rwfile,spn_data,YSPNT);
	
	
  fclose(fp);
    
	file->Write();
	file->Close();
        
    } // found data file
    } // found input LookUpTable
    } // found input file OFTinput
  
   
	
	return 0;
} /* end main */



double sig_compt(double E)
     
{
  /* sigma = 2 pi r_0^2 *[((1+gamma)/gamma^2 *(2(1+gamma)/1+2gamma - ln(1+2gamma)/gamma)) + ln(1+2gamma)/2gamma - 1+3gamma/(1+2gamma)^2] */
  // fits geant data very well   
  
  double temp;
  double temp0;
  double temp1;
  double temp2;
  double temp3;
  double gamma;
  
  temp0 = 1e4*2*PI*SQ(r_0)*Z_ge; /* in cm2/atom */
  gamma = E/mec2;
  
  temp1 =1+gamma;
  temp2= 1 +(2*gamma);
  temp3= 1+(3*gamma);
  temp = (temp1/SQ(gamma)) * (((2*temp1)/temp2) - (log(temp2)/gamma));
  temp=temp +(log(temp2)/(2*gamma) - (temp3/SQ(temp2))); 
  temp=temp*temp0;
  
  return temp ;
  
}

double sig_abs(double E)
     
{
  
  /* sigma = 4 alpha^4*sqrt(2)*Z^5 phi_0*(E/mec2)^(-7/2),phi_0=8/3 pi r_0^2 */
  /* sigma abs K+L = 9/8* sigma_Kshell */
  
  double temp;
  double gamma;
  double hnu_k;
  
  
  hnu_k = SQ(Z_ge - 0.03)*mec2*SQ(alpha)/2;
  gamma = CB(E/mec2)*CB(E/mec2)*E/mec2;
  gamma = sqrt(gamma);
  
  temp = 4*CB(alpha)*alpha*1.4142*6.651e-25*CB(Z_ge)*SQ(Z_ge);
  
  temp = sqrt(E/mec2) * temp/gamma; /* en cm2/atom */
  
  // not well suited for energies below 20 keV 
  // removed the 1.125 factor and added sqrt(E/mec2) to fit data 
  
  if(E < 0.025) {
    
    temp = 2.2*pow((hnu_k/E),2.6666)*6.3e-18/SQ(Z_ge);
    if(E<0.0111)
      temp = temp/8.5;
    
  }
  return temp;
  
}


double range_process(double sig)
{
  
  /* SIGMA MACRO = 1/lambda = 6.022e23 * rho_ge/A_ge * sigma (cm2/atom) */
  
  
  double temp;
  
  temp = (sig * N_av * rho_ge) / A_ge;
  temp = 1/(temp); /* in cm */
  
  return temp ;
  
  
  
}

double proba(double range, double distance)
     
{
  double temp;
  double nlambda;
  
  nlambda = distance/range;
  temp = exp(-nlambda);
  
  return temp;
  
  
  
}

void swap(double v[],int m, int l)
{
  double temp;
  
  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
  
}

void swapi(int v[],int m, int l)
{
  int temp;
  
  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
  
}

double sig_pair(double E)
{/* fitted in 1e-24cm2 units for Ge from atomic data tables 1970 7 page 590 */
  
  double temp;
  
  temp = 0.792189*log(E+0.948261-1.1332*E+0.15567*SQ(E));
  
  if(E<1.022)
    temp = 0;
  
  if((E < 1.15) & (E >= 1.022))
    temp=(1-((1.15-E)/0.129))*7.55e-28;
  
  
  return temp*1e-24;
}
double err_cos(double xa, double xb, double xc, double ya, double yb, double yc, double za, double zb, double zc, double rab, double rbc, double err)
     /* error on AB.BC/AB x BC */
{
  double prod;
  double dcosdxa,dcosdxb,dcosdxc;
  double dcosdya,dcosdyb,dcosdyc;
  double dcosdza,dcosdzb,dcosdzc;
  
  double temp;
  
  prod = (xa-xb)*(xb-xc) + (ya-yb)*(yb-yc)+ (za-zb)*(zb-zc);
  
  
  dcosdxa = -(xc-xb)/(rab*rbc) -0.5*(2*xa-2*xb)*prod/(CB(rab)*rbc);  
  dcosdxb = (xa-2*xb+xc)/(rab*rbc) -0.5*(2*xb-2*xa)*prod/(CB(rab)*rbc) +0.5*(2*xc-2*xb)*prod/(CB(rbc)*rab);
  dcosdxc = (xb-xa)/(rab*rbc) - 0.5*(2*xc-2*xb)*prod/(CB(rbc)*rab);
  
  dcosdya = -(yc-yb)/(rab*rbc) -0.5*(2*ya-2*yb)*prod/(CB(rab)*rbc);
  dcosdyb = (ya-2*yb+yc)/(rab*rbc) -0.5*(2*yb-2*ya)*prod/(CB(rab)*rbc) +0.5*(2*yc-2*yb)*prod/(CB(rbc)*rab);
  dcosdyc = (yb-ya)/(rab*rbc) - 0.5*(2*yc-2*yb)*prod/(CB(rbc)*rab);

  dcosdza = -(zc-zb)/(rab*rbc) -0.5*(2*za-2*zb)*prod/(CB(rab)*rbc);
  dcosdzb = (za-2*zb+zc)/(rab*rbc) -0.5*(2*zb-2*za)*prod/(CB(rab)*rbc) +0.5*(2*zc-2*zb)*prod/(CB(rbc)*rab);  
  dcosdzc = (zb-za)/(rab*rbc) - 0.5*(2*zc-2*zb)*prod/(CB(rbc)*rab);

  temp = SQ(dcosdxa)+SQ(dcosdxb)+SQ(dcosdxc);
  
  temp=temp+SQ(dcosdya)+SQ(dcosdyb)+SQ(dcosdyc);
  
  temp=temp+SQ(dcosdza)+SQ(dcosdzb)+SQ(dcosdzc);
  
  temp=sqrt(temp);
  
  temp=temp*err;
  


  return temp;
  } 

int packpoints(int number, double posx[], double posy[], double posz[], double posxrel[], double posyrel[], double poszrel[], double transx[], double transy[], double transz[], double energy[],int intnb[], int intorder[],int nseg[], int ndet[], double resolution)
{
    int i,j,n,l,jp[10000],jjp[10000],ip[10000],iip[10000];
    double rpack,esum;
    
    l=0;
    for (i=0;i<number;i++) {
        for (j=i+1;j<number;j++) {
            
            rpack= sqrt(SQ(posx[i]-posx[j])+SQ(posy[i]-posy[j])+SQ(posz[i]-posz[j]));
            
            if(rpack < resolution && ndet[i]==ndet[j]) {
                l++;
                jp[l]=j;
                jjp[l]=j;
                ip[l]=i;
                iip[l]=i;
            }
        }
    }
    
    /* check if couples have already been packed by previous couples */
    
    for (i=1;i<=l;i++) {
        for(j=i+1;j<=l;j++) {
            if(ip[i] == ip[j]){ //&& ip[i] != -1) {
                for(n=j+1;n<=l;n++) {
                    if(ip[n] == jp[i] && jp[n]==jp[j]) {
                        
                        
                        iip[n] = -1;
                        jjp[n]= -1;
                        
                        
                    }
                    if(ip[n] == jp[j] && jp[n] == jp[i]) {
                        
                        iip[n] = -1;
                        jjp[n] = -1;
                    }
                    
                }
                
            }
            
        }
        
    }
    
    for(n=1;n<=l;n++) {
        if(iip[n]==-1) ip[n] = -1;
        if(jjp[n]==-1) jp[n]=-1;
        
    }
    
    
    
    for (j=1;j<=l;j++) {
        
        
        
        
        for (i=0;i<number;i++) {
            
            if(ip[j] == i && jp[j] != i) {
                
                esum=energy[i]+energy[jp[j]];
                
                /* new position pondered by energie */
                
                
                posx[i]=((posx[i]*energy[i])+(posx[jp[j]]*energy[jp[j]]))/esum;
                posy[i]=((posy[i]*energy[i])+(posy[jp[j]]*energy[jp[j]]))/esum;
                posz[i]=((posz[i]*energy[i])+(posz[jp[j]]*energy[jp[j]]))/esum;

		//#ifdef POSRESMAP
                posxrel[i]=((posxrel[i]*energy[i])+(posxrel[jp[j]]*energy[jp[j]]))/esum;
                posyrel[i]=((posyrel[i]*energy[i])+(posyrel[jp[j]]*energy[jp[j]]))/esum;
                poszrel[i]=((poszrel[i]*energy[i])+(poszrel[jp[j]]*energy[jp[j]]))/esum;
		//#endif
		
		energy[i]=esum; /* put 2 energies into 1 */
                
                
                swap(energy,number-1,jp[j]); /* put unused energy at the end of the list */
                swap(posx,number-1,jp[j]);
                swap(posy,number-1,jp[j]);
                swap(posz,number-1,jp[j]);
                swap(posxrel,number-1,jp[j]);
                swap(posyrel,number-1,jp[j]);
                swap(poszrel,number-1,jp[j]);
                swap(transx,number-1,jp[j]);
                swap(transy,number-1,jp[j]);
                swap(transz,number-1,jp[j]);
                swapi(intnb,number-1,jp[j]);
                swapi(intorder,number-1,jp[j]);
                swapi(nseg,number-1,jp[j]);
                swapi(ndet,number-1,jp[j]);
                
                
                
                for (n=j+1;n<=l;n++) { 	    
                    if(ip[n] == jp[j]) { /* if the one just packed needs to be packed */
                        
                        ip[n]= i;
                        
                    }
                    
                    if(jp[n] == jp[j]) { /* if the one just packed needs to be packed */
                        jp[n]= i;
                    }
                    if(ip[n] == number-1) {/* if end of list needs to be packed */
                        
                        
                        ip[n] = jp[j];
                        
                    }
                    
                    if(jp[n]== number-1) {/* if end of list needs to be packed */
                        
                        jp[n]=jp[j];
                    }
                    
                }
                number-=1;  /* decrement the number of interactions */	
            }	
            
        }      
    }
    
    
    return number;
    
    
    
}


void smearpoints(int number, double energy[], double posx[], double posy[], double posz[],double errorfc[],double step, int nbsteps, double uncertainty)
{
  int i,j;


  double tru,err,ener_u,ener_p;

 for (i=0;i<number;i++) {

   ener_u = sqrt(1 + energy[i]*3.7)/2.35; // fwhm/2.35 in keV
   //ener_p = uncertainty*sqrt(0.1/energy[i]); // fwhm/2.35 in cm
   ener_p = uncertainty*( (2.7+6.2*sqrt(0.1/energy[i]))/4.66 ); //from NIMA638
   
   tru=drand48();
    if (tru < errorfc[0])
	err = step;
      for(j=0;j<=nbsteps;j++){
	if(errorfc[j] > tru ){
	  err = j*step;
	  break;
	}
      }      

      tru=drand48();
    
      if(tru <= 0.5)err=-err;
  
      
      energy[i] = energy[i]+(err*ener_u)/1000; // add error in MeV
      if(energy[i]<0)energy[i] = fabs((err*ener_u)/1000);
   
      tru = drand48();
      // printf("random = %f \n",tru);
      
      if (tru < errorfc[0])
	err = step;
      for(j=0;j<=nbsteps;j++){
	if(errorfc[j] > tru ){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48();
    
      if(tru<=0.5) /* if random <=0.5 dx=>-dx */
	err=-err;      
      posx[i]=posx[i]+ener_p*err;  
      
      tru = drand48();
      if (tru < errorfc[0])
	err = step;
      for(j=0;j<=nbsteps;j++){
	
	if(errorfc[j] > tru){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48();
     
      if(tru<=0.5)
	err=-err;      
      posy[i]=posy[i]+ener_p*err;
      
      tru=drand48();
      if (tru < errorfc[0])
	err = step;
      for(j=0;j<nbsteps;j++){
	
	if(errorfc[j] > tru){
	  err = j*step;
	  break;
	}
      }      
      tru=drand48();
     
      if(tru<=0.5)
	err=-err;      
      posz[i]=posz[i]+ener_p*err;
      
    }


	



}


void wsec ( char file_l[] , int matrix_l[][MAXYSPN] , int yspn_l )
{
	int        matseg_l [XSPN];
	FILE       *fp_l;
	int        i_l , ii_l;
	
	
	
	printf ( "writing .spn file : %s\n" , file_l );
	
	if ( ( fp_l = fopen ( file_l , "w" )) == NULL )
	{
		printf ( "error : opening .spn file : %s\n" , file_l );
		exit ( 1 );
	}
	
	
	for ( ii_l = 1 ; ii_l <= yspn_l ; ii_l++ )
	{
		for ( i_l = 0 ; i_l <= XSPN-1 ; i_l++ )
			matseg_l[i_l] = (int) matrix_l[i_l][ii_l-1];
		fwrite ( matseg_l , sizeof (matseg_l), 1 , fp_l );
	}
	
	fclose ( fp_l );
	
	printf ( "end   writing .spn file : %s\n" , file_l );
}



