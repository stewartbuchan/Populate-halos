#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stack>
#include <list>
#include <vector>
#include <iterator>
#include <random>
#include <cmath>


using namespace std;


/*cosmological parameters*/
double omega_m=0.25;
double omega_l=1.-omega_m;
double h=0.73;
double spectral_index=1.0;
double omega_baryon=0.045;
double sigma_8=0.9;
double omega_0=1.;
double domega=0.1;
double dc_0=1.69;
double G=4.32; //  km^2 * kpc^3 * Mpc^-2 * s^-2 * Msun^-1
double forb=0.;
bool includemerger=true;
bool includeshutdown=false;
bool includediskregrowth=true;
bool includeclumpy_accretion=false;
bool includebathtub_model_accretion=false;
bool includegasdissipation=false;
bool includeconcentration=false;
bool includeexpansion=true;
//bool includeDMsubstructure=false;
bool reinitialiseDM=false;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0,1.0);

class satillite_halo;
typedef std::list<satillite_halo> listgals;  
typedef std::list<listgals> listlistgals;



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> interpol(std::list<double> x1, std::list<double> x2, std::list<double> y1)
{
//linear interpolation; x1 y1 is library and x2 is input.

cout << "asdfasf  " <<  x1.size() << "  " << x2.size() << "  " << y1.size() << endl; 

std::list<double> y2;
std::list<double>::iterator itx1=x1.begin(),itx2=x2.begin(),ity1=y1.begin();
double prex,prey;
while(itx2!=x2.end()){
cout <<"blub_blub  " << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;
while(*itx1 > *itx2&&itx2!=x2.end()){cout <<  *itx1 << "  " << *ity1 << "  " << *itx2 << endl;y2.push_back(((*itx2-prex)/(*itx1-prex)*(*ity1-prey))+prey);itx2++;}
prex=*itx1;prey=*ity1;
itx1++;ity1++;
}
cout << "sdfsfg  " << y2.front()<<endl;
return y2;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double Mdm_to_Mstellar(double logM, double logZ)
{
// Z = 1+z
double Z=pow(10,logZ);
float m10=11.59, m11=1.195, n10=0.0351, n11=-0.0257, beta10=1.376, beta11=-0.826, gamma10=0.608, gamma11=0.329;
double logM1=m10+m11*(Z-1)/Z;
double n=n10+n11*(Z-1)/Z;
double beta=beta10+beta11*(Z-1)/Z;
double gamma=gamma10+gamma11*(Z-1)/Z;
double logm=logM+log10(2*n)-log10(pow(10,gamma*(logM-logM1))+pow(10,-beta*(logM-logM1)));
//cout <<logm<<endl;
return logm;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double get_Mgas(double logMstar, double logZ,double maxgas)
{
//Stewart et al 2009 apj 702:307
double Z=pow(10,logZ);
double alpha=0.59*pow(Z,0.45);
double logMgas=logMstar+log10(0.04)-alpha*(logMstar-log10(4.5*pow(10,11)));
//cout << "GAS!  " << logMgas << "  " << maxgas<<endl;
if(logMgas>maxgas){/*cout <<"ops"<<endl;*/ logMgas=maxgas;}
return logMgas;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



/*

double get_Rvir(double Mhalo,double Z){
//Guo et al 2010
//Z=log(1+Z) Mhalo=log(Mhalo)
double H_Z=(h*h*10000)*(omega_m*pow(10.,3.*Z)+omega_l);
return pow(pow(10,Mhalo)*G/100./H_Z,(1./3.));
}
*/


double get_Rvir(double Mhalo,double Z){

double omegaZ=omega_m*pow(10,3.*Z)/(omega_m*pow(10,3*Z)+omega_l);
double dc=10.*M_PI*M_PI+82*(1.-omegaZ)-39.*(1.-omegaZ)*(1.-omegaZ);
double alz=31.*pow(omega_m/omegaZ*dc/(18.*M_PI*M_PI),(-1./3.))*(7./pow(10,Z));
return (0.7/h)*alz*pow(10,(Mhalo-12.)/3.)/1000;

}


double get_Rdisk(double Mdisk,double Z){
return (-1)-0.4*Z+0.14*Mdisk+(0.39-0.14)*log10(1+pow(10,(Mdisk-10.59988)));
}

double get_Rbulge(double M1,double R1,double M2,double R2, double fc,double gas){
//double c=0.5;
//double f=forb;
double f0=0.25;

if(M2<11){R2*=2;}

double AA,BB;
if(R1<-3){AA=0.;}
else{AA=pow(10.,(2.*M1-R1));}
if(R2<-3){BB=0.;}
else{BB=pow(10.,(2.*M2-R2));}
//cout << "testR  " << AA << "  " << BB << "  " << endl;
double p=AA+BB+fc*pow(10.,(M1+M2))/(pow(10.,R1)+pow(10.,R2));
double r=(pow(10.,M1)+pow(10.,M2))*(pow(10.,M1)+pow(10.,M2))/p;
if(includegasdissipation){
double fgas=pow(10.,gas)/(pow(10.,M1)+pow(10.,M2));
r=r/(1.+fgas/f0);
}
return log10(r);
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_T_dynamical_friction(double Mh,double Ms,double Z){
//using formula in Boylan & kolchin 2008
double A=0.9,B=1.,C=0.6,D=0.1; //revised values from McCavana et al 2012
double T_dynamical=(0.1*9.78/(h*sqrt(omega_m*pow(10,3.*Z)+omega_l)));  // UNITS!!!!
cout << "Td  " << Mh << "  " << Ms << "  " << Z << "  " << T_dynamical<<endl;
double M_ratio=pow(10,(Mh-Ms));
double coeff=A*T_dynamical*pow(M_ratio,B)/log(1+(M_ratio));
double eta=0.23*distribution(generator)+0.5;
//cout << eta<<endl;
if(eta>=1.){eta=0.999;}
if(eta<=0.){eta=0.001;}
double eps=sqrt(1-eta*eta);
double rvir=get_Rvir(Mh,Z)/1000;
double rper=rvir*pow(eta,2.17);
double rc=rper/(1.-eps);
//cout << "TDF" << "  " <<T_dynamical <<endl;
double T_dynamical_friction=coeff*exp(C*eta)*pow(rc/rvir,D);
cout << "TDF" << "  " << Mh << "  " << Ms << "  " << pow(10,Z)-1. << "  " <<  rc/rvir << "  " << T_dynamical_friction << "  " << T_dynamical << "  " << T_dynamical_friction/T_dynamical <<endl;

return T_dynamical_friction;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_time_from_redshift(double Z){
//returns time in Gyrs from log(1+Z)
if(omega_m==1){return 9.777505969*(2./3./h)*(pow(10.,(-Z*3./2.)));}
else{return 9.777505969*(2./3/h/sqrt(1.-omega_m))*asinh(sqrt((1.-omega_m)/omega_m)/pow(10,(Z*(3./2))));}
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_M_halo_shutdown(double logZ){

//see Cattaneo+06 for details

double logZc=log10(3.+1);
double logMshock=12.;
double K;

if(logZ<logZc){K=1.3*(pow(10,logZ)-pow(10,logZc));}
else{K=0;}

return logMshock+K;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


// sub halo class
class satillite_halo {
double mass; //mass of DM in halo
double redshift; //redshift of impact
double Tdy; //dynamical time - from impact to merger
double stellar_bulge; //mass of stars in bulge
double stellar_disk; //mass of stars in disk
double gas; //mass of gas in disk - assuming no gas in bulge
double Rhl_disk; //half light radius of disk
double Rhl_bulge; //hald light redius of bulge
bool had_major_merger;
bool shutdown;
bool hadextememerger;
bool cgt4;
float test;
double mu; // mass ratio of merger


public:
satillite_halo(double M, double Z,double new_stellar_bulge, double new_stellar_disk, double new_gas){mass=M;gas=gas;redshift=Z;stellar_bulge=new_stellar_bulge;stellar_disk=new_stellar_disk;gas=new_gas;shutdown=false;had_major_merger=false;hadextememerger=false;cgt4=false;test=0.1;}
double return_mass(){return mass;}
void put_DM(double new_DM){mass=new_DM;}
void add_DM(double new_DM){mass=log10(pow(10,mass)+pow(10,new_DM));}
void subtract_DM(double new_DM){mass=log10(pow(10,mass)-pow(10,new_DM));}
double return_gas(){return gas;}
void put_gas(double new_gas){gas=new_gas;}
void add_gas(double add_mass){gas=log10(pow(10,gas)+pow(10,add_mass));}
double return_stellar_bulge() {return stellar_bulge;}
void add_stellar_bulge(double add_mass){stellar_bulge=log10(pow(10,stellar_bulge)+pow(10,add_mass));}
double return_stellar_disk() {return stellar_disk;}
void add_stellar_disk(double new_stellar){stellar_disk=log10(pow(10,stellar_disk)+pow(10,new_stellar));}
double return_stellar_total(){return(log10(pow(10,stellar_disk)+pow(10,stellar_bulge)));}
double return_baryons() {return log10(pow(10,stellar_bulge)+pow(10,stellar_disk)+pow(10,gas));}

double return_Rhl_bulge() {return Rhl_bulge;}
double return_Rhl_disk() {return Rhl_disk;}
void put_Rhl_bulge(double new_Rhl_bulge){Rhl_bulge=new_Rhl_bulge;}
void put_Rhl_disk(double new_Rhl_disk){Rhl_disk=new_Rhl_disk;}
double return_Rhl(){return (Rhl_bulge*pow(10,stellar_bulge)+Rhl_disk*(pow(10,stellar_disk)+pow(10,gas)))/(pow(10,stellar_bulge)+pow(10,stellar_disk)+pow(10,gas));}


double return_redshift() {return redshift;}
void put_Z(double Z){redshift=Z;}
//double return_Tdy() {return 10000000000000000.;}
//double return_Tdy() {return 0.;}
double return_Tdy() {return Tdy;}
double put_Tdy(double nTdy){Tdy=nTdy;}
bool return_had_major_merger(){return had_major_merger;}
bool return_shutdown(){return shutdown;}
bool return_extrememerger(){return hadextememerger;}
bool put_had_major_merger(){had_major_merger=true;}
void put_shutdown(){shutdown=true;}
void put_extrememerger(){hadextememerger=true;}
void put_mu(double new_mu){mu=new_mu;}
double return_mu(){return mu;}
bool return_cgt4(){return cgt4;}
void put_cgt4(){cgt4=true;}

void output_vars(){cout<< mass << "  " << redshift << "  " << gas << "  " << stellar_bulge << "  " << stellar_disk << "  " << return_stellar_total() << "  " << Rhl_disk << "  " << Rhl_bulge << "  " << return_Rhl() << "  " << had_major_merger << "  " << shutdown << "  " << hadextememerger << "  " << cgt4 << endl;}
};

// accreted mass class class
class accreted_mass {
double mass;
double redshift;
bool addmass;
public:
accreted_mass(double M, double Z,bool read_add){mass=M;redshift=Z;addmass=read_add;}
double return_mass() {return mass;}
double return_redshift() {return redshift;}
double return_addmass(){return addmass;}
};

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


satillite_halo clumpyaccretion(satillite_halo cent,double Z1, double Z2,double Z,double gas_frac_crit){

double k=1;

if(pow(10,(cent.return_gas()-cent.return_baryons()))>gas_frac_crit){
double dgasdt=k*25.*pow(10,cent.return_stellar_disk()-11.)*pow(10,1.5*(Z-log10(3))); //rate of change of gas mass in solar mass/yr

/////check this!!! which frame is dt in? and does it matter? 
double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod where gas is accreting in Gyr

double dgas=log10(dgasdt*dt)+9; //log(change in gas mass) in log(solar mass)
if(dgas>cent.return_gas()){dgas=cent.return_gas();} //if there isn't enough gas in disk; use it all up!!

//cent.put_Rhl_bulge(log10(pow(10,(cent.return_Rhl_bulge()+cent.return_stellar_bulge()))/(pow(10,cent.return_stellar_bulge())+pow(10,dgas))));

cout << "CLUMPS  " << pow(10,cent.return_redshift())-1. << "  " << cent.return_mass() << "  " <<  cent.return_stellar_disk() << "  " << cent.return_gas() << "  " << dgas << "  " << dgas/cent.return_stellar_disk() <<endl;

//double  rb_old=cent.return_Rhl_bulge();

double rclump=cent.return_Rhl_disk()*M_PI/4.*cent.return_gas()/cent.return_baryons();

if(cent.return_stellar_bulge()<3.){cent.put_Rhl_bulge(rclump);}
else{cent.put_Rhl_bulge(get_Rbulge(cent.return_stellar_bulge(),cent.return_Rhl_bulge(),dgas,rclump, 2.0/0.5,cent.return_gas()));}

/*
double logaa=2.*cent.return_stellar_total();
double bb=pow(10,(2.*cent.return_stellar_bulge()-cent.return_Rhl_bulge()));
double cc=pow(10,(2.*cent.return_stellar_disk()-cent.return_Rhl_disk()));
double dd=pow(10,cent.return_stellar_bulge()+cent.return_stellar_disk())/(pow(10,cent.return_Rhl_bulge())+pow(10,cent.return_Rhl_disk()));
cent.put_Rhl_bulge(logaa-log10(bb+1.*cc+(2./0.5)*dd));
*/

cent.add_stellar_bulge(dgas);
cent.add_gas(-dgas);
//cout << "clumpyaccretion  " << dgasdt << "  " << dgas << "  " << rb << "  " << cent.return_Rhl_bulge() << "  " <<  cent.return_mass() << "  " << Z << endl;
}

return cent;

}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


void get_ssfr(satillite_halo cent, double Z1,double Z2){
double logA11=log10(0.0324)+11, alpha=3.45, beta=-0.35; 
double dmdt=pow(10,logA11+(beta+1)*(cent.return_stellar_disk()-11)+alpha*Z1);
double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod in Gyr //again check!!
double dm=log10(dmdt*dt);
cent.add_gas(-dm);
cent.add_stellar_disk(dm);
cent.put_Rhl_disk(get_Rdisk(cent.return_stellar_total(),Z1));
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

void bathtub_model_accretion(satillite_halo cent, double Z1,double Z2){
double beta=0.14,alpha=2.4,mu=0.54,epsilon=0.02;
double fga=1., eta=1.;
double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod where gas is accreting in Gyr
double td=0.071*get_time_from_redshift(Z1);

double Mdotb=80+pow(10,(1+beta)*cent.return_mass()-12)*pow(10,mu*Z1)/3*omega_baryon/omega_m;
double Mdotsf=cent.return_gas()*epsilon/td;
double Mdotg=fga*Mdotb-(mu+eta)*Mdotsf;
double Mdots=(1.-fga)*Mdotb+mu*Mdotsf;

cent.add_gas(log10(dt*pow(10,9)*Mdotg));
cent.add_stellar_bulge(log10(dt*pow(10,9)*Mdots));
cout << "bathtub_model_accretion  " << Mdots << "  "  << cent.return_mass() << "  " << Z1 <<endl;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double starburst(satillite_halo central, satillite_halo satillite){
if(satillite.return_baryons()>central.return_baryons()){return log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))+log10(0.54)+0.7*(central.return_baryons()-satillite.return_baryons());}
else{return log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))+log10(0.54)+0.7*(satillite.return_baryons()-central.return_baryons());}
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

void disk_regowth(satillite_halo central, satillite_halo satillite, double *gas, double *stellar_bulge,double *r_bulge, double *r_disk){

double newstars=starburst(central,satillite);
cout << "gas  " << *gas << endl;
if(newstars>=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))){newstars=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()));*gas=-300;}
else{*gas=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())-pow(10,newstars));}
*stellar_bulge=log10(pow(10,central.return_baryons())+pow(10,satillite.return_baryons())-pow(10,*gas));
*r_bulge=get_Rbulge(central.return_baryons(),central.return_Rhl(),satillite.return_baryons(),satillite.return_Rhl(),forb/0.5,*gas);  // // need to think of a better way of doing this!!!
r_disk=r_bulge; // // and this line!!
cout << "gas  " << *gas << "  " << newstars  << endl;

}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double get_concentration(double logM, double logZ){

double z=pow(10,logZ)-1.;
double a=0.520+(0.905-0.520)*exp(-0.617*pow(z,1.21));
double b=-0.101+0.026*z;

double logc200=a+b*(logM-12.);
return pow(10,logc200);
}



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


// // // this bit does the initializing either as disks/bulges or finds a galaxy from the catallogue \\  \\  \\

void find_galaxy(double DM, double Z, listlistgals *mergercat, std::stack<satillite_halo> *stack_Shalos){

	
	bool onlydisk=true;
	bool onlybulge=false;
	bool findfromcatalogue=false;
	

	listlistgals::iterator itZ=mergercat->begin();
	listgals::iterator itH;
    	std::list<double> tmpDM,tmpG,tmpB,tmpD,newDM,tmpRD,tmpRB;
	double stars;
	bool foundhalo=false;
	cout << "new_satillite  " << DM << "  " << Z<<endl;
	
	std::advance(itZ,(int)(Z*100));
	cout << "ok1\n";
//	cout << DM << "  " << Z << "  " << mergercat->size() << "  " << itZ->size() << "  " << itZ->back().return_mass() << "  " << itZ->front().return_mass() << endl;
	
	if(findfromcatalogue){
	if(itZ->size()>6){ //maybe getrid of this condition?
	if(itZ->back().return_mass()<DM){
	cout << "ok2  " << itZ->back().return_mass() << "\n";
	if(itZ->front().return_mass()>DM){
	cout << "ok3  " << itZ->front().return_mass() << "\n";	
		foundhalo=true;
		itH=itZ->begin();
		tmpDM.clear();tmpG.clear();tmpB.clear();tmpD.clear();tmpRB.clear();tmpRD.clear();newDM.clear();
		newDM.push_front(DM);
		cout << itZ->size() << "  " << mergercat->size() << endl;
		while(itH!=itZ->end()){
	//cout<<"FOUND IT !!!!!!!!!!!!!!!!!!!\n";
		tmpDM.push_back(itH->return_mass());tmpG.push_back(itH->return_gas());tmpD.push_back(itH->return_stellar_disk());tmpB.push_back(itH->return_stellar_bulge());tmpRB.push_back(itH->return_Rhl_bulge());tmpRD.push_back(itH->return_Rhl_disk());
//	cout <<  "ok4\n";			
		itH++;
		}
		cout <<"ok!\n";
		stack_Shalos->push(satillite_halo(DM,Z,interpol(tmpDM,newDM,tmpB).back(),interpol(tmpDM,newDM,tmpD).back(),interpol(tmpDM,newDM,tmpG).back()));
		stack_Shalos->top().put_Rhl_bulge(interpol(tmpDM,newDM,tmpRB).back()); stack_Shalos->top().put_Rhl_disk(interpol(tmpDM,newDM,tmpRD).back());

		cout << "we chose one!!  "; stack_Shalos->top().output_vars();
	}
	else if(DM-5<itZ->front().return_mass()){cout << "LARGER  " << DM -itZ->front().return_mass() <<endl; stack_Shalos->push(itZ->front());stack_Shalos->top().put_DM(DM);stack_Shalos->top().put_Z(Z);}
	else{cout<<"AHHHHHHHHH\n";itH=itZ->begin();while(itH!=itZ->end()){itH->output_vars();itH++;} exit(EXIT_FAILURE);}
	}}}
	
		
	if(!foundhalo || onlydisk){
	stars=Mdm_to_Mstellar(DM,Z);
	stack_Shalos->push(satillite_halo(DM,Z,-300,stars,get_Mgas(stars,Z,DM+log10(omega_baryon/omega_m))));
	stack_Shalos->top().put_Rhl_disk(get_Rdisk(stars,Z));
	stack_Shalos->top().put_Rhl_bulge(-300);
	
	}


/*	
	if(oblybulge){
	stars=Mdm_to_Mstellar(DM,Z);
	stack_Shalos->push(satillite_halo(DM,Z,stars,-300,get_Mgas(stars,Z,DM+log10(omega_baryon/omega_m))));
	stack_Shalos->top().put_Rhl_disk(get_Rdisk(stars,Z));
	}

*/

	if(includeshutdown){if(get_M_halo_shutdown(Z)<DM){stack_Shalos->top().put_shutdown();}}

	cout << "added  ";
	stack_Shalos->top().output_vars();
}



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



satillite_halo merge(satillite_halo central, satillite_halo satillite,double Z){

double DM=log10(pow(10,central.return_mass())+pow(10,satillite.return_mass()));
double stellar_bulge=central.return_stellar_bulge(),stellar_disk=central.return_stellar_disk(),gas=central.return_gas();
double mu=satillite.return_mu();
double r_bulge=-300,r_disk=-300;
bool major_merger=false;

if(includemerger){
	if(mu<0.3){
		
		// //minor merger//
		cout << "minor  "<<mu << "  " <<endl;
		satillite.output_vars();
		central.output_vars();
		gas=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()));
		stellar_bulge=log10(pow(10,central.return_stellar_bulge())+pow(10,satillite.return_stellar_bulge())+pow(10,satillite.return_stellar_disk()));
		stellar_disk=central.return_stellar_disk();
		double newstars=starburst(central,satillite);gas=log10(pow(10,gas)-newstars);stellar_disk=log10(pow(10,stellar_disk)+newstars);
		r_bulge=get_Rbulge(central.return_stellar_bulge(), central.return_Rhl_bulge(), satillite.return_stellar_total(),satillite.return_Rhl(),forb/0.5,gas);
		r_disk=get_Rdisk(log10(pow(10,stellar_disk)+pow(10,stellar_bulge)),Z);
	}
	else{
		// //major merger//
		major_merger=true;
		cout << "MAGOR!!!!!  "<<mu << "  " <<satillite.return_redshift() << "  " << satillite.return_mass() << "  " << central.return_mass() << "  " << satillite.return_stellar_total() << "  " << central.return_stellar_total() << "  " << satillite.return_gas() << "  " << central.return_gas() << "  " << satillite.return_Rhl() << "  " << central.return_Rhl() <<endl;
		cout<< "major  " << central.return_mass() << "  " << central.return_baryons() << "  " << satillite.return_mass() << "  " << satillite.return_baryons() << "  " << satillite.return_redshift() << endl;
		double newstars=starburst(central,satillite);
		cout << "gas  " << gas << endl;
		if(newstars>=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))){newstars=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()));gas=-300;}
		else{gas=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())-pow(10,newstars));}
		stellar_bulge=log10(pow(10,central.return_baryons())+pow(10,satillite.return_baryons())-pow(10,gas));
		stellar_disk=-300;
		r_bulge=get_Rbulge(central.return_baryons(),central.return_Rhl(),satillite.return_baryons(),satillite.return_Rhl(),forb/0.5,gas);  // // need to think of a better way of doing this!!!
		r_disk=-300; // // and this line!!
		cout << "gas  " << gas << "  " << newstars  << endl;
	}
}
satillite_halo new_cent(DM,Z,stellar_bulge,stellar_disk,gas);
new_cent.put_Rhl_bulge(r_bulge);
new_cent.put_Rhl_disk(r_disk);
if(satillite.return_extrememerger() || central.return_extrememerger()){new_cent.put_extrememerger();}
if(central.return_shutdown()){new_cent.put_shutdown();}
if(central.return_cgt4()){new_cent.put_cgt4();}
if(major_merger || central.return_had_major_merger()){new_cent.put_had_major_merger();}
cout << "new cent  "; new_cent.output_vars();
return new_cent;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



//mergeing
std::stack<satillite_halo> merge_all_galaxies(std::stack<satillite_halo> *stack_Shalo, std::stack<accreted_mass> *stack_acc,std::stack<accreted_mass> *stack_mpro){

std::stack<satillite_halo> merged_halos;
satillite_halo cent(-300,-300,-300,-300,-300);
double top_redshift,next_redshift, MDM=0.0, Mstellar=0.0,T;
std::list<satillite_halo> to_be_merged;

cent=stack_Shalo->top();
MDM=cent.return_mass();
Mstellar=cent.return_stellar_total();
stack_Shalo->pop();

while(! stack_acc->empty()){////

//cout << "SDFDSF  " << stack_Shalo->top().return_mass() <<endl;// << "  "  << stack_Shalo->top().return_redshift() << "  " << stack_acc->top().return_mass() << "   "<< stack_acc->top().return_redshift()<<endl;

	top_redshift=(*stack_acc).top().return_redshift();
	accreted_mass tmp(stack_acc->top().return_mass(),stack_acc->top().return_redshift(),stack_acc->top().return_addmass());
	stack_acc->pop();
	if(! stack_acc->empty()){next_redshift=(*stack_acc).top().return_redshift();}
	else{next_redshift=0;}
	stack_acc->push(tmp);

	T=get_time_from_redshift(next_redshift); //check this

	if(!stack_Shalo->empty()){
		while(! stack_Shalo->empty() && stack_Shalo->top().return_redshift()>=next_redshift){
			to_be_merged.push_front(stack_Shalo->top());
			if(MDM<to_be_merged.begin()->return_mass()-1.){to_be_merged.begin()->put_extrememerger();}
			to_be_merged.begin()->put_Tdy(T+get_T_dynamical_friction(MDM, to_be_merged.begin()->return_mass(),top_redshift));
			to_be_merged.begin()->put_mu(pow(10,(to_be_merged.begin()->return_baryons()-cent.return_baryons())));
			cout << "mudfgdfg  " << to_be_merged.begin()->return_mu() << "  " << to_be_merged.begin()->return_baryons() << "  " <<cent.return_baryons() <<endl;
			stack_Shalo->pop();
		}
	}


	//cycle through list of to be merged//
	std::list<satillite_halo>::iterator it = to_be_merged.begin();
	while (it!= to_be_merged.end()){
		if(it->return_Tdy()<T){
			cent=merge(cent,*it,top_redshift);
			MDM=log10(pow(10.,MDM)+pow(10.,it->return_mass()));
			//cout << "M stellar  " << it->return_baryons() << endl;
			Mstellar=log10(pow(10.,Mstellar)+pow(10.,it->return_stellar_total()));
			it = to_be_merged.erase(it);
		}
		else{
			if(! it->return_shutdown()){
				if(includeshutdown && get_M_halo_shutdown(next_redshift)<it->return_mass()){it->put_shutdown();}
			else{get_ssfr(*it,top_redshift,next_redshift);}
			}
			++it;
		}
	}
	cout << "MDMtest  " << MDM << "  " << stack_acc->top().return_mass() << "  " << stack_acc->top().return_addmass() <<endl;

	if(reinitialiseDM){MDM=stack_mpro->top().return_mass();cent.put_DM(MDM);}
	else{
		if(stack_acc->top().return_addmass()){cent.add_DM(stack_acc->top().return_mass());MDM=log10(pow(10.0,MDM)+pow(10.0,stack_acc->top().return_mass()));}
		else{cent.subtract_DM(stack_acc->top().return_mass());MDM=log10(pow(10.0,MDM)-pow(10.0,stack_acc->top().return_mass()));}
	}

	cent.put_Z(top_redshift);
	cout << "MDMtest2  " << MDM << "  " << stack_acc->top().return_mass() << "  " << stack_acc->top().return_addmass() <<endl;
	if(includeshutdown){
		if(! cent.return_shutdown()){
			if(get_M_halo_shutdown(top_redshift)<cent.return_mass()){
				cout <<"shutdown!!\n";cent.put_shutdown();
			}
		}
	}

	if(! cent.return_had_major_merger() && ! cent.return_shutdown()){
		cout << "reinitialise\n";
		if(cent.return_stellar_total() < Mdm_to_Mstellar(cent.return_mass(),top_redshift)){cent.add_stellar_disk(log10(pow(10,Mdm_to_Mstellar(cent.return_mass(),top_redshift))-pow(10,cent.return_stellar_total())));}
		if(cent.return_gas() < get_Mgas(cent.return_stellar_total(),top_redshift,cent.return_mass()+log10(omega_baryon/omega_m))){cent.put_gas(get_Mgas(cent.return_stellar_total(),top_redshift,cent.return_mass()+log10(omega_baryon/omega_m)));}
		if(includeclumpy_accretion){cout << "clumpy_before  " << cent.return_Rhl_bulge()<<endl;cent=clumpyaccretion(cent,top_redshift,next_redshift,top_redshift,0.);cout << "clumpy_after  " << cent.return_Rhl_bulge()<<endl;}
		if(includebathtub_model_accretion){bathtub_model_accretion(cent,top_redshift,next_redshift);}
	}
	
	if(includeexpansion){if(! cent.return_cgt4()){if(get_concentration(cent.return_mass(), top_redshift)>4){cent.put_cgt4();cent.put_Rhl_bulge(2*cent.return_Rhl_bulge());cout << "EXPANSION  " << top_redshift << "  " << cent.return_mass() << "  "  << cent.return_baryons() << endl;}}} //include mass loss?!

	//double Mdisk=log10(pow(10,cent.return_gas())+pow(10,cent.return_stellar_disk()));
	double rdisk=get_Rdisk(cent.return_stellar_total(),top_redshift);
	cent.put_Rhl_disk(rdisk);
	cout << "merge_tree  " << cent.return_shutdown()  << "  ";
	cent.output_vars();
	merged_halos.push(cent);
	stack_acc->pop();
	stack_mpro->pop();
} /////////


cout << "end  " << stack_Shalo->empty() << "  " <<  to_be_merged.empty() << "\n";

cout << "central ";
cent.output_vars();
//cout << "ok\n";
if(! stack_Shalo->empty()){while (! stack_Shalo->empty()){cout << "satillite  ";stack_Shalo->top().output_vars();stack_Shalo->pop();}}
if(! stack_acc->empty()){cout << "acc "; while(! stack_acc->empty()){cout << stack_acc->top().return_mass() << "  ";stack_acc->pop();}cout << endl;}
if(! to_be_merged.empty()){while(! to_be_merged.empty()){cout<<"merge   ";to_be_merged.front().output_vars();to_be_merged.pop_front();}} 

cout << "length  " << merged_halos.size() <<endl;
return merged_halos;
} 



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


 
// main program



int main(int argc,char *argv[]) {
double fileDM,fileZ,stars;
double endtreeDM,endtreeZ;

int i;

cout << "cosmology   " << omega_m << "  " << h << "  " << omega_baryon << "  " << sigma_8 << "  " <<endl;
cout << "mergers?  " << includemerger << endl;
cout << "halo shutdown? " << includeshutdown <<endl;
cout << "disk regrowth?  " << includediskregrowth << endl;
cout << "clumpy accretion?  " << includeclumpy_accretion <<endl;
cout << "bathtub model accretion?  " << includebathtub_model_accretion << endl;
cout << "gas dissipation?  " << includegasdissipation << endl;
cout << "forb  " << forb << endl;

listgals emp;
listlistgals mergercat;
while(i<=6000){mergercat.push_back(emp);i++;}
cout <<"mergetcat!!!  " << mergercat.size()<<endl;
listlistgals::iterator itZ;
listgals::iterator itH;
std::stack<satillite_halo> tmpmergedstack;

//stack of satilite halos
std::stack<satillite_halo> stack_Shalos;
//stack of accreted DM
std::stack<accreted_mass> stack_acc;
//stack of mass of the main progenitor
std::stack<accreted_mass> stack_mpro;


    // get merger tree filename fron argument //
    // decode arguments
    if(argc <2) {
        printf("You must provide at least one argument\n");
        exit(0);
    }
   // report settings
   string filename=argv[1];
   cout << filename << endl;
int k=-1;
//while (k<600){
cout<<"RESET\n";
k+=1;
ifstream myfile;
myfile.open (filename);
string line;
int j;
i=0;
string tmp;



while ( getline (myfile,line))
    {
	std::istringstream ss(line);
        std::istream_iterator<std::string> begin(ss), end;
        std::vector<std::string> words(begin, end); 

	if(words[0]=="START"){i+=1;}
	else if(words[0]=="mpro"){
		stack_mpro.push(accreted_mass((double)atof(words[1].c_str()),(double)atof(words[2].c_str()),1));
		j=0;}
	else if(words[0]=="1"){
		if(words[1]!="acc"){
			if(j==0){endtreeDM=(double)atof(words[2].c_str());endtreeZ=(double)atof(words[1].c_str());}
			else{
				fileDM=(double)atof(words[2].c_str());
				fileZ=(double)atof(words[1].c_str());
				find_galaxy(fileDM,fileZ,&mergercat,&stack_Shalos);
			}
		j+=1;
		}
		else if(words[1]=="acc"){
			if(reinitialiseDM){stack_acc.push(accreted_mass(-300,(double)atof(words[2].c_str()),1));}
			else{stack_acc.push(accreted_mass((double)atof(words[3].c_str()),(double)atof(words[2].c_str()),(bool)atof(words[4].c_str())));}
		}
	}


	if(words[0]=="START") {
		find_galaxy(endtreeDM,endtreeZ,&mergercat,&stack_Shalos);
		
		cout << "NEW\n";
		tmpmergedstack=merge_all_galaxies(&stack_Shalos,&stack_acc,&stack_mpro);
		cout <<"DONE\n";
		
		while(!tmpmergedstack.empty()){
			if(tmpmergedstack.top().return_mass()==tmpmergedstack.top().return_mass()){
				itZ=mergercat.begin();
				std::advance(itZ, (int)(tmpmergedstack.top().return_redshift()*100.));
				itH=itZ->begin();
				while(itH->return_mass()>tmpmergedstack.top().return_mass()){itH++;}
				itZ->insert(itH,tmpmergedstack.top());
			}
			tmpmergedstack.pop();
		}
		while(!tmpmergedstack.empty()){cout<<tmpmergedstack.top().return_redshift()<<endl;tmpmergedstack.pop();}
		cout << "done merging"<< endl;

		while(! stack_Shalos.empty()){cout << "Shalo sdaf  " << stack_Shalos.top().return_mass()<<"  " << stack_Shalos.top().return_redshift() <<endl; stack_Shalos.pop();}
		while(! stack_acc.empty()){cout << "acc sdaf " << stack_acc.top().return_mass()<<"  " << stack_acc.top().return_redshift() <<endl; stack_acc.pop();}
	}
}
myfile.close();

return(0);

}
