//Activated Region Fitting C-source.
//Wouter D. Weeda.
//Libary of functions for ARF


#include<R.h>
#include<Rmath.h>
#include<R_ext/Utils.h>


void gauss(double *theta, int *np, int *dimx, int *dimy, int *dimz, double *gx)
{

	int reg,x,y,z,p;
	double f,theta_x,theta_y,theta_z,sig_x,sig_y,sig_z,sig_xy,sig_xz,sig_yz,det_sig,dif_x,dif_y,dif_z;


	//theta 1,2,3 = x,y,z coordinates
	//theta 4,5,6 = sd's of x,y,z
	//theta 7,8,9 = corr xy, xz, yz
	//theta 10    = amplitude

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				f=0; //f becomes the sum of all regions  (zeroed every region)

				for(reg=0;reg<(*np);reg=reg+10) {
					//parameter coordinates
					theta_x=theta[reg+0];
					theta_y=theta[reg+1];
					theta_z=theta[reg+2];

					//sigma matrix
					sig_x=pow(theta[reg+3],2);
					sig_xy=theta[reg+6]*theta[reg+4]*theta[reg+3];
					sig_xz=theta[reg+7]*theta[reg+5]*theta[reg+3];
					sig_y=pow(theta[reg+4],2);
					sig_yz=theta[reg+8]*theta[reg+4]*theta[reg+5];
					sig_z=pow(theta[reg+5],2);


					//determinant of sigma
					det_sig=sig_x*sig_y*sig_z-sig_x*sig_yz*sig_yz-sig_xy*sig_xy*sig_z+sig_xy*sig_xz*sig_yz+sig_xz*sig_xy*sig_yz-sig_xz*sig_y*sig_xz;
					if(det_sig < 0) det_sig=0;

					//(x-pc)
					dif_x=(x-theta_x);
					dif_y=(y-theta_y);
					dif_z=(z-theta_z);

					//add to f gaussian value for each region
					f=f+theta[reg+9]*(1/(pow(sqrt(2*M_PI),3)*sqrt(det_sig)))*exp(-.5*(dif_x*(dif_x*(sig_y*sig_z-sig_yz*sig_yz)+dif_y*(sig_yz*sig_xz-sig_xy*sig_z)+dif_z*(sig_xy*sig_yz-sig_y*sig_xz))/det_sig+dif_y*(dif_x*(sig_xz*sig_yz-sig_z*sig_xy)+dif_y*(sig_x*sig_z-sig_xz*sig_xz)+dif_z*(sig_xy*sig_xz-sig_x*sig_yz))/det_sig+dif_z*(dif_x*(sig_xy*sig_yz-sig_xz*sig_y)+dif_y*(sig_xz*sig_xy-sig_yz*sig_x)+dif_z*(sig_x*sig_y-sig_xy*sig_xy))/det_sig));

				}

				gx[p]=f; //set output vector to sum of gaussians
				p++;

			}
		}
	}

	void R_CheckUserInterrupt(void);

}

void ssqgauss(double *theta, double *dat, double *W, int *brain, int *np, int *dimx, int *dimy, int *dimz, double *ss)
{

	int reg,x,y,z,p;
	double f,g,theta_x,theta_y,theta_z,sig_x,sig_y,sig_z,sig_xy,sig_xz,sig_yz,det_sig,dif_x,dif_y,dif_z;

	//theta 1,2,3 = x,y,z coordinates
	//theta 4,5,6 = sd's of x,y,z
	//theta 7,8,9 = corr xy, xz, yz
	//theta 10    = amplitude

	p=0;
	g=0e0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				f=0; //f becomes the sum of all regions  (zeroed every region)

				if(brain[p]!=0) {

					for(reg=0;reg<(*np);reg=reg+10) {

						//parameter coordinates
						theta_x=theta[reg+0];
						theta_y=theta[reg+1];
						theta_z=theta[reg+2];

						//sigma matrix
						sig_x=pow(theta[reg+3],2);
						sig_xy=theta[reg+6]*theta[reg+4]*theta[reg+3];
						sig_xz=theta[reg+7]*theta[reg+5]*theta[reg+3];
						sig_y=pow(theta[reg+4],2);
						sig_yz=theta[reg+8]*theta[reg+4]*theta[reg+5];
						sig_z=pow(theta[reg+5],2);


						//determinant of sigma
						det_sig=sig_x*sig_y*sig_z-sig_x*sig_yz*sig_yz-sig_xy*sig_xy*sig_z+sig_xy*sig_xz*sig_yz+sig_xz*sig_xy*sig_yz-sig_xz*sig_y*sig_xz;
						if(det_sig < 0) det_sig=0;

						//(x-pc)
						dif_x=(x-theta_x);
						dif_y=(y-theta_y);
						dif_z=(z-theta_z);

						//add to f gaussian value for each region
						f=f+theta[reg+9]*(1/(pow(sqrt(2*M_PI),3)*sqrt(det_sig)))*exp(-.5*(dif_x*(dif_x*(sig_y*sig_z-sig_yz*sig_yz)+dif_y*(sig_yz*sig_xz-sig_xy*sig_z)+dif_z*(sig_xy*sig_yz-sig_y*sig_xz))/det_sig+dif_y*(dif_x*(sig_xz*sig_yz-sig_z*sig_xy)+dif_y*(sig_x*sig_z-sig_xz*sig_xz)+dif_z*(sig_xy*sig_xz-sig_x*sig_yz))/det_sig+dif_z*(dif_x*(sig_xy*sig_yz-sig_xz*sig_y)+dif_y*(sig_xz*sig_xy-sig_yz*sig_x)+dif_z*(sig_x*sig_y-sig_xy*sig_xy))/det_sig));

					}
				}

				//sum (data-model)^2 over voxels and weight
				g=g+pow((dat[p]-f),2)*(1/W[p]);
				p++;

			}
		}
	}

	void R_CheckUserInterrupt(void);

	//set ss to g
	ss[0]=g;
}

void ssqdata(double *dat, double *W, int *brain, int *n, double *ss)
{

	int i;
	double g;

	g=0e0;
	for(i=0;i<(*n);i++) {
		if(brain[i]!=0) {
			g=g+pow((dat[i]),2)*(1/W[i]);
		}
	}

	ss[0]=g;
}

void simplegauss(double *theta, int *np, int *dimx, int *dimy, int *dimz, double *gx)
{

	int reg,x,y,z,p;
	double f,theta_x,theta_y,theta_z,sig_x,det_sig,dif_x,dif_y,dif_z;

	//theta 1,2,3 = x,y,z coordinates
	//theta 4,5,6 = sd's of x,y,z
	//theta 7,8,9 = corr xy, xz, yz
	//theta 10    = amplitude

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				f=0e0; //f becomes the sum of all regions  (zeroed every region)

				for(reg=0;reg<(*np);reg=reg+5) {
					//parameter coordinates
					theta_x=theta[reg+0];
					theta_y=theta[reg+1];
					theta_z=theta[reg+2];

					//sigma matrix
					sig_x=pow(theta[reg+3],2);

					//determinant of sigma
					det_sig=pow(sig_x,3);
					if(det_sig < 0) det_sig=0;

					//(x-pc)
					dif_x=pow((x-theta_x),2);
					dif_y=pow((y-theta_y),2);
					dif_z=pow((z-theta_z),2);

					//add to f gaussian value for each region
					f=f+theta[reg+4]*(1/(pow(sqrt(2*M_PI),3)*sqrt(det_sig)))*exp(-.5*(dif_x/sig_x+dif_y/sig_x+dif_z/sig_x));

				}

				gx[p]=f; //set output vector to sum of gaussian
				p++;

			}
		}
	}

	void R_CheckUserInterrupt(void);
}


void simplessqgauss(double *theta, double *dat, double *W, int *brain, int *np, int *dimx, int *dimy, int *dimz, double *ss)
{

	int reg,x,y,z,p;
	double f,g,theta_x,theta_y,theta_z,sig_x,det_sig,dif_x,dif_y,dif_z;

	//theta 1,2,3 = x,y,z coordinates
	//theta 4,5,6 = sd's of x,y,z
	//theta 7,8,9 = corr xy, xz, yz
	//theta 10    = amplitude

	p=0;
	g=0e0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				f=0e0; //f becomes the sum of all regions  (zeroed every region)

				if(brain[p]!=0) {

					for(reg=0;reg<(*np);reg=reg+5) {
						//parameter coordinates
						theta_x=theta[reg+0];
						theta_y=theta[reg+1];
						theta_z=theta[reg+2];

						//sigma matrix
						sig_x=pow(theta[reg+3],2);

						//determinant of sigma
						det_sig=pow(sig_x,3);
						if(det_sig < 0) det_sig=0;

						//(x-pc)
						dif_x=pow((x-theta_x),2);
						dif_y=pow((y-theta_y),2);
						dif_z=pow((z-theta_z),2);

						//add to f gaussian value for each region
						f=f+theta[reg+4]*(1/(pow(sqrt(2*M_PI),3)*sqrt(det_sig)))*exp(-.5*(dif_x/sig_x+dif_y/sig_x+dif_z/sig_x));
					}
				}

				//sum (data-model)^2 over voxels and weight
				g=g+pow((dat[p]-f),2)*(1/W[p]);
				p++;

			}
		}
	}

	void R_CheckUserInterrupt(void);

	//set ss to g
	ss[0]=g;
}


void innerSWdiag(int *n, int *p, int *runs, char **fnderiv, char **fnresid, char **fnweight, double *B)
{

	int i,j,k,l,m, Brow, Bcol;
	FILE *fderiv, *fresid, *fweight;
	double *Fv,*Ft, *Rv, *Wv, *Mr, s, Bm[*p][*p];

	Fv = (double *) R_alloc(*n,sizeof(double));
	Ft = (double *) R_alloc(*n,sizeof(double));
	Rv = (double *) R_alloc(*n,sizeof(double));
 	Wv = (double *) R_alloc(*n,sizeof(double));
 	Mr = (double *) R_alloc(*n,sizeof(double));

 	fderiv=fopen(*fnderiv,"r"); //NxP derivs vector (n incr. fastest)
 	fresid=fopen(*fnresid,"r"); //Nxrun residual vector (n incr fastest)
 	fweight=fopen(*fnweight,"r"); //N vector of weights

 	fread(Wv,sizeof(double),*n,fweight);
 	fclose(fweight);

 	for(m=0;m<(*n);m++) {
 		*(Mr+m)=0e0;
 	}

	for(l=0;l<(*runs);l++) { // make mean residual vectors

		fseek(fresid,sizeof(double)*((l)**n),SEEK_SET);
		fread(Rv,sizeof(double),*n,fresid);

		for(m=0;m<(*n);m++) {
			*(Mr+m)=*(Mr+m)+((1/pow((double) *runs,2))**(Rv+m)**(Rv+m));
		}
	}

	for(m=0;m<(*n);m++) { //weight mean residuals
		*(Mr+m)=*(Mr+m)/(*(Wv+m)**(Wv+m));
	}


	fclose(fresid);

 	// OFF DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p-1);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop

			*(Ft+j)=*(Fv+j)**(Mr+j);

		} // end R loop


		for(Bcol=(Brow+1);Bcol<(*p);Bcol++) { //BVEC COL LOOP
			fseek(fderiv,sizeof(double)*(Bcol**n),SEEK_SET);
			fread(Fv,sizeof(double),*n,fderiv);

			s=0;
			for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
				s=s+Fv[i]*Ft[i];

			}

			Bm[Brow][Bcol]=s;
			Bm[Bcol][Brow]=s;

		}
	}

	//DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop

			*(Ft+j)=*(Fv+j)**(Mr+j);

		} // end R loop

		s=0;
		for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
			s=s+Fv[i]*Ft[i];

		}

		Bm[Brow][Brow]=s;

	}


	k=0;
	for(Bcol=0;Bcol<(*p);Bcol++) {
		for(Brow=0;Brow<(*p);Brow++) {
			*(B+k)=Bm[Brow][Bcol];
			k++;
		}
	}

	fclose(fderiv);

}



void innerSWfull(int *n, int *p, int *runs, char **fnderiv, char **fnresid, char **fnweight, char **fnmeanresid, double *B)
{

	int i,j,k,l, Brow, Bcol;
	FILE *fderiv, *fresid, *fweight, *fmeanres;
	double *Fv,*Ft, *Rv, *Wv, *Mr, s, Bm[*p][*p], *mrv;

	Fv = (double *) R_alloc(*n,sizeof(double));
	Ft = (double *) R_alloc(*n,sizeof(double));
	Rv = (double *) R_alloc(*runs**n,sizeof(double));
 	Wv = (double *) R_alloc(*n,sizeof(double));
 	Mr = (double *) R_alloc((*n),sizeof(double));
 	mrv = (double *) R_alloc((int) 1,sizeof(double));

 	fderiv=fopen(*fnderiv,"r"); //NxP derivs vector (n incr. fastest)
 	fresid=fopen(*fnresid,"r"); //Nxrun residual vector (n incr fastest)
 	fweight=fopen(*fnweight,"r"); //N vector of weights

 	fread(Wv,sizeof(double),*n,fweight);
 	fclose(fweight);

 	fread(Rv,sizeof(double),*runs**n,fresid);
 	fclose(fresid);

 	//make mean residuals
	fmeanres = fopen(*fnmeanresid,"w");
	for(i=0;i<(*n);i++) {
		for(j=0;j<(*n);j++) {
			mrv[0]=0e0;
			for(l=0;l<(*runs);l++) { //run loop
				mrv[0] = mrv[0] + ((1/pow((double) *runs,2))**(Rv+i+l**n)**(Rv+j+l**n));
			} //end run loop
			fwrite(mrv,sizeof(double),1,fmeanres);
		} //end row loop
	} // end column loop

 	void R_CheckUserInterrupt(void);

 	fclose(fmeanres);
 	fmeanres = fopen(*fnmeanresid,"r");

 	// OFF DIAGONAL B, LOOP
 	for(Brow=0;Brow<(*p-1);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			fseek(fmeanres,sizeof(double)*(j**n),SEEK_SET);
			fread(Mr,sizeof(double),*n,fmeanres);

			s=0e0;
			for(i=0;i<(*n);i++) {  // R matrix row loop
				s=s+*(Fv+i)*(*(Mr+i)/(*(Wv+j)**(Wv+i)));
			}

			*(Ft+j)=s;

		} // end R loop

		for(Bcol=(Brow+1);Bcol<(*p);Bcol++) { //BVEC COL LOOP
			fseek(fderiv,sizeof(double)*(Bcol**n),SEEK_SET);
			fread(Fv,sizeof(double),*n,fderiv);

			s=0e0;
			for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
				s=s+Fv[i]*Ft[i];

			}

			Bm[Brow][Bcol]=s;
			Bm[Bcol][Brow]=s;
			void R_CheckUserInterrupt(void);
			//Rprintf("Done %d %d\n",Brow,Bcol);

		}

	}

 	fclose(fmeanres);
 	fmeanres = fopen(*fnmeanresid,"r");

 	//DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			fseek(fmeanres,sizeof(double)*(j**n),SEEK_SET);
			fread(Mr,sizeof(double),*n,fmeanres);

			s=0e0;
			for(i=0;i<(*n);i++) {  // R matrix row loop
				s=s+*(Fv+i)*(*(Mr+i)/(*(Wv+j)**(Wv+i)));
			}

			*(Ft+j)=s;

		} // end R loop

		s=0e0;
		for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
			s=s+Fv[i]*Ft[i];
		}

		Bm[Brow][Brow]=s;
		void R_CheckUserInterrupt(void);
		//Rprintf("Done %d %d\n",Brow,Brow);
	}

	k=0;
	for(Bcol=0;Bcol<(*p);Bcol++) {
		for(Brow=0;Brow<(*p);Brow++) {
			*(B+k)=Bm[Brow][Bcol];
			k++;
		}
	}


 	fclose(fderiv);

}

void innerSWfast(int *n, int *p, int *runs, char **fnderiv, char **fnresid, char **fnweight, double *B)
{

	int i,j,k,l,m, Brow, Bcol;
	FILE *fderiv, *fresid, *fweight;
	double *Fv,*Ft, *Rv, *Wv, *Mr, s, Bm[*p][*p];


	//Rprintf("band1: %d\n",*band);

	Fv = (double *) R_alloc(*n,sizeof(double));
	Ft = (double *) R_alloc(*n,sizeof(double));
	Rv = (double *) R_alloc(*n,sizeof(double));
 	Wv = (double *) R_alloc(*n,sizeof(double));
 	Mr = (double *) R_alloc((*n)*(*n),sizeof(double));

 	fderiv=fopen(*fnderiv,"r"); //NxP derivs vector (n incr. fastest)
 	fresid=fopen(*fnresid,"r"); //Nxrun residual vector (n incr fastest)
 	fweight=fopen(*fnweight,"r"); //N vector of weights

 	fread(Wv,sizeof(double),*n,fweight);
 	fclose(fweight);

	for(m=0;m<(*n);m++) {
		for(j=0;j<(*n);j++) {
			*(Mr+m+j**n)=0e0;
		}
 	}

 	//make mean residuals
	for(l=0;l<(*runs);l++) { //run loop

		fseek(fresid,sizeof(double)*((l)**n),SEEK_SET);
		fread(Rv,sizeof(double),*n,fresid);

		for(i=0;i<(*n);i++) {
			for(j=0;j<(*n);j++) {
				*(Mr+j+i**n)=*(Mr+j+i**n)+((1/pow((double) *runs,2))**(Rv+i)**(Rv+j));
			} //end column loop
		} //end row loop
 	} // end run loop

 	//Rprintf("Made mean residual band vectors\n");


 	// OFF DIAGONAL B, LOOP
 	for(Brow=0;Brow<(*p-1);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n
			s=0e0;
			for(i=0;i<(*n);i++) {  // R matrix row loop
				s=s+*(Fv+i)*(*(Mr+j+i**n)/(*(Wv+j)**(Wv+i)));
			}

			*(Ft+j)=s;

		} // end R loop

		for(Bcol=(Brow+1);Bcol<(*p);Bcol++) { //BVEC COL LOOP
			fseek(fderiv,sizeof(double)*(Bcol**n),SEEK_SET);
			fread(Fv,sizeof(double),*n,fderiv);

			s=0e0;
			for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
				s=s+Fv[i]*Ft[i];

			}

			Bm[Brow][Bcol]=s;
			Bm[Bcol][Brow]=s;
			void R_CheckUserInterrupt(void);
			//Rprintf("Done %d %d\n",Brow,Bcol);

		}

	}


	//DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n
			s=0e0;
			for(i=0;i<(*n);i++) {  // R matrix row loop
				s=s+*(Fv+i)*(*(Mr+j+i**n)/(*(Wv+j)**(Wv+i)));
			}

			*(Ft+j)=s;

		} // end R loop


		s=0e0;
		for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
			s=s+Fv[i]*Ft[i];
		}

		Bm[Brow][Brow]=s;
		void R_CheckUserInterrupt(void);
		//Rprintf("Done %d %d\n",Brow,Brow);
	}



	k=0;
	for(Bcol=0;Bcol<(*p);Bcol++) {
		for(Brow=0;Brow<(*p);Brow++) {
			*(B+k)=Bm[Brow][Bcol];
			k++;
		}
	}


 	fclose(fderiv);
	fclose(fresid);

}

void innerSWbw(int *n, int *p, int *runs, int *bw, int *escapevar, int *Lv, char **fnderiv, char **fnresid, char **fnweight, char **fnmeanresid, double *B)
{

	int i,j,k,l, Brow, Bcol,loc;
	FILE *fderiv, *fresid, *fweight, *fmeanres;
	double *Fv,*Ft, *Rv, *Wv, *Mr, s, Bm[*p][*p], *mrv;

	Fv = (double *) R_alloc(*n,sizeof(double));
	Ft = (double *) R_alloc(*n,sizeof(double));
	Rv = (double *) R_alloc(*runs**n,sizeof(double));
 	Wv = (double *) R_alloc(*n,sizeof(double));
 	Mr = (double *) R_alloc((*n),sizeof(double));
 	mrv = (double *) R_alloc((int) 1,sizeof(double));

 	fweight=fopen(*fnweight,"r"); //N vector of weights
 	fread(Wv,sizeof(double),*n,fweight);
 	fclose(fweight);

 	fresid=fopen(*fnresid,"r"); //Nxrun residual vector (n incr fastest)
 	fread(Rv,sizeof(double),*runs**n,fresid);
 	fclose(fresid);

  	fderiv=fopen(*fnderiv,"r"); //NxP derivs vector (n incr. fastest)

 	//make mean residuals
	fmeanres = fopen(*fnmeanresid,"w");
	for(i=0;i<(*n);i++) {
		for(j=0;j<(*n);j++) {
			mrv[0]=0e0;
			for(l=0;l<(*runs);l++) { //run loop
				mrv[0] = mrv[0] + ((1/pow((double) *runs,2))**(Rv+i+l**n)**(Rv+j+l**n));
			} //end run loop
			fwrite(mrv,sizeof(double),1,fmeanres);
		} //end row loop
	} // end column loop

 	//void R_CheckUserInterrupt(void);
 	//Rprintf("Finished MEANRES\n");

 	fclose(fmeanres);
 	fmeanres = fopen(*fnmeanresid,"r");

 	// OFF DIAGONAL B, LOOP
 	for(Brow=0;Brow<(*p-1);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			fseek(fmeanres,sizeof(double)*(j**n),SEEK_SET);
			fread(Mr,sizeof(double),*n,fmeanres);

			s=0e0;
			for(i=0;i<(*bw);i++) {  // R matrix row loop
				loc = Lv[j+i**n];
				//Rprintf("loc: %d\n",loc);
				if(loc!=*escapevar)	s=s+*(Fv+loc)*(*(Mr+loc)/(*(Wv+j)**(Wv+loc)));
			}

			*(Ft+j)=s;

		} // end R loop

		for(Bcol=(Brow+1);Bcol<(*p);Bcol++) { //BVEC COL LOOP
			fseek(fderiv,sizeof(double)*(Bcol**n),SEEK_SET);
			fread(Fv,sizeof(double),*n,fderiv);

			s=0e0;
			for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
				s=s+Fv[i]*Ft[i];

			}

			Bm[Brow][Bcol]=s;
			Bm[Bcol][Brow]=s;
			//void R_CheckUserInterrupt(void);
			//Rprintf("Done %d %d\n",Brow,Bcol);

		}

	}

 	fclose(fmeanres);
 	fmeanres = fopen(*fnmeanresid,"r");

 	//DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			fseek(fmeanres,sizeof(double)*(j**n),SEEK_SET);
			fread(Mr,sizeof(double),*n,fmeanres);

			s=0e0;
			for(i=0;i<(*bw);i++) {  // R matrix row loop
				loc = Lv[j+i**n];
				if(loc!=*escapevar)	s=s+*(Fv+loc)*(*(Mr+loc)/(*(Wv+j)**(Wv+loc)));
			}
			*(Ft+j)=s;

		} // end R loop

		s=0e0;
		for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
			s=s+Fv[i]*Ft[i];

		}

		Bm[Brow][Brow]=s;
		//void R_CheckUserInterrupt(void);
		//Rprintf("Done %d %d\n",Brow,Brow);
	}

	k=0;
	for(Bcol=0;Bcol<(*p);Bcol++) {
		for(Brow=0;Brow<(*p);Brow++) {
			*(B+k)=Bm[Brow][Bcol];
			k++;
		}
	}

	fclose(fmeanres);
 	fclose(fderiv);

}


void innerSWbwfast(int *n, int *p, int *runs, int *bw, int *escapevar, int *Lv, char **fnderiv, char **fnresid, char **fnweight, char **fnmeanresid, double *B)
{

	int i,j,k,l, Brow, Bcol,loc;
	FILE *fderiv, *fresid, *fweight; // *fmeanres;
	//double *Fv,*Ft, *Rv, *Wv, *Mr, s, Bm[*p][*p], *mrv;
	double *Fv,*Ft, *Rv, *Wv, s, Bm[*p][*p], *mrv;

	Fv = (double *) R_alloc(*n,sizeof(double));
	Ft = (double *) R_alloc(*n,sizeof(double));
	Rv = (double *) R_alloc(*runs**n,sizeof(double));
 	Wv = (double *) R_alloc(*n,sizeof(double));
 	//Mr = (double *) R_alloc((*n),sizeof(double));
 	mrv = (double *) R_alloc((*n)**bw,sizeof(double));
 	//mrv = (double *) R_alloc((int) 1,sizeof(double));

 	fweight=fopen(*fnweight,"r"); //N vector of weights
 	fread(Wv,sizeof(double),*n,fweight);
 	fclose(fweight);

 	fresid=fopen(*fnresid,"r"); //Nxrun residual vector (n incr fastest)
 	fread(Rv,sizeof(double),*runs**n,fresid);
 	fclose(fresid);

  	fderiv=fopen(*fnderiv,"r"); //NxP derivs vector (n incr. fastest)

 	//make mean residuals
	//fmeanres = fopen(*fnmeanresid,"w");
	for(i=0;i<(*n);i++) {
		for(j=0;j<(*bw);j++) {
			for(l=0;l<(*runs);l++) { //run loop
				loc = Lv[i+j**n];
				//Rprintf("loc%d, i%d, j%d, l%d ijn%d\n",loc,i,j,l,i+j**n);
				if(loc!=*escapevar)	*(mrv+i+j**n) = *(mrv+i+j**n) + ((1/pow((double) *runs,2))**(Rv+i+l**n)**(Rv+loc+l**n));
				else *(mrv+i+j**n)=0;
			} //end run loop

			//fwrite(mrv,sizeof(double),1,fmeanres);
		} //end row loop
	} // end column loop

 	//void R_CheckUserInterrupt(void);
 	//Rprintf("Finished MEANRES\n");

 	//fclose(fmeanres);
 	//fmeanres = fopen(*fnmeanresid,"r");

 	// OFF DIAGONAL B, LOOP
 	for(Brow=0;Brow<(*p-1);Brow++) { //BVEC ROW LOOP

 		//Rprintf("brow: %d\n",Brow);
 		//void R_CheckUserInterrupt(void);

 		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			//Rprintf("j%d\n",j);
			//void R_CheckUserInterrupt(void);
			s=0e0;
			for(i=0;i<(*bw);i++) {  // R matrix row loop
				loc = Lv[j+i**n];
				if(loc!=*escapevar)	s=s+*(Fv+loc)*(*(mrv+j+i**n)/(*(Wv+j)**(Wv+loc)));

			}

			*(Ft+j)=s;

		} // end R loop

		for(Bcol=(Brow+1);Bcol<(*p);Bcol++) { //BVEC COL LOOP
			fseek(fderiv,sizeof(double)*(Bcol**n),SEEK_SET);
			fread(Fv,sizeof(double),*n,fderiv);

			s=0e0;
			for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
				s=s+Fv[i]*Ft[i];

			}

			Bm[Brow][Bcol]=s;
			Bm[Bcol][Brow]=s;
			//void R_CheckUserInterrupt(void);
			//Rprintf("Done %d %d\n",Brow,Bcol);

		}

	}


 	//DIAGONAL B, LOOP
	for(Brow=0;Brow<(*p);Brow++) { //BVEC ROW LOOP

		fseek(fderiv,sizeof(double)*(Brow**n),SEEK_SET);
		fread(Fv,sizeof(double),*n,fderiv);

		for(j=0;j<(*n);j++) { // R matrix column loop ) j = 1..n

			s=0e0;
			for(i=0;i<(*bw);i++) {  // R matrix row loop
				loc = Lv[j+i**n];
				if(loc!=*escapevar)	s=s+*(Fv+loc)*(*(mrv+j+i**n)/(*(Wv+j)**(Wv+loc)));
			}
			*(Ft+j)=s;

		} // end R loop

		s=0e0;
		for(i=0;i<(*n);i++) { // (t(F)%*%R)%*%F LOOP
			s=s+Fv[i]*Ft[i];

		}

		Bm[Brow][Brow]=s;
		//void R_CheckUserInterrupt(void);
		//Rprintf("Done %d %d\n",Brow,Brow);
	}

	k=0;
	for(Bcol=0;Bcol<(*p);Bcol++) {
		for(Brow=0;Brow<(*p);Brow++) {
			*(B+k)=Bm[Brow][Bcol];
			k++;
		}
	}

 	fclose(fderiv);

}


void approxHessian(int *np, int *n, char **dffile, char **wfile, double *hessian) {

	int r,c,j,k;
	double hess[*np][*np], *cdf, *W, *rdf;
	FILE *fderiv, *fweight;

	cdf = (double *) R_alloc(*n,sizeof(double));
	rdf = (double *) R_alloc(*n,sizeof(double));
	W = (double *) R_alloc(*n,sizeof(double));

	fweight=fopen(*wfile,"r");
 	fread(W,sizeof(double),*n,fweight);
 	fclose(fweight);

	fderiv=fopen(*dffile,"r");

	//diag loop
	for(r=0;r<(*np);r++) {

		fseek(fderiv,sizeof(double)*(r**n),SEEK_SET);
		fread(rdf,sizeof(double),*n,fderiv);

		hess[r][r] = 0e0;

		for(j=0;j<(*n);j++) {
			hess[r][r]=hess[r][r]+(rdf[j]*(rdf[j]/W[j]));
		}

	}


	//off diag loop
	for(r=0;r<(*np-1);r++) {

		fseek(fderiv,sizeof(double)*(r**n),SEEK_SET);
		fread(rdf,sizeof(double),*n,fderiv);

		for(c=r+1;c<(*np);c++) {

			fseek(fderiv,sizeof(double)*(c**n),SEEK_SET);
			fread(cdf,sizeof(double),*n,fderiv);

			hess[r][c] = 0e0;
			hess[c][r] = 0e0;

			for(j=0;j<(*n);j++) {
				hess[r][c]=hess[r][c]+(rdf[j]*(cdf[j]/W[j]));
				hess[c][r]=hess[c][r]+(rdf[j]*(cdf[j]/W[j]));

			}

		}

	}


	k=0;
	for(c=0;c<(*np);c++) {
		for(r=0;r<(*np);r++) {
			*(hessian+k)=2*hess[r][c];
			k++;
		}
	}

	fclose(fderiv);

}

void dfgaussFile(int *np, int *brain, int *dimx, int *dimy, int *dimz, double *thetavec, char **filename)
{

	void dftheta0();
	void dftheta1();
	void dftheta2();
	void dftheta3();
	void dftheta4();
	void dftheta5();
	void dftheta6();
	void dftheta7();
	void dftheta8();
	void dftheta9();

	int i,n=(*dimx)*(*dimy)*(*dimz), reg,p;
	double *grad, *theta;
	grad = (double *) R_alloc((n),sizeof(double));
	theta = (double *) R_alloc((10),sizeof(double));

	p=0;
	for(i=0;i<n;i++) p=p + brain[i];

	FILE *f;
	f=fopen(*filename,"w");

	for(reg=0;reg<(*np);reg=reg+10) {

		for(i=0;i<10;i++) {
			*(theta+i)=*(thetavec+reg+i);
		}


		dftheta0(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta1(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta2(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta3(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta4(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta5(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta6(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta7(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta8(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);
		dftheta9(theta,brain,dimx,dimy,dimz,grad);
		fwrite(grad,sizeof(double),p,f);

	}

	fclose(f);
}

void dfgauss(int *np, int *brain, int *dimx, int *dimy, int *dimz, double *thetavec, double *derivs)
{

	void dftheta0();
	void dftheta1();
	void dftheta2();
	void dftheta3();
	void dftheta4();
	void dftheta5();
	void dftheta6();
	void dftheta7();
	void dftheta8();
	void dftheta9();

	int i,j,n=(*dimx)*(*dimy)*(*dimz), reg,p;
	double *theta,*grad;
	p=0;
	for(i=0;i<n;i++) p = p + brain[i];

	theta = (double *) R_alloc((10),sizeof(double));
	grad = (double *) R_alloc((p),sizeof(double));


	for(reg=0;reg<(*np);reg=reg+10) {

		for(i=0;i<10;i++) {
			*(theta+i)=*(thetavec+reg+i);
		}

		dftheta0(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(0*p)+(p*reg)]=grad[j];
		}

		dftheta1(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(1*p)+(p*reg)]=grad[j];
		}

		dftheta2(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(2*p)+(p*reg)]=grad[j];
		}

		dftheta3(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(3*p)+(p*reg)]=grad[j];
		}

		dftheta4(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(4*p)+(p*reg)]=grad[j];
		}

		dftheta5(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(5*p)+(p*reg)]=grad[j];
		}

		dftheta6(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(6*p)+(p*reg)]=grad[j];
		}

		dftheta7(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(7*p)+(p*reg)]=grad[j];
		}

		dftheta8(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(8*p)+(p*reg)]=grad[j];
		}

		dftheta9(theta,brain,dimx,dimy,dimz,grad);
		for(j=0;j<p;j++) {
			derivs[j+(9*p)+(p*reg)]=grad[j];
		}

	}

}


void dfssq(int *np, int *brain, int *dimx, int *dimy, int *dimz, double *thetavec, double *data, double *model, double *weights, double *ssqgrad)
{

	void dftheta0();
	void dftheta1();
	void dftheta2();
	void dftheta3();
	void dftheta4();
	void dftheta5();
	void dftheta6();
	void dftheta7();
	void dftheta8();
	void dftheta9();

	int i,j,n=(*dimx)*(*dimy)*(*dimz), reg,p;
	double *grad, *theta;
	grad = (double *) R_alloc((n),sizeof(double));
	theta = (double *) R_alloc((10),sizeof(double));

	for(reg=0;reg<(*np);reg=reg+10) {

		for(i=0;i<10;i++) {
			*(theta+i)=*(thetavec+reg+i);
		}

		dftheta0(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+0)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+0)=*(ssqgrad+reg+0)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+0)=*(ssqgrad+reg+0)*-2;

		dftheta1(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+1)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+1)=*(ssqgrad+reg+1)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+1)=*(ssqgrad+reg+1)*-2;

		dftheta2(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+2)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+2)=*(ssqgrad+reg+2)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+2)=*(ssqgrad+reg+2)*-2;

		dftheta3(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+3)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+3)=*(ssqgrad+reg+3)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+3)=*(ssqgrad+reg+3)*-2;

		dftheta4(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+4)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+4)=*(ssqgrad+reg+4)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+4)=*(ssqgrad+reg+4)*-2;

		dftheta5(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+5)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+5)=*(ssqgrad+reg+5)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+5)=*(ssqgrad+reg+5)*-2;

		dftheta6(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+6)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+6)=*(ssqgrad+reg+6)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+6)=*(ssqgrad+reg+6)*-2;

		dftheta7(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+7)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+7)=*(ssqgrad+reg+7)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+7)=*(ssqgrad+reg+7)*-2;

		dftheta8(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+8)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+8)=*(ssqgrad+reg+8)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+8)=*(ssqgrad+reg+8)*-2;

		dftheta9(theta,brain,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+9)=0e0;
		p=0;
		for(j=0;j<n;j=j+1) {
			if(brain[j]!=0) {
				*(ssqgrad+reg+9)=*(ssqgrad+reg+9)+((1/(*(weights+j)))**(grad+p)*(*(data+j)-(*(model+j))));
				p++;
			}
		}
		*(ssqgrad+reg+9)=*(ssqgrad+reg+9)*-2;

	}


}

void dfsimplegauss(int *np, int *dimx, int *dimy, int *dimz, double *thetavec, double *derivs)
{

	void dfsgtheta0();
	void dfsgtheta1();
	void dfsgtheta2();
	void dfsgtheta3();
	void dfsgtheta4();

	int i,j,n=(*dimx)*(*dimy)*(*dimz), reg;
	double *grad, *theta;
	grad = (double *) R_alloc((n),sizeof(double));
	theta = (double *) R_alloc((5),sizeof(double));


	for(reg=0;reg<(*np);reg=reg+5) {

		for(i=0;i<5;i++) {
			*(theta+i)=*(thetavec+reg+i);
		}

		dfsgtheta0(theta,dimx,dimy,dimz,grad);
		for(j=0;j<n;j++) derivs[j+(0*n)+(n*reg)]=grad[j];

		dfsgtheta1(theta,dimx,dimy,dimz,grad);
		for(j=0;j<n;j++) derivs[j+(1*n)+(n*reg)]=grad[j];

		dfsgtheta2(theta,dimx,dimy,dimz,grad);
		for(j=0;j<n;j++) derivs[j+(2*n)+(n*reg)]=grad[j];

		dfsgtheta3(theta,dimx,dimy,dimz,grad);
		for(j=0;j<n;j++) derivs[j+(3*n)+(n*reg)]=grad[j];

		dfsgtheta4(theta,dimx,dimy,dimz,grad);
		for(j=0;j<n;j++) derivs[j+(4*n)+(n*reg)]=grad[j];


	}

}

void dfsimplessq(int *np, int *dimx, int *dimy, int *dimz, double *thetavec, double *data, double *model, double *weights, double *ssqgrad)
{

	void dfsgtheta0();
	void dfsgtheta1();
	void dfsgtheta2();
	void dfsgtheta3();
	void dfsgtheta4();


	int i,j,n=(*dimx)*(*dimy)*(*dimz), reg;
	double *grad, *theta;
	grad = (double *) R_alloc((n),sizeof(double));
	theta = (double *) R_alloc((5),sizeof(double));


	for(reg=0;reg<(*np);reg=reg+5) {

		for(i=0;i<5;i++) {
			*(theta+i)=*(thetavec+reg+i);
		}

		dfsgtheta0(theta,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+0)=0e0;
		for(j=0;j<n;j=j+1) *(ssqgrad+reg+0)=*(ssqgrad+reg+0)+((1/(*(weights+j)))**(grad+j)*(*(data+j)-(*(model+j))));
		*(ssqgrad+reg+0)=*(ssqgrad+reg+0)*-2;

		dfsgtheta1(theta,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+1)=0e0;
		for(j=0;j<n;j=j+1) *(ssqgrad+reg+1)=*(ssqgrad+reg+1)+((1/(*(weights+j)))**(grad+j)*(*(data+j)-(*(model+j))));
		*(ssqgrad+reg+1)=*(ssqgrad+reg+1)*-2;

		dfsgtheta2(theta,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+2)=0e0;
		for(j=0;j<n;j=j+1) *(ssqgrad+reg+2)=*(ssqgrad+reg+2)+((1/(*(weights+j)))**(grad+j)*(*(data+j)-(*(model+j))));
		*(ssqgrad+reg+2)=*(ssqgrad+reg+2)*-2;

		dfsgtheta3(theta,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+3)=0e0;
		for(j=0;j<n;j=j+1) *(ssqgrad+reg+3)=*(ssqgrad+reg+3)+((1/(*(weights+j)))**(grad+j)*(*(data+j)-(*(model+j))));
		*(ssqgrad+reg+3)=*(ssqgrad+reg+3)*-2;

		dfsgtheta4(theta,dimx,dimy,dimz,grad);
		*(ssqgrad+reg+4)=0e0;
		for(j=0;j<n;j=j+1) *(ssqgrad+reg+4)=*(ssqgrad+reg+4)+((1/(*(weights+j)))**(grad+j)*(*(data+j)-(*(model+j))));
		*(ssqgrad+reg+4)=*(ssqgrad+reg+4)*-2;

	}


}



void dftheta0(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[33] = *(theta+9) * (1/(pow(sqrt(2 * M_PI),3) * sqrt(ev[29])));
					ev[35] = x - *(theta+0);
					ev[39] = ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[51] = ev[39] + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]);
					ev[56] = ev[20] * ev[10] - ev[7] * ev[15];
					ev[69] = ev[47] - ev[27];
					ev[82] = exp(-0.5 * (ev[35] * ev[51]/ev[29] + ev[40] * (ev[35] * ev[56] + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]))/ev[29] + ev[46] * (ev[35] * ev[69] + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]))/ev[29]));

					*(grad+p) = ev[33] * (ev[82] * (0.5 * (ev[46] * ev[69]/ev[29] + (ev[40] * ev[56]/ev[29] + (ev[39] + ev[51])/ev[29]))));

					p++;
				}
				a++;
			}
		}
	}

}

void dftheta1(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[33] = *(theta+9) * (1/(pow(sqrt(2 * M_PI),3) * sqrt(ev[29])));
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[43] = ev[10] * ev[20] - ev[15] * ev[7];
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[61] = ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]);
					ev[65] = ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[61] + ev[46] * (ev[21] - ev[11]);
					ev[72] = ev[24] - ev[10] * ev[4];
					ev[82] = exp(-0.5 * (ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * ev[43] + ev[46] * (ev[47] - ev[5] * ev[20]))/ev[29] + ev[40] * ev[65]/ev[29] + ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * ev[72] + ev[46] * (ev[6] - ev[16]))/ev[29]));

					*(grad+p) = ev[33] * (ev[82] * (0.5 * (ev[46] * ev[72]/ev[29] + ((ev[61] + ev[65])/ev[29] + ev[35] * ev[43]/ev[29]))));

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta2(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[33] = *(theta+9) * (1/(pow(sqrt(2 * M_PI),3) * sqrt(ev[29])));
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[49] = ev[47] - ev[5] * ev[20];
					ev[63] = ev[21] - ev[11];
					ev[76] = ev[46] * (ev[6] - ev[16]);
					ev[77] = ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[76];
					ev[82] = exp(-0.5 * (ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) +  ev[46] * ev[49])/ev[29] + ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * ev[63])/ev[29] + ev[46] * ev[77]/ev[29]));

					*(grad+p) = ev[33] * (ev[82] * (0.5 * ((ev[76] + ev[77])/ev[29] + (ev[40] * ev[63]/ev[29] + ev[35] * ev[49]/ev[29]))));

					p++;
				}
				a++;
			}
		}
	}
}


void dftheta3(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[14] = *(theta+6) * *(theta+4);
					ev[15] = ev[14] * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[19] = *(theta+7) * *(theta+5);
					ev[20] = ev[19] * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[88] = ev[14] * ev[10];
					ev[95] = 2 * *(theta+3);
					ev[96] = ev[95] * ev[5];
					ev[98] = ev[95] * ev[10];
					ev[103] = ev[14] * ev[15] + ev[15] * ev[14];
					ev[108] = ev[14] * ev[20] + ev[15] * ev[19];
					ev[113] = ev[19] * ev[15] + ev[20] * ev[14];
					ev[116] = ev[19] * ev[5];
					ev[120] = ev[96] * ev[7] - ev[98] * ev[10] - ev[103] * ev[7] + ev[108] * ev[10] + ev[113] * ev[10] - (ev[116] * ev[20] + ev[27] * ev[19]);
					ev[122] = pow(ev[29],2);

					*(grad+p)= -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[40] * (ev[10] * ev[19] - ev[14] * ev[7]) + ev[46] * (ev[88] - ev[5] * ev[19]))/ev[29] - ev[52] * ev[120]/ev[122] + (ev[40] * (ev[35] * (ev[19] * ev[10] - ev[7] * ev[14]) + ev[40] * (ev[95] * ev[7] - (ev[19] * ev[20] + ev[20] * ev[19])) + ev[46] * (ev[108] - ev[98]))/ev[29] - ev[66] * ev[120]/ev[122]) + (ev[46] * (ev[35] * (ev[88] - ev[116]) + ev[40] * (ev[113] - ev[10] * ev[95]) + ev[46] * (ev[96] - ev[103]))/ev[29] - ev[78] * ev[120]/ev[122])))) + *(theta+9) * (ev[3] * (0.5 * (ev[120] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta4(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[84] = 2 * *(theta+4);
					ev[86] = *(theta+8) * *(theta+5);
					ev[93] = *(theta+6) * *(theta+3);
					ev[100] = ev[93] * ev[10] + ev[15] * ev[86];
					ev[107] = ev[4] * ev[84];
					ev[109] = ev[4] * ev[86];
					ev[116] = ev[93] * ev[15] + ev[15] * ev[93];
					ev[119] = ev[93] * ev[20];
					ev[124] = ev[20] * ev[93];
					ev[129] = ev[20] * ev[84];
					ev[131] = ev[107] * ev[7] - (ev[109] * ev[10] + ev[11] * ev[86]) - ev[116] * ev[7] + (ev[119] * ev[10] + ev[21] * ev[86]) + (ev[124] * ev[10] + ev[24] * ev[86]) - ev[129] * ev[20];
					ev[133] = pow(ev[29],2);

					*(grad+p) = -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[35] * (ev[84] * ev[7] - (ev[86] * ev[10] + ev[10] * ev[86])) + ev[40] * (ev[86] * ev[20] - ev[93] * ev[7]) + ev[46] * (ev[100] - ev[84] * ev[20]))/ev[29] - ev[52] * ev[131]/ev[133] + (ev[40] * (ev[35] * (ev[20] * ev[86] - ev[7] * ev[93]) + ev[46] * (ev[119] - ev[109]))/ev[29] - ev[66] * ev[131]/ev[133]) + (ev[46] * (ev[35] * (ev[100] - ev[129]) + ev[40] * (ev[124] - ev[86] * ev[4]) + ev[46] * (ev[107] - ev[116]))/ev[29] - ev[78] * ev[131]/ev[133])))) + *(theta+9) * (ev[3] * (0.5 * (ev[131] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta5(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[9] = *(theta+8) * *(theta+4);
					ev[10] = ev[9] * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[84] = 2 * *(theta+5);
					ev[92] = *(theta+7) * *(theta+3);
					ev[99] = ev[15] * ev[9];
					ev[107] = ev[4] * ev[9];
					ev[114] = ev[15] * ev[92];
					ev[119] = ev[92] * ev[15];
					ev[124] = ev[92] * ev[5];
					ev[128] = ev[6] * ev[84] - (ev[107] * ev[10] + ev[11] * ev[9]) - ev[16] * ev[84] + (ev[114] * ev[10] + ev[21] * ev[9]) + (ev[119] * ev[10] + ev[24] * ev[9]) - (ev[124] * ev[20] + ev[27] * ev[92]);
					ev[130] = pow(ev[29],2);

					*(grad+p)= -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[35] * (ev[5] * ev[84] - (ev[9] * ev[10] + ev[10] * ev[9])) + ev[40] * (ev[9] * ev[20] + ev[10] * ev[92] - ev[15] * ev[84]) + ev[46] * (ev[99] - ev[5] * ev[92]))/ev[29] - ev[52] * ev[128]/ev[130] + (ev[40] * (ev[35] * (ev[92] * ev[10] + ev[20] * ev[9] - ev[84] * ev[15]) + ev[40] * (ev[4] * ev[84] - (ev[92] * ev[20] + ev[20] * ev[92])) + ev[46] * (ev[114] - ev[107]))/ev[29] - ev[66] * ev[128]/ev[130]) + (ev[46] * (ev[35] * (ev[99] - ev[124]) + ev[40] * (ev[119] - ev[9] * ev[4]))/ev[29] - ev[78] * ev[128]/ev[130])))) + *(theta+9) * (ev[3] * (0.5 * (ev[128] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta6(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[84] = *(theta+4) * *(theta+3);
					ev[85] = ev[84] * ev[10];
					ev[92] = ev[84] * ev[20];
					ev[96] = ev[84] * ev[15] + ev[15] * ev[84];
					ev[99] = ev[20] * ev[84];
					ev[101] = ev[92] * ev[10] - ev[96] * ev[7] + ev[99] * ev[10];
					ev[103] = pow(ev[29],2);

					*(grad+p) = -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[46] * ev[85] - ev[40] * (ev[84] * ev[7]))/ev[29] - ev[52] * ev[101]/ev[103] + (ev[40] * (ev[46] * ev[92] - ev[35] * (ev[7] * ev[84]))/ev[29] - ev[66] * ev[101]/ev[103]) + (ev[46] * (ev[35] * ev[85] + ev[40] * ev[99] - ev[46] * ev[96])/ev[29] - ev[78] * ev[101]/ev[103])))) + *(theta+9) * (ev[3] * (0.5 * (ev[101] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta7(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[84] = *(theta+5) * *(theta+3);
					ev[92] = ev[15] * ev[84];
					ev[94] = ev[84] * ev[15];
					ev[97] = ev[84] * ev[5];
					ev[101] = ev[92] * ev[10] + ev[94] * ev[10] - (ev[97] * ev[20] + ev[27] * ev[84]);
					ev[103] = pow(ev[29],2);

					*(grad+p)= -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[40] * (ev[10] * ev[84]) - ev[46] * (ev[5] * ev[84]))/ev[29] - ev[52] * ev[101]/ev[103] + (ev[40] * (ev[35] * (ev[84] * ev[10]) - ev[40] * (ev[84] * ev[20] + ev[20] * ev[84]) + ev[46] * ev[92])/ev[29] - ev[66] * ev[101]/ev[103]) + (ev[46] * (ev[40] * ev[94] - ev[35] * ev[97])/ev[29] - ev[78] * ev[101]/ev[103])))) + *(theta+9) * (ev[3] * (0.5 * (ev[101] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta8(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[3] = pow(sqrt(2 * M_PI),3);
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[31] = ev[3] * sqrt(ev[29]);
					ev[33] = *(theta+9) * (1/ev[31]);
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[52] = ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]));
					ev[66] = ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]));
					ev[78] = ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]));
					ev[82] = exp(-0.5 * (ev[52]/ev[29] + ev[66]/ev[29] + ev[78]/ev[29]));
					ev[84] = *(theta+4) * *(theta+5);
					ev[92] = ev[15] * ev[84];
					ev[98] = ev[4] * ev[84];
					ev[104] = ev[21] * ev[84] - (ev[98] * ev[10] + ev[11] * ev[84]) + ev[24] * ev[84];
					ev[106] = pow(ev[29],2);

					*(grad+p) = -(ev[33] * (ev[82] * (0.5 * (ev[35] * (ev[40] * (ev[84] * ev[20]) - ev[35] * (ev[84] * ev[10] + ev[10] * ev[84]) + ev[46] * ev[92])/ev[29] - ev[52] * ev[104]/ev[106] + (ev[40] * (ev[35] * (ev[20] * ev[84]) - ev[46] * ev[98])/ev[29] - ev[66] * ev[104]/ev[106]) + (ev[46] * (ev[35] * ev[92] - ev[40] * (ev[84] * ev[4]))/ev[29] - ev[78] * ev[104]/ev[106])))) + *(theta+9) * (ev[3] * (0.5 * (ev[104] * (1/sqrt(ev[29]))))/pow(ev[31],2)) * ev[82]);

					p++;
				}
				a++;
			}
		}
	}
}

void dftheta9(double *theta, int *brain, int *dimx, int *dimy, int *dimz, double *grad) {

    double ev[199];
    int x,y,z,p,a;

    p=0;
    a=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {

				if(brain[a]!=0) {
					ev[4] = pow(*(theta+3),2);
					ev[5] = pow(*(theta+4),2);
					ev[6] = ev[4] * ev[5];
					ev[7] = pow(*(theta+5),2);
					ev[10] = *(theta+8) * *(theta+4) * *(theta+5);
					ev[11] = ev[4] * ev[10];
					ev[15] = *(theta+6) * *(theta+4) * *(theta+3);
					ev[16] = ev[15] * ev[15];
					ev[20] = *(theta+7) * *(theta+5) * *(theta+3);
					ev[21] = ev[15] * ev[20];
					ev[24] = ev[20] * ev[15];
					ev[27] = ev[20] * ev[5];
					ev[29] = ev[6] * ev[7] - ev[11] * ev[10] - ev[16] * ev[7] + ev[21] * ev[10] + ev[24] * ev[10] - ev[27] * ev[20];
					ev[32] = 1/(pow(sqrt(2 * M_PI),3) * sqrt(ev[29]));
					ev[35] = x - *(theta+0);
					ev[40] = y - *(theta+1);
					ev[46] = z - *(theta+2);
					ev[47] = ev[15] * ev[10];
					ev[82] = exp(-0.5 * (ev[35] * (ev[35] * (ev[5] * ev[7] - ev[10] * ev[10]) + ev[40] * (ev[10] * ev[20] - ev[15] * ev[7]) + ev[46] * (ev[47] - ev[5] * ev[20]))/ev[29] + ev[40] * (ev[35] * (ev[20] * ev[10] - ev[7] * ev[15]) + ev[40] * (ev[4] * ev[7] - ev[20] * ev[20]) + ev[46] * (ev[21] - ev[11]))/ev[29] + ev[46] * (ev[35] * (ev[47] - ev[27]) + ev[40] * (ev[24] - ev[10] * ev[4]) + ev[46] * (ev[6] - ev[16]))/ev[29]));

					*(grad+p) = ev[32] * ev[82];

					p++;
				}
				a++;
			}
		}
	}
}


void dfsgtheta0(double *theta,int *dimx, int *dimy, int *dimz, double *grad) {
	double ev[199],theta1,theta2,theta3,theta4,theta5;
	int x,y,z,p;

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {
				theta1=*(theta+0);
				theta2=*(theta+1);
				theta3=*(theta+2);
				theta4=*(theta+3);
				theta5=*(theta+4);

				ev[7] = theta5/(pow(sqrt(2 * M_PI),3) * sqrt(pow(theta4,6)));
				ev[9] = x - theta1;
				ev[11] = pow(theta4,2);
				ev[22] = exp(-0.5 * (pow(ev[9],2)/ev[11] + pow((y - theta2),2)/ev[11] + pow((z - theta3),2)/ev[11]));
				*(grad+p) = ev[7] * (ev[22] * (0.5 * (2 * ev[9]/ev[11])));
				p++;
			}
		}
	}
}

void dfsgtheta1(double *theta,int *dimx, int *dimy, int *dimz, double *grad) {
	double ev[199],theta1,theta2,theta3,theta4,theta5;
	int x,y,z,p;

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {
				theta1=*(theta+0);
				theta2=*(theta+1);
				theta3=*(theta+2);
				theta4=*(theta+3);
				theta5=*(theta+4);

				ev[7] = theta5/(pow(sqrt(2 * M_PI),3) * sqrt(pow(theta4,6)));
				ev[9] = x - theta1;
				ev[11] = pow(theta4,2);
				ev[13] = y - theta2;
				ev[22] = exp(-0.5 * (pow(ev[9],2)/ev[11] + pow((y - theta2),2)/ev[11] + pow((z - theta3),2)/ev[11]));
				*(grad+p) =	ev[7] * (ev[22] * (0.5 * (2 * ev[13]/ev[11])));
				p++;
			}
		}
	}
}

void dfsgtheta2(double *theta,int *dimx, int *dimy, int *dimz, double *grad) {
	double ev[199],theta1,theta2,theta3,theta4,theta5;
	int x,y,z,p;

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {
				theta1=*(theta+0);
				theta2=*(theta+1);
				theta3=*(theta+2);
				theta4=*(theta+3);
				theta5=*(theta+4);

				ev[7] = theta5/(pow(sqrt(2 * M_PI),3) * sqrt(pow(theta4,6)));
				ev[9] = x - theta1;
				ev[11] = pow(theta4,2);
				ev[13] = y - theta2;
				ev[17] = z - theta3;
				ev[22] = exp(-0.5 * (pow(ev[9],2)/ev[11] + pow((y - theta2),2)/ev[11] + pow((z - theta3),2)/ev[11]));
				*(grad+p) = ev[7] * (ev[22] * (0.5 * (2 * ev[17]/ev[11])));
				p++;
			}
		}
	}
}


void dfsgtheta3(double *theta,int *dimx, int *dimy, int *dimz, double *grad) {
	double ev[199],theta1,theta2,theta3,theta4,theta5;
	int x,y,z,p;

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {
				 theta1=*(theta+0);
				 theta2=*(theta+1);
				 theta3=*(theta+2);
				 theta4=*(theta+3);
				 theta5=*(theta+4);

				 ev[3] = pow(sqrt(2 * M_PI),3);
				 ev[4] = pow(theta4,6);
				 ev[6] = ev[3] * sqrt(ev[4]);
				 ev[7] = theta5/ev[6];
				 ev[10] = pow((x - theta1),2);
				 ev[11] = pow(theta4,2);
				 ev[14] = pow((y - theta2),2);
				 ev[18] = pow((z - theta3),2);
				 ev[22] = exp(-0.5 * (ev[10]/ev[11] + ev[14]/ev[11] + ev[18]/ev[11]));
				 ev[24] = 2 * theta4;
				 ev[26] = pow(ev[11],2);
				 *(grad+p) = ev[7] * (ev[22] * (0.5 * (ev[18] * ev[24]/ev[26] + (ev[14] * ev[24]/ev[26] + ev[10] * ev[24]/ev[26])))) - theta5 * (ev[3] * (0.5 * (6 * pow(theta4,5) * 1/sqrt(ev[4]))))/pow(ev[6],2) * ev[22];
				 p++;
			}
		}
	}
}

void dfsgtheta4(double *theta,int *dimx, int *dimy, int *dimz, double *grad) {
	//double ev[199],theta1,theta2,theta3,theta4,theta5;
	double ev[199],theta1,theta2,theta3,theta4;
	int x,y,z,p;

	p=0;
	for(z=1;z<(*dimz+1);z++) {
		for(y=1;y<(*dimy+1);y++) {
			for(x=1;x<(*dimx+1);x++) {
				 theta1=*(theta+0);
				 theta2=*(theta+1);
				 theta3=*(theta+2);
				 theta4=*(theta+3);
				 //theta5=*(theta+4);

				 ev[6] = pow(sqrt(2 * M_PI),3) * sqrt(pow(theta4,6));
				 ev[11] = pow(theta4,2);
				 ev[22] = exp(-0.5 * (pow((x - theta1),2)/ev[11] + pow((y - theta2),2)/ev[11] + pow((z - theta3),2)/ev[11]));
				 *(grad+p) = 1/ev[6] * ev[22];
				 p++;
			}
		}
	}
}




void simTS(double *model, double *mb, double *snr, int *tslen, int *numvox, double *signal, double *weight, double *errors)
{

	int c,r,p;
	//double *tsvec,sigma,tssum,varsum,signal_var,signal_mean,signal_sd,signal_se,signal_weight,sq_tslen;
	double *tsvec,tssum,varsum,signal_var,signal_mean,signal_sd,signal_se,signal_weight,sq_tslen;

	tsvec = (double *) R_alloc(*tslen,sizeof(double));

	//sigma = *mb / *snr;
	sq_tslen = sqrt((double) *tslen - 1);

	GetRNGstate();

	p=0;
	for(r=0;r<(*numvox);r++) {

		tssum=0e0;
		for(c=0;c<(*tslen);c++) {
			tsvec[c]=model[r]+errors[p]*sq_tslen;
			tssum = tssum + tsvec[c];
			p++;
		}

		signal_mean = tssum/(*tslen);

		varsum=0e0;
		for(c=0;c<(*tslen);c++) {
			varsum = varsum + pow(tsvec[c]-signal_mean,2);
		}

		signal_var = varsum/(*tslen - 1);

		signal_sd = sqrt(signal_var);
		signal_se = signal_sd / sq_tslen;
		signal_weight = pow(signal_se,2);

		signal[r]=signal_mean;
		weight[r]=signal_weight;

	}

	PutRNGstate();

}



