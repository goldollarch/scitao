/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation of 
wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2000-2006 IAPCM, Beijing, China.  Written by
Chen jingyuan.  For comments or questions about this software,
please contact the author at jingyuan_chen@yahoo.com.cn.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

/* 
    Based on the software LightPipes, there are some modifications to accord 
with the Scilab/Scicos environment. the author thanks Dr. Gleb Vdovin 
for developing the excellent software package LightPipes.
*/

/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include "lightPipes.h"

typedef struct {
    double r,i;
} fcomplex;

fcomplex Cadd();
fcomplex Csub();
fcomplex Cmul();
fcomplex Complex();
fcomplex Conjg();
fcomplex Cdiv();
fcomplex Csqrt();
fcomplex RCmul();
double Cabs();

void all_mem(), free_mem(), error_print(), elim();

fcomplex * a, * b, * c, * alpha, * beta, * u, * p, * u1, * u2;
double * refr, * absorb;

FIELD lp_steps( FIELD field,  double z, int nstep, 
		char *rfname, char *afname, char *ofname, int n_st )
{
    long ik, ik1, ikij,  ikij1,  ikij_1, iki1j, iki_1j;
    double  delta, delta2, Pi4lz,A,AA,band_pow,K,dx, dist;
    fcomplex uij, uij1, uij_1, ui1j, ui_1j, medium;
    int i,j,jj, ii, istep, i_left, i_right;

    FILE * fr;
    double * int1;
    double int2;
    double * phase1;
    double phase2;
	
	all_mem(field);

    K=2.*Pi/field.lambda;

	if((strstr(rfname, "void")) != NULL ) {
		printf(" void refractive coefficient file, skipping \n");
		goto l1;
	}
	if((fr=fopen(rfname,"r"))==NULL){
		printf(" error opening file %s",rfname );
		exit(1);
	} 

    ik=0;  
    for (i=1;i<= field.number; i ++){
        for (j=1;j <= field.number;j ++){
			double fi;
			if ((fscanf(fr,"%le",&fi))==EOF){
				printf(" reading the refractive indices:\
					end of input file reached, exiting\n");
				exit(1);
			}
	      
			refr[ik]=fi;
			ik ++;
		}      
	}
	
	fclose(fr);

l1: 
	
	if((strstr(afname, "void")) != NULL ) {
		printf(" void absorption coefficient file, skipping \n");
		goto l2;
	}

	if((fr=fopen(afname,"r"))==NULL){
		printf(" error opening file %s\n",afname);
		exit(1);
	} 

    ik=0;  
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
			double fi;
			if ((fscanf(fr,"%le",&fi))==EOF){
				printf(" absorption  end of input file reached, exiting\n");
                exit(1);
			}
			absorb[ik]=fi;
			ik ++;
		}
    }
	
	fclose(fr);

l2:
	if(ofname!=NULL)
		if((fr=fopen(ofname,"w"))==NULL){
			printf(" error create file %s\n",ofname );
			exit(1);
		} 

    int1=(double *) calloc(field.number+2, sizeof(double));
    if(int1 == NULL){
		printf("Allocation error in int1\n");
		exit(1);
    }

    phase1=(double *) calloc(field.number+2, sizeof(double));
    if(phase1 == NULL){
		printf("Allocation error in phase1\n");
		exit(1);
    }

	dx=field.size/(field.number-1.);
	n_st=1;

	/* the arrays for refraction and absorption are finally formed */

    if (field.double1 !=0.) {
        printf(" can not be applied in spherical\
					   coordinates,\nuse CONVERT first\n");
        exit(1);
    }
    
    z=z/2.;
    Pi4lz = 4.*Pi/field.lambda/z;
    delta=field.size/((double)(field.number-1.));
    delta2 = delta*delta;
    A=0.;
	
	/* absorption at the borders is described here */
    AA= -10./z/nstep; /* total absorption */
    band_pow=2.;   /* profile of the absorption border, 2=quadratic*/
	
	/* width of the absorption border */
    i_left=field.number/2+1-0.4*field.number;
    i_right=field.number/2+1+0.4*field.number;
	/* end absorption */
  

    for (i=1; i<=field.number; i++){
      u2[i].r=u2[i].i = 0.;
	  /*      a[i].r=b[i].r= -1./delta2; */
      a[i].i=b[i].i=0.;
      a[i].r=b[i].r= -1./delta2;
    }

	medium.r=0.;
    dist =0.;
    
	/*  Main  loop, steps here */
    for(istep = 1; istep <= nstep ; istep ++){
		dist=dist + 2.*z;

		/*  Elimination in the direction i, halfstep  */

		ik=0;
		
		for (i=1; i<= field.number; i++){
			for(j=1; j<= field.number; j++){
				double cab, sab, fi, cc;
				fi= 0.25*K*z*(refr[ik]-1.);
				cab=cos(fi);
				sab=sin(fi);
				cc=field.real[ik]*cab-field.imaginary[ik]*sab;
				field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
				field.real[ik]=cc;
				ik++;
			}
		}

		for(jj=2; jj<= field.number-2; jj += 2){
			j=jj;
			for (i=2; i<=field.number-1; i++){
				ikij=(j-1)*field.number+i-1;
				ikij1=(j)*field.number+i-1;
				ikij_1=(j-2)*field.number+i-1;

				uij=Complex(field.real[ikij], field.imaginary[ikij]);
				uij1=Complex(field.real[ikij1], field.imaginary[ikij1]);
				uij_1=Complex(field.real[ikij_1], field.imaginary[ikij_1]);

				p[i] = RCmul(-2., uij);
				p[i] = Cadd(p[i], uij1);
				p[i] = Cadd(p[i], uij_1);
				p[i] = RCmul(-1./delta2, p[i]);
				p[i] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[i]);
			}

			for (i=1; i<=field.number; i++){
				ik=(j-1)*field.number + i -1;
				medium.i= -2.*Pi*absorb[ik]/field.lambda;
				c[i].r=-2./delta2;
				c[i].i= Pi4lz;
				c[i].i += medium.i;  
				/* absorption borders are formed here */
				if( i <= i_left){
					int iii=i_left-i;
					c[i].i -= (AA*K)*pow((double) iii/ ((double)(i_left)),band_pow);
				} 
				if( i >= i_right){ 
					int iii =i-i_right;
					int im=field.number-i_right;
					c[i].i -= (AA*K)*pow((double) iii/ ((double)(im)),band_pow);
				}
				
				/* end absorption */
			}
			
			elim(field);

			for (i=1; i<=field.number; i++){
				ikij_1=(j-2)*field.number+i-1;
				field.real[ikij_1]=u2[i].r;
				field.imaginary[ikij_1]=u2[i].i;
				u1[i].r=u[i].r;
				u1[i].i=u[i].i;
			}

			j=jj+1;
			
			for (i=2; i<=field.number-1; i++){
				ikij=(j-1)*field.number+i-1;
				ikij1=(j)*field.number+i-1;
				ikij_1=(j-2)*field.number+i-1;

				uij=Complex(field.real[ikij], field.imaginary[ikij]);
				uij1=Complex(field.real[ikij1], field.imaginary[ikij1]);
				uij_1=Complex(field.real[ikij_1], field.imaginary[ikij_1]);
				
				p[i] = RCmul(-2., uij);
				p[i] = Cadd(p[i], uij1);
				p[i] = Cadd(p[i], uij_1);
				p[i] = RCmul(-1./delta2, p[i]);
				p[i] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[i]);
			}
			
			for (i=1; i<=field.number; i++){
				ik=(j-1)*field.number + i -1;
				medium.i= -2.*Pi*absorb[ik]/field.lambda;
				c[i].r=-2./delta2;
				c[i].i= Pi4lz;
				c[i].i += medium.i;   

				/* absorption borders are formed here */
				if( i <= i_left){
					int iii=i_left-i;
					c[i].i -= (AA*K)*pow((double) iii/ ((double)(i_left)),band_pow);
				} 
				
				if( i >= i_right){ 
					int iii =i-i_right;
					int im=field.number-i_right;
					c[i].i -= (AA*2.*K)*pow((double) iii/ ((double)(im)),band_pow);
				}
				
				/* end absorption */
			}
			
			elim(field);

			for (i=1; i<=field.number; i++){
				ikij=(j-2)*field.number+i-1;
				field.real[ikij]=u1[i].r;
				field.imaginary[ikij]=u1[i].i;
				
				u2[i].r=u[i].r;
				u2[i].i=u[i].i;
			}
		}

		for (i=1; i<=field.number; i++){
			ikij=(field.number-1)*field.number+i-1;
			field.real[ikij]=u2[i].r;
			field.imaginary[ikij]=u2[i].i; 
		} 
		
		ik=0;
		for (i=1; i<= field.number; i++){
			for(j=1; j<= field.number; j++){
				double cab, sab, fi, cc;
				fi= 0.5*K*z*(refr[ik]-1.);
				cab=cos(fi);
				sab=sin(fi);
				cc=field.real[ik]*cab-field.imaginary[ik]*sab;
				field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
				field.real[ik]=cc;
				ik++;
			}
		}

		/* Elimination in the j direction is here, halfstep */

		for (i=1; i<=field.number; i++){
			u2[i].r=u2[i].i = 0.;
		}

		for(ii=2; ii<= field.number-2; ii += 2){
			i=ii;
			for (j=2; j<=field.number-1; j++){
				ikij=(j-1)*field.number+i-1;
				iki1j=(j-1)*field.number+i;
				iki_1j=(j-1)*field.number+i-2;

				uij=Complex(field.real[ikij], field.imaginary[ikij]);
				ui1j=Complex(field.real[iki1j], field.imaginary[iki1j]);
				ui_1j=Complex(field.real[iki_1j], field.imaginary[iki_1j]);
				
				p[j] = RCmul(-2., uij);
				p[j] = Cadd(p[j], ui1j);
				p[j] = Cadd(p[j], ui_1j);
				p[j] = RCmul(-1./delta2, p[j]);
				p[j] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[j]);
			}
			
			for (j=1; j<=field.number; j++){
				ik=(j-1)*field.number + i -1;
				medium.i= -2.*Pi*absorb[ik]/field.lambda;
				c[j].r=-2./delta2;
				c[j].i= Pi4lz;
				c[j].i += medium.i;  
				/* absorption borders are formed here */
				if( j <= i_left){
					int iii=i_left-j;
					c[j].i -= (AA*K)*pow((double) iii/ ((double)(i_left)),band_pow);
				} 
				
				if( j >= i_right){ 
					int iii =j-i_right;
					int im=field.number-i_right;
					c[j].i -= (AA*K)*pow((double) iii/ ((double)(im)),band_pow);
				}
				
				/* end absorption */
			}
			
			elim(field);
			
			for (j=1; j<=field.number; j++){
				iki_1j=(j-1)*field.number+i-2;
				field.real[iki_1j]=u2[j].r;
				field.imaginary[iki_1j]=u2[j].i;
				u1[j].r=u[j].r;
				u1[j].i=u[j].i;      
			}

			i=ii+1;
			for (j=2; j<=field.number-1; j++){
				ikij=(j-1)*field.number+i-1;
				iki1j=(j-1)*field.number+i;
				iki_1j=(j-1)*field.number+i-2;
				uij=Complex(field.real[ikij], field.imaginary[ikij]);
				ui1j=Complex(field.real[iki1j], field.imaginary[iki1j]);            ui_1j=Complex(field.real[iki_1j], field.imaginary[iki_1j]);

				p[j] = RCmul(-2., uij);
				p[j] = Cadd(p[j], ui1j);
				p[j] = Cadd(p[j], ui_1j);
				p[j] = RCmul(-1./delta2, p[j]);
				p[j] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[j]);
			}
			
			for (j=1; j<=field.number; j++){
				ik=(j-1)*field.number + i -1;
				medium.i= -2.*Pi*absorb[ik]/field.lambda;
				c[j].r=-2./delta2;
				c[j].i= Pi4lz;
				c[j].i += medium.i;  

				/* absorption borders are formed here */
				if( j <= i_left){
					int iii=i_left-j;
					c[j].i -= (AA*K)*pow((double) iii/ ((double)(i_left)),band_pow);
				} 
                
				if( j >= i_right){ 
					int iii =j-i_right;
					int im=field.number-i_right;
					c[j].i -= (AA*K)*pow((double) iii/ ((double)(im)),band_pow);
				}
				
				/* end absorption */
			}
			
			elim(field);

			for (j=1; j<=field.number; j++){
				ikij=(j-1)*field.number+i-2;
				field.real[ikij]=u1[j].r;
				field.imaginary[ikij]=u1[j].i;
				u2[j].r=u[j].r;
				u2[j].i=u[j].i;
			}
		}

		for (j=1; j<=field.number; j++){
			ikij=(j-1)*field.number+i-1;
			field.real[ikij]=u2[j].r;
			field.imaginary[ikij]=u2[j].i;
		}
    
		/* end j */

		if ( (ofname!=NULL) && (istep/n_st == (float) istep / (float) n_st ) )
		{
			/* writing the intensity into arrays  */
			i=field.number/2+1;
			for (j=1;j<=field.number;j += 1){ 
				ik1=(i-1)*field.number+j-1;
				jj=j-1; 
				int1[jj]=field.real[ik1]*field.real[ik1]+field.imaginary[ik1] *field.imaginary[ik1]; 
				phase1[jj]=phase(field.imaginary[ik1],field.real[ik1]);
			}

			j=field.number/2+1;
			for (i=1;i<=field.number;i += 1){
				double cc;
				ik1=(i-1)*field.number+j-1;
				jj=i-1;
				cc=dx*(i-field.number/2-1);
				int2=field.real[ik1]*field.real[ik1]+field.imaginary[ik1] *field.imaginary[ik1]; 
				phase2=phase(field.imaginary[ik1], field.real[ik1]);
				fprintf(fr," %e %e %e %e %e %e\n", cc, int1[jj], int2, phase1[jj], phase2, dist);
			}
		}
		
		if ( (ofname!=NULL) && (istep/n_st == (float) istep / (float) n_st ) )
		{
			fprintf(fr,"\n");
		}
	}

	ik=0;
    for (i=1; i<= field.number; i++){
		for(j=1; j<= field.number; j++){
			double cab, sab, fi, cc;
            fi=0.25*K*z*(refr[ik]-1.);
            cab=cos(fi);
            sab=sin(fi);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
            ik++;
		}
	}

    if(ofname!=NULL) fclose(fr);

    free_mem();

    free(int1);
	free(phase1);

	return ( field );

}


/////////////////////////////////////////////////////////////////////////

void elim(FIELD field)
{ 
    int i;
    fcomplex cc;

    /* initial condition, everything is going to be zero at the edge */
    alpha[2].r=0.;
    alpha[2].i=beta[2].r=beta[2].i=0.;

    alpha[field.number].r=0.;
    alpha[field.number].i=beta[field.number].r=beta[field.number].i=0.;

    /* forward elimination */
    for(i=2;i<=field.number-2; i++){
        cc=Csub(c[i],Cmul(a[i],alpha[i]));
        alpha[i+1]=Cdiv(b[i],cc);
        beta[i+1]=Cdiv(Cadd(p[i],Cmul(a[i],beta[i])),cc);
    }

    i=field.number;
    cc=Csub(c[i],Cmul(a[i],alpha[i]));
    beta[field.number+1]=Cdiv(Cadd(p[i],Cmul(a[i],beta[i])),cc);

    /* edge amplitude =0 */
    u[field.number]=beta[field.number+1];

    /* backward elimination        */
    for(i=field.number-1; i>=1; i--){
        u[i]=Cadd(Cmul(alpha[i+1],u[i+1]),beta[i+1]);

    }
}

void all_mem(FIELD field)
{ 
	int i,j; 
	long ik;

    /*    Allocatiom of the  memory for arrays: */
	refr=(double *) calloc( (field.number+2)*(field.number+2), sizeof(double) );
	absorb=(double *) calloc( (field.number+2)*(field.number+2), sizeof(double) ); 

	ik=0;
    for(i=1; i<= field.number+2; i++)	{
		for(j=1; j<= field.number+2; j++)	  {
		  refr[ik]=1.;
		  absorb[ik]=0.;
		  ik ++;
      }
	}

	a=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	p=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	b=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	c=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	u=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	u1=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	u2=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	alpha=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
	beta=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) );
}

void free_mem()
{
    /* freeing  the memory        */
    free(a);
    free(b);
    free(c);
    free(u);
    free(u1);
    free(u2);
    free(alpha);
    free(beta);
    free(p);
    free(absorb);
    free(refr);

}

////////////////////////////////////////////////////////
/* complex number handling is here, not optimal ... : */

fcomplex Cadd(fcomplex a,fcomplex b)
{
    fcomplex c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplex Csub(fcomplex a,fcomplex b)
{
    fcomplex c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplex Cmul(fcomplex a,fcomplex b)
{
    fcomplex c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplex Complex(double re,double im)
{
    fcomplex c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplex Conjg(fcomplex z)
{
    fcomplex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplex Cdiv(fcomplex a,fcomplex b)
{
    fcomplex c;
    double den;

    den=b.r*b.r+b.i*b.i;

    if( den != 0.){
        c.r= (a.r*b.r +a.i*b.i)/den;
        c.i= (a.i*b.r - a.r*b.i)/den;
    }
    else{
        exit(1);
    }

    return c;
}

double Cabs(fcomplex z)
{
    double x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    }
    else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplex Csqrt(fcomplex z)
{
    fcomplex c;
    double x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    }
    else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        }
        else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        }
        else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplex RCmul(double x,fcomplex a)
{
    fcomplex c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}
