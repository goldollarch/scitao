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
    Based on the software LightPipes, there are some changes to accord 
with the Scilab/Scicos environment. the author thanks Dr. Gleb Vdovin 
for developing the excellent software package LightPipes.
*/

#ifndef LIGHTPIPES_H
#define LIGHTPIPES_H

#include "fftn.h"
#include "pipes.h"

#include "string.h"

void fft3(int,FIELD);

void forvard(double,FIELD);
int fresnl(double,double*,double*);
void fresnel(double,FIELD);

void lens(double,double,double,FIELD);

FIELD lp_begin(double,double,int,int);

FIELD lp_absorber(FIELD,double);
FIELD lp_cros_out(FIELD,char*,double);

FIELD lp_circ_ap(FIELD,double,double,double);
FIELD lp_circ_screen(FIELD,double,double,double);

FIELD lp_gauss(FIELD,double,double,double,double);
FIELD lp_gauss_screen(FIELD,double,double,double,double);

FIELD lp_rect_ap(FIELD,double,double,double,double,double);
FIELD lp_rect_screen(FIELD,double,double,double,double,double);

FIELD lp_random(FIELD,char*,double,unsigned int);

FIELD lp_file_pgm(FIELD,char*,double,int,int,double);
FIELD lp_file_ps(FIELD,char*,double,int,double);
FIELD lp_file_int(FIELD,char*,int,double,int);
FIELD lp_file_pha(FIELD,char*,int,double,int);

FIELD lp_forvard(FIELD,double);
FIELD lp_fresnel(FIELD,double);
double *lp_forward(FIELD,double,int,double);

FIELD lp_zernike(FIELD,int,int,double,double);
FIELD lp_l_amplif(FIELD,double,double,double);

FIELD lp_lens(FIELD,double,double,double);
FIELD lp_lens_fresn(FIELD,double,double);
FIELD lp_lens_forvard(FIELD,double,double);

FIELD lp_normal(FIELD);
FIELD lp_convert(FIELD);

FIELD lp_tilt(FIELD,double,double );
FIELD lp_tor_lens(FIELD,double,double,double,double);

FIELD lp_pip_fft(FIELD,int);

double *lp_interp1(FIELD,double,int,double,double,double,double);
double *lp_interpol(FIELD,double,int,double,double,double,double);
double *lp_b_split(FIELD,double );

FIELD *b_split(FIELD,double );

FIELD lp_b_mix(FIELD,FIELD);

FIELD lp_fil_ter(FIELD,char*,char*,char*,int);
FIELD lp_steps(FIELD,double,int,char*,char*,char*,int);

double *lp_strehl(FIELD);
double *lp_field_int( double*,int);
double *lp_field_pha( double*);

double *unf3(double,double,int,double*);
double *unf4(double,int,double*);
void c_scilab(char*,char*);

#endif
