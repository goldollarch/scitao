mode(-1);
//
// -------------------------------------------------------------------------
// scitao - Scilab/Scicos Adaptive Optics tooolbox
//
// Copyright (C) 2006  IAPCM , Beijing, China.  Written by
// Chen jingyuan.  For comments or questions about this software,
// please contact the author at jingyuan_chen@yahoo.com.cn.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// -------------------------------------------------------------------------
//
//  Phase reconstruction from a hologram.

m=1;cm=1e-2*m;mm=1e-3*m;nm=1e-9*m;
rad=1;mrad=1e-3*rad;

sz=10*mm;lambda=1000*nm;N=200;
nz=7; mz=1;R=7.2*mm; A=10;theta=2*mrad;

//  First we construct an aberrated beam using Zernike aberration:
Fab=begin(sz,lambda,N);
Fab=zernike(Fab,nz,mz,R,A);
Phase=field_pha(Fab); pha=unf3(1,sz,N,Phase);

// Here we generate the reference beam, mix it with the aberrated beam
// to obtain the interferogram. We tilt the aberrated beam to get vertical fringes:
Fab=tilt(Fab,theta,0);

Fref=begin(sz,lambda,N); Fmixed=b_mix(Fref,Fab);
Iinterferogram=field_int(Fmixed);
ampl=field_int(Fmixed,1);

// To reconstruct the phase of the aberrated beam we first substitude the
// hologram in the (reference) beam:
Frec=begin(sz,lambda,N); 
[header,r,imaginary]=field_contents(Frec);
Frec=create_field(header,ampl,imaginary);

// Then we Fourier transform the field and shift the spectral domain:
Frec=pip_fft(Frec,1); FrecShifted=interpol(Frec,sz,N,1*mm,0,0,1);

tmp_f=interpol(Frec,4e-3);file_ps(tmp_f,"fft_shi",-1,200,14);
tmp_f=rect_ap(tmp_f,15e-4,1);file_ps(tmp_f,"fft_sf",-1,200,14);

// Next we filter the beam in the Fourier domain and apply a reverse Fourier
// transformation:
FrecShifted=rect_ap(FrecShifted,1.5*mm,1,0,0,0);
FrecShifted=pip_fft(FrecShifted,-1);

// Finally we calculate the phase and unwrap it
pha=field_pha(FrecShifted); Phaserec=unf3(-1,sz,N,pha);

xbasc();drawlater();
subplot(1,3,1);grayplot(1:N,1:N,Iinterferogram);
subplot(1,3,2);surf(Phase);subplot(1,3,3);surf(Phaserec);
drawnow();