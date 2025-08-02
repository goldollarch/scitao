
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

// Rotational shearing interferometer.
// In the rotational shear interferometer the beam interferes with a copy of 
// itself rotated by the angle  around optical axis. Interferometer is useful 
// for detecting of asymmetrical aberrations. 

xbasc();
drawlater();
ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);

N=200;
x=1:N;y=1:N;

F=begin(0.04,0.5e-6,N);
H=circ_ap(F,0.01);

// Rotational interferometer with  defocus aberration
E=zernike(H, 2, 0, 0.01, 10);
G=forvard(E,0.5);
[S1,S2]=b_split(G);
S3=interpol(S2,0.04,N,0.,0.,180.,1.);
S=b_mix(S1,S3);
K=forvard(S,0.5);
I=field_int(K);
subplot(2,2,1);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('defocus aberration')
xtitle(Str);

// Rotational interferometer with  coma
E=zernike(H,3,1,0.01,10);
G=forvard(E,0.5);
[S1,S2]=b_split(G);
S3=interpol(S2,0.04,N,0.,0.,180.,1.);
S=b_mix(S1,S3);
K=forvard(S,0.5);
I=field_int(K);
subplot(2,2,2);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('coma aberration')
xtitle(Str);

// Rotational interferometer with astigmatism
E=zernike(H,2,2,0.01,10);
G=forvard(E,0.5);
[S1,S2]=b_split(G);
S3=interpol(S2,0.04,N,0.,0.,90.,1.);
S=b_mix(S1,S3);
K=forvard(S,0.5);
I=field_int(K);
subplot(2,2,3);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('astigmatism aberration')
xtitle(Str);

// Rotational interferometer with high order aberration
E=zernike(H,7,3,0.01,10);
G=forvard(E,0.5);
[S1,S2]=b_split(G);
S3=interpol(S2,0.04,N,0.,0.,90.,1.);
S=b_mix(S1,S3);
K=forvard(S,0.5);
I=field_int(K);
subplot(2,2,4);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('high order aberration')
xtitle(Str);

drawnow();

