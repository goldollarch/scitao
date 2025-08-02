
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

//	We form here the distribution of the refractive index absorbtion coefficients
//	in the file "refractive" and "absorptive"

exec("distribution.sci");
distribution(100,1.5,400,1e-3);

// 	propagation of a tilted Gaussian beam beam in a quadratic medium: 
//	the steps back are implemented with the 
//	second operator to show the reversibility of the
//	operator:

N=100;
f1=begin(1e-3,1e-6,N);
f2=gauss(f1,113e-6);
f3=tilt(f2,1e-3,0);

f4=steps(f3,1e-3,1000,"refractive","absorptive","out_1",10);
out1=field_int(f4);
//steps_out_plot("out_1",100,100,2);

f5=steps(f4,-1e-3,1000,"refractive","absorptive","out_2",10);
out2=field_int(f5);
  
// here is the  propagation of a non-axial (shifted) gauss beam:

f2=gauss(f1,113e-6,2e-4);
f3=steps(f2,1e-3,1000,"refractive","absorptive","out_3",10);
out3=field_int(f3);

// here we have two non-axial beams:

f3=gauss(f1,113e-6,-2e-4);
f4=b_mix(f2,f3);
f5=steps(f4,1e-3,1000,"refractive","absorptive","out_4",10);
out4=field_int(f5);

xbasc();
drawlater();

subplot(2,2,1);
grayplot(1:N,1:N,out1,strf="030",rect=[1,1,N,N]);

subplot(2,2,2);
grayplot(1:N,1:N,out2,strf="030",rect=[1,1,N,N]);

subplot(2,2,3);
grayplot(1:N,1:N,out3,strf="030",rect=[1,1,N,N]);

subplot(2,2,4);
grayplot(1:N,1:N,out4,strf="030",rect=[1,1,N,N]);

drawnow();

