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

//	a model of an unstable resonator
//	with N=10, M=2, 
//	See D.B. Rench, Applied  Optics 13, 2546...2561 (1974)
//	All dimensions are taken from the table on p. 2552
//	in this reference

//	theoretically the initial field distribution may be random 
//	but we'll use a good plane wave because it won't
//	converge fast with a random distribution 

//     here we use the spherical coordinates

F=begin(7,3e-4,100);
for i=1:20
   G=rect_ap(F,5.48);
   E=l_amplify(G,1e-4,1e4,1);
   D=lens_fresnel(E,-1e4,1e4);
   C=rect_ap(D,10.96);
   B=l_amplify(C,1e-4,1e4,1);
   A=lens_fresnel(B,2e4,1e4);
   S=interpol(A,7);
end
F=convert(A);
out=rect_screen(F,5.48);
field_plot(out);
