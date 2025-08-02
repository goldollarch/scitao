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

//     The program to form a quadratic 
//     distribution of the refractive index;
//     nc is the grid  sampling
//     an0 is the axis value of refractive index
//     an1 is n1 in the formula for refractive index: 
//     n=sqrt(an0^2-an0*an1*r^2) where 
//     r is the radial coordinate
//     xmax is the grid size

function [] = distribution(nc,an0,an1,xmax)
	fr=mopen('refractive','w');
	fa=mopen('absorptive','w');
	dx=xmax/(nc-1.);
	n2=nc/2+1;
	for i=1:nc
		ii=i-n2
		x=dx*ii
		for j=1:nc
			jj=j-n2
			y=jj*dx
			r2=x*x+y*y
			mfprintf(fr,'%e\n',sqrt(an0*an0 -an1*an0*r2));
			mfprintf(fa,'%e\n',-sqrt(an0*an0 -an1*an0*r2));
		end
		mfprintf(fr,'\n'); 
		mfprintf(fa,'\n'); 
	end
	mclose(fr);
	mclose(fa);
endfunction 
