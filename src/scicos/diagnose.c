/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation 
of wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2005-2006 IAPCM, Beijing, China.  Written by
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
with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include "optics_scicos.h"

void diagnose(scicos_block *block,int flag)
{
	double t;
	int *ipar;

	char str[20],*fname;
	int job = 1;
	FILE *fr;

    int i,j, n2;
    double sum, sum1r, sum1i, sum1, sum2, dx,dx2, s_p, x,y, x_c, y_c;
    double sum1x, sum1y;
    long ik1;

	FIELD field;

	/*   ipar[1]= lfil : file name length  */
	/*   ipar[2:1+lfil]:  character codes for file name */

    ipar=block->ipar;
	--ipar;

	F2C(cvstr)(&(ipar[1]),&(ipar[2]), str,&job,strlen(str));
	str[ipar[1]] = '\0';

	t=get_scicos_time();

	if (flag == 1) {

		field=block_inptr_field(block,0);

		fname=current_filename(str,t,".txt");
		fr=fopen(fname,"w");
 
		fprintf(fr,"Grid size: %e\n", field.size); 
		fprintf(fr,"Grid sampling: %d\n", field.number); 
		fprintf(fr,"Wave length: %e\n\n", field.lambda); 
		fprintf(fr,"time: %e\n", t); 

		dx =field.size/(field.number);
		dx2 = dx*dx;
		n2=field.number/2+1;

		/* Calculating the power */
		sum=sum1r=sum1i=sum2=0.;
		ik1=0;
		for (i=1;i<=field.number ;i++){
			for (j=1;j<=field.number ;j++){
				s_p=(field.real[ik1]*field.real[ik1]+ 
					field.imaginary[ik1]*field.imaginary[ik1]);
				sum2 += s_p;
				sum += sqrt(s_p);
				sum1r += field.real[ik1];
				sum1i += field.imaginary[ik1];
				ik1++;
			}
		}

		sum1=(sum1r*sum1r+sum1i*sum1i);
		
		if (sum == 0) {
			sciprint("Strehl: Zero beam power, program terminated\n");
			exit(1);
		}

		fprintf(fr,"energy : %e\n", sum2*dx2);
		fprintf(fr,"Strehl ratio : %e\n",sum1/sum/sum);

		/* Calculating the center of gravity: */
		sum=sum1r=sum1i=sum2=0.;
		ik1=0;
		for (i=1;i<=field.number ;i++){
			y=(i-n2)*dx;
			for (j=1;j<=field.number ;j++){
				x=(j-n2)*dx;
				sum2=field.real[ik1]*field.real[ik1]+
					field.imaginary[ik1]*field.imaginary[ik1];
				sum1r += sum2*x;
				sum1i += sum2*y;
				sum += sum2;
				ik1++;
			}
		}

		x_c=sum1r/sum;
		y_c=sum1i/sum;

		fprintf(fr,"Center of gravity: x= %e y= %e\n", x_c, y_c);

		/* Calculating moments of the distribution */
		sum1r=sum1x=sum1y=0.;
		ik1=0;
		for (i=1;i<=field.number ;i++)
		{
			double y_y_c;
			y=(i-n2)*dx;
			y_y_c=y-y_c;

			for (j=1;j<=field.number ;j++)
			{
				double temp_int, x_x_c;
				x=(j-n2)*dx;
				x_x_c=x-x_c;
				temp_int = field.real[ik1]*field.real[ik1] +
					field.imaginary[ik1]*field.imaginary[ik1];
				sum1r += temp_int*(x_x_c*x_x_c+y_y_c*y_y_c);
				sum1x += temp_int*(x_x_c*x_x_c);
				sum1y += temp_int*(y_y_c*y_y_c);
				ik1++;
			}
		}

		fprintf(fr,"Standard deviation:  S_r=%e S_x= %e S_y= %e\n",
			sqrt(sum1r/sum), sqrt(sum1x/sum), sqrt(sum1y/sum));
		fclose(fr);

		field_block_outptr(field,block,0);

	}

}
