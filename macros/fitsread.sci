//
// -------------------------------------------------------------------------
// scitao - Scilab/Scicos Adaptive Optics tooolbox
//
// Copyright (C) 2006  IAPCM , Beijing, China.  Written by
// gan guangyong.  For comments or questions about this file,
// please contact the author at .
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

function data = fitsread(varargin)
    data = []; 

		[filename,extension,index,raw] = parseInputs(varargin);
		info = fitsinfo(filename);
			
    select extension, 
     case 'primary'     then
		  data = readprimary(info,raw),
		 case 'ascii'       then
		//  data = readasciitable(info,index,raw),
		warning('Problem reading extension ascii. The codes are not tested.');		
		 case 'binary'      then
		//  data = readbinarytable(info,index,raw),
		warning("Problem reading extension binary. The codes are not tested.");			
		 case 'image'       then
		  data = readimage(info,index,raw),
		 case 'unknown'     then
		  data = readunknown(info,index,raw),
		 else               
		  error('EXTENSION must be ''Primary'', ''Table'', ''BinTable'', or ''Unknown''.'),
		end    //select

endfunction    //fitsread

function [filename,extension,index,raw] = parseInputs(varargin)
			//Verify inputs are correct
			//error(nargchk(1,4,nargin));
			right_num = length(varargin(1));
			if right_num<1. | right_num>4. then
			    error('Number of argument is wrong!');
			end    
			
			filename = varargin(1)(1);
			extension = 'primary' ;    // varargin(1)(2)
			index = 1.;
			raw = 0.;

			allStrings = ['primary','image','table','bintable','unknown','raw'];
	    k = 2.;  //for
			 while k<=right_num  then 
			  if length(varargin(1)(k)) > 1   then
			    lit_var = convstr(varargin(1)(k),'l');	
					  select lit_var,   
					   case 'primary'     then			   
					    extension = 'primary';
					    if (k+1)<=right_num & (length(varargin(1)(k+1))==1.)  then
                  warning('File only one primary.');
					    end					    
					   case 'bintable'    then
					    extension = 'binary';
					    if (k+1)<=right_num & (length(varargin(1)(k+1))==1.)  then
					      index  = varargin(1)(k+1);
					    end
					   case 'image'       then
					    extension = 'image';
					    if (k+1)<=right_num & (length(varargin(1)(k+1))==1.)  then
					      index  = varargin(1)(k+1);
					    end
					   case 'table'       then
					    extension = 'ascii';
					    if (k+1)<=right_num & (length(varargin(1)(k+1))==1.)  then
					      index  = varargin(1)(k+1);
					    end
					   case 'unknown'     then
					    extension = 'unknown';
					    if  (k+1)<=right_num & (length(varargin(1)(k+1))==1.)  then
					      index  = varargin(1)(k+1);
					    end
					   case 'raw'         then
					    raw = 1.;
					  end  //select

			  else
			    //Don't allow fitsread(filename,idx);
			    if k==2.    then
			      error('The extension index IDX must follow the extension name EXTNAME.');
			    end   //if
			  end   //if
			  k=k+1.;
			end   //while
				
endfunction   //   parseInputs

function data = readprimary(info,raw)
			//Read data from primary data 
			data = [];
			msg = 'Error reading file.  File may be an invalid FITS file or may be corrupt.';
			
			if info.PrimaryData.DataSize == 0. then
			  return;
			end
			
			startpoint = info.PrimaryData.Offset;
			
			//Data will be scaled by scale values BZERO, BSCALE if they exist
			bscale = info.PrimaryData.Slope;
			bzero = info.PrimaryData.Intercept;
			nullvals = info.PrimaryData.MissingDataValue;
			
			fid = mopen(info.Filename,'rb');
			if fid == -1   then
			  error(msg);
			end
			mseek(startpoint,fid,'set');
			
			data = mget(prod(info.PrimaryData.Size), precision2scilab(info.PrimaryData.DataType),fid);
			mclose(fid);
			if length(data) < prod(info.PrimaryData.Size)  then
			  warning('Problem reading primary data. Data has been truncated.');
			else
			  //Data is stored in column major order so the first two dimensions must be permuted.
			  data = permute(matrix(data,info.PrimaryData.Size),[2 1 3:length(info.PrimaryData.Size)]);
			  //Scale data and replace undefined data with NaN by default
			  if ~raw & ~isempty(nullvals)  then
			    data(data==nullvals) = NaN;
			  end
			  if ~raw then
			    data = double(data)*bscale+bzero;
			  end
			end
endfunction   //END READFITSPRIMARY

function data = readimage(info,index,raw)
			//Read data from image extension
			
			data = [];
			msg = 'Error reading file.  File may be an invalid FITS file or may be corrupt.';
			
			if ~length(info.Image)   then
			  error('File does not contain any Image Extensions.');
			elseif length(info.Image) < index   then
			  error('File only contains Image extensions.');
			end
			
			if info.Image(index).DataSize == 0.  then
			  //No data
			  return;
			end
			
			//Data will be scaled by scale values BZERO, BSCALE if they exist
			bscale = info.Image(index).Slope;
			bzero = info.Image(index).Intercept;
			nullvals = info.Image(index).MissingDataValue;
			
			startpoint = info.Image(index).Offset;
			
			fid = mopen(info.Filename,'r');
			if fid==-1  then
			  error(msg);
			end
			mseek(startpoint,fid,'set');
			data = mgeti(prod(info.Image(index).Size),precision2scilab(info.Image(index).DataType),fid);
			mclose(fid);
			if length(data)<prod(info.Image(index).Size)   then
			  warning('Problem reading image data. Data has been truncated.');
			else
			  //Data is stored in column major order so the first two dimensions must be permuted.
			  data = permute(matrix(data,info.Image(index).Size),[2 1 3:length(info.Image(index).Size)]);
			  //Scale data and replace undefined data with NaN by default
			  if ~raw & ~isempty(nullvals)   then
			    data(data==nullvals) = NaN;
			  end
			end
endfunction  //END READFITSIMAGE


function data = readbinarytable(info,index,raw)
			//Read data from binary table
			
			data = [];
			msg = 'Error reading binary table extension.  File may be an invalid FITS file or may be corrupt.';
			
			if ~length(info.BinaryTable)   then
			  error('File does not contain any Binary Table Extensions.');
			  return;
			elseif length(info.BinaryTable)< index   then
			  error('File only contains  Binary Table extensions.');
			  return;
			end
			
			if info.BinaryTable(index).DataSize == 0.   then
			  // No data
			  return;
			end
			
			tscal = info.BinaryTable(index).Slope;
			tzero = info.BinaryTable(index).Intercept;
			nullvals = info.BinaryTable(index).MissingDataValue;
			
			startpoint = info.BinaryTable(index).Offset;
			
			fid = mopen(info.Filename,'r');
			if fid==-1   then
			  error(msg);
			end
			mseek(startpoint,fid,'set');
			data=list();
			m=1;
	while m	<=	info.BinaryTable(index).NFields
	    data(m)=list();
	    m=m+1;
	end    
		
			//Read data. Take care of complex data and scaling.
			i=1;j=1;
			while i<=info.BinaryTable(index).Rows
			  while j<=info.BinaryTable(index).NFields
			        
			    //Field has no data if size is zero
			    if info.BinaryTable(index).FieldSize(j) == 0.   then
			      data(j)(i) = zeros(0,1);
			      break;
			    end
			    
			    //[precision, cmplx] = strtok(info.BinaryTable(index).FieldPrecision{j});
			    precision = msscanf(info.BinaryTable(index).FieldPrecision(j),'%s %*s');
			    cmplx = msscanf(info.BinaryTable(index).FieldPrecision(j),'%*s %s');
			    
			    if isempty(cmplx)    then
			      fielddata = mget(info.BinaryTable(index).FieldSize(j),precision2scilab(precision),fid);
			    else
			      fielddata = mget([info.BinaryTable(index).FieldSize(j),2],precision2scilab(precision),fid);
			    end
	    
			    if ~raw & ~isempty(nullvals(j))    then
			      fielddata(fielddata==nullvals(j)) = %nan;
			    end
			    
			    if ~isempty(cmplx)   then
			      fielddata = fielddata(:,1)+fielddata(:,2)*%i;
			    end
		    j=j+1;
			  end
			  i=i+1;
			end
			mclose(fid);
endfunction     //End READFITSBINARYTABLE


function data = readunknown(info,index,raw)
			//Read data from unknown data 
			
			data = [];
			msg = 'Error reading file.  File may be an invalid FITS file or may be corrupt.';
			
			if length(info.Unknown) == 0.  then
			  error('File does not contain any Unknown Extensions.');
			elseif length(info.Unknown)< index  then
			  error('File only contains Unknown extensions.');
			end
			
			if info.Unknown(index).DataSize == 0.  then
			  return;
			end
			
			startpoint = info.Unknown(index).Offset;
			
			//Data will be scaled by scale values BZERO, BSCALE if they exist
			bscale = info.Unknown(index).Slope;
			bzero = info.Unknown(index).Intercept;
			nullvals = info.Unknown(index).MissingDataValue;
			
			fid = mopen(info.Filename,'r');
			if fid==-1  then
			  error(msg);
			end
			mseek(startpoint,fid,'set');
			data = mget(prod(info.Unknown(index).Size),precision2scilab(info.Unknown(index).DataType),fid);
			mclose(fid);
			if length(data) < prod(info.Unknown(index).Size)  then
			  warning('Problem reading data. Data has been truncated.');
			else
			  data = permute(matrix(data,info.Unknown(index).Size),[2 1 3:length(info.Unknown(index).Size)]);
			  
			  if ~raw & ~isempty(nullvals)   then
			    data(data==nullvals) = NaN;
			  end
			  //Scale data
			  if ~raw   then
			    data = double(data)*bscale+bzero;
			  end
			end
endfunction   //END READUNKNOWN


function sci_data = precision2scilab(precision)
    sci_data = 'cb';
	 select precision,
	  case 'uint8' then  sci_data = 'cb',
	  case 'int16' then  sci_data = 'sb',    //ib
	  case 'int32' then  sci_data = 'lb',
	  case 'single' then  sci_data = 'fb',
	  case 'double' then  sci_data = 'db';
	  
	  case 'char' then  sci_data = 'ucb',   
	  case 'bit16' then  sci_data = 'usb',
	  case 'single complex' then  sci_data = 'fb',
	  case 'double complex' then  sci_data = 'db';	  
	  else return,
	 end
endfunction    //precision2scilab