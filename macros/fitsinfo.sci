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

function info = fitsinfo(filename)
 info =[];
	//Open file 
	[fid,err] = mopen(filename,'rb');
	if fid == -1 then
	   error('Unable to open this file.');
	end
 	[d,ierr] = fileinfo(filename);  
 	        x= getdate(d(6));
 	        datecc= [x(1) x(2) x(6) x(7) x(8)];

    info.Filename=filename;
    info.FileModDate=datecc;      	 
    info.FileSize=d(1);
    info.Contents=[];
    info.PrimaryData=[];
    info.AsciiTable=list();      	 
    info.BinaryTable=list();
    info.Image=list();
    info.Unknown=list();

	[info.PrimaryData,headertype,datasize,atEOF,extensions] = readheaderunit(fid);
	
		if 'SIMPLE' <> info.PrimaryData.Keywords(2)(1)  then
		  warning('File is not a standard FITS file.  \nFITSINFO will attempt to determine the contents.');
		else
		  info.Contents = 'Primary  ';
		end	
		
		asciicount = 1;
		binarycount = 1;
		imagecount = 1;
		unknowncount = 1;
		
			if extensions  then
			  while ~atEOF
			    //Data must be multiple of 2880 bytes
			    if datasize <> 0   then
			      bseek = 2880 - modulo(datasize,2880)+ datasize;
			      mseek(bseek,fid,'cur');
			    end     //if
			    [extensioninfo,headertype,datasize,atEOF,extensions] = readheaderunit(fid);
			    if ~length(extensioninfo) then
			      mclose(fid);
			      break;
			    end
			    
			    select headertype,
			     case 'ascii'   then
			      info.AsciiTable = extensioninfo; 
			      info.Contents = info.Contents + 'ASCII Table  ';
			      asciicount = asciicount+1;
			     case 'binary'  then
			      info.BinaryTable(binarycount) = extensioninfo;
			      info.Contents = info.Contents + 'Binary Table  ';
			      info.BinaryTable(binarycount).ExtensionOffset = info.BinaryTable(binarycount).Offset+ info.BinaryTable(binarycount).RowSize*info.BinaryTable(binarycount).Rows;
			      binarycount = binarycount+1;
			     case 'image'   then
			      info.Image(imagecount) = extensioninfo;
			      info.Contents = info.Contents + 'Image  ';
			      imagecount = imagecount+1;
			     case 'unknown' then
			      info.Unknown(unknowncount) = extensioninfo;
			      info.Contents = info.Contents + 'Unknown  ';
			      unknowncount = unknowncount+1;
			    end  //select
			  end   //while
			end    //if

	
	//Close file
  mclose(fid);
endfunction


function [info,headertype,datasize,atEOF,extend] = readheaderunit(fid)
	 extend = 0.; 
	 atEOF = 0.;  
	 datasize = 5.;
	 headertype = '';
	 info = list();
	 extension = 0.;
	 
	 endFound = 0.;

	card = mgeti(80,'uc',fid);
	readcount = length(card);
		if meof(fid) then
	      return;
	  elseif readcount <> 80  then
	        mclose(fid);
	        error('FITS card image was not 80 bytes. FITS file may be invalid or corrupt.');
	  end
	  //end 
   
   card = char(card);
			try
			  [keyword,value,comment] = parsecard(card);
			catch
			  mclose(fid);
			  error('Unexpected FITS card image format. FITS file may be invalid or corrupt.');
			end  
			
			select convstr(keyword,'l'), 
			 case 'simple'    then
			  headertype = 'primary';
			  info = tlist(['pri','DataType','Size','DataSize','MissingDataValue','Intercept','Slope','Offset','Keywords'], ..
			               "",[],0.,[],[],[],[],list(["keywords","value","comment"]));
			 case 'xtension'   then
					  value_tt= msscanf(value,'%s');     //valu = msscanf(value,'%f')
					  select convstr(value_tt,'l'), 
					   case 'table'   then
					    headertype = 'ascii'; 
					    info = tlist(['asc','Rows','RowSize','NFields','FieldFormat','FieldPrecision','FieldWidth','FieldPos','DataSize','MissingDataValue','Intercept','Slope','Offset','Keywords'], ..
					                 [],[],[],list(),list(),[],[],[],list(),[],[],[],list(["keywords","value","comment"]));
					   case 'bintable'  then
					    headertype = 'binary';
					    info = tlist(['bin','Rows','RowSize','NFields','FieldFormat','FieldPrecision','FieldSize','DataSize','MissingDataValue','Intercept','Slope','Offset','ExtensionSize','ExtensionOffset','Keywords'], ..
					                  [],[],[],list(),list(),[],[],[],[],[],[],[],[],list(["keywords","value","comment"]));
					   case 'image'    then
					    headertype = 'image';
					    info = tlist(['ima','DataType','Size','DataSize','Offset','MissingDataValue','Intercept','Slope','Keywords'], ..
					                 "",[],0.,[],[],[],[],list(["keywords","value","comment"]));
					   else     
					    headertype = 'unknown';
					    info = tlist(['unk','DataType','Size','DataSize','PCOUNT','GCOUNT','Offset','MissingDataValue','Intercept','Slope','Keywords'], ..
					                 "",[],0.,[],[],[],[],[],[],list(["keywords","value","comment"]));
					  end  //select
			end  //select

    info.Keywords($+1)=[keyword,value,comment];		
    while ~endFound    //~atEOF keyword <> "END"
       card = char(mgeti(80,'uc',fid));
       try
	    		[keyword,value,comment] = parsecard(card);
	     catch
	    		mclose(fid);
	    		error('Unexpected FITS card image format. FITS file may be invalid or corrupt.');
	  	 end
	  	 
			  select headertype,
			   case 'primary'   then     
			    [info,endFound,extension] = knownprimarykeywords(info,keyword,value);
			   case 'ascii'     then
			    [info,endFound] = knowntablekeywords(info,keyword,value);
			   case 'binary'    then
			    [info,endFound] = knownbintablekeywords(info,keyword,value);
			   case 'image'     then
			    [info,endFound] = knownimagekeywords(info,keyword,value);
			   else
			    [info,endFound] = knowngenerickeywords(info,keyword,value);
			  end	 //select
			  
			    if keyword <> "" then
               info.Keywords($+1)=[keyword,value,comment];	
          end     
				  if extension   then
				    extend = 1.;
				  end
				  atEOF = meof(fid);
    end	 //while

			info.DataSize = 20.; 
			select headertype,
			 case 'primary'   then
				  if ~isempty(info.Size)    then
				    info.DataSize = prod(info.Size)*precision2bitdepth(info.DataType)/8;
				  end
			 case 'image'     then
			  info.DataSize = prod(info.Size)*precision2bitdepth(info.DataType)/8;
			 case 'ascii'	    then
			  info.DataSize = info.RowSize*info.Rows;
			 case 'binary'	  then
			  info.DataSize = info.RowSize*info.Rows;
			 case 'unknown'	  then
			  if ~isempty(info.Size)     then
			    info.DataSize = (precision2bitdepth(info.DataType)*info.GCOUNT*(info.PCOUNT + prod(info.Size)))/8;
			  end
			end	
					
				if headertype == "binary"   then
			  datasize = info.DataSize+info.ExtensionSize;
			else
			  datasize = info.DataSize;
			end	
			
			curr_pos = mtell(fid);
			if curr_pos <= (36*80)   then 
			    info.Offset = curr_pos;
			    mseek(2880,fid,'set');			    
			else 
					if modulo(curr_pos,80) <> 0  then
					    warning('The header did not end with a full card! Check the file!');
					    info.Offset = curr_pos;  
					else
					    info.Offset = curr_pos; 
					end 

			      bsek = 2880 - modulo(curr_pos,2880) + curr_pos;
			      mseek(bsek,fid,'set');			
			end			
//  mclose(fid);		   be careful!
  	
endfunction	



function [keyword, value, comment] = parsecard(card)
   
  ca = strsplit(card,8);
  keyword = stripblanks(ca(1));
  card = ca(2);
  
  if keyword == '' then
      value = 0;
      comment="";
      return;
  end    
  
  if keyword == 'COMMENT'  then
     value = "";
      comment = card;
     return;
  end
  
  // Seperate Value / Comment
  slashidx = strindex(card,'/');
  tempcard = strsubst(card,'''''','''');
  
  quoteidx = strindex(tempcard,'''');
		if isempty(slashidx) then
		  value = card;
		  comment = '';
		else
			  if length(quoteidx)>1 then
			  
			        if slashidx(1) > quoteidx(1) & slashidx(1) < quoteidx(2)  then
			           slashidx = slashidx(slashidx>quoteidx(2));
			               if isempty(slashidx) then
			                     value = card;
			                     comment = '';
			               else
			                     sl = strsplit(card,(slashidx(1)-1));
			                     value = sl(1);
			                     comment = sl(2);
			               end
			        else
			            sl = strsplit(card,(slashidx(1)-1));
			            value = sl(1);
			            comment = sl(2);
			        end
			  
			  else
			    sl = strsplit(card,(slashidx(1)-1));
			    value = sl(1);
			    comment = sl(2);
	  		end
    end 

   if ~isempty(comment) & length(comment)>=2 then
       co = strsplit(comment,2);
       comment = co(2); 
   end
   
   //value
   va = strsplit(value,2);
   if va(1) == '= ' then
        value = va(2);
         value = stripblanks(value);
        if msscanf(value,'%c') == '''' then
            value = stripblanks(value);
            val = strsplit(value,[1,length(value)-1]);            
            if val(3) == '''' then
                 value = val(2);
            else 
            end    //Complex number in form of (spaces intentional) is not consider. gan
        else 
       end
   
   else  
      value = stripblanks(value);
   end

endfunction   //parsecard

function [info,endFound,extensions] = knownprimarykeywords(info,keyword,value)
		endFound = 0.;
		extensions = 0.;
		
		if keyword == "" | keyword =="COMMENT"  then
		     info = info;
		     return;
		end     
		
		valu = msscanf(value,'%f'); 
		lit_key = convstr(keyword,'l');		       		
		select lit_key,
		 case 'bitpix'    then		  info.DataType =  bitdepth2precision(valu);
		 case 'naxis'     then		  info.Intercept = 0.,		  info.Slope = 1.;
		 case 'extend'    then		  extensions = 1.;
		 case 'bzero'     then		  info.Intercept = valu;
		 case 'bscale'    then		  info.Slope = valu;
		 case 'blank'     then		  info.MissingDataValue = valu;
		 case 'end'       then		  endFound = 1.;    
		 else   
			  //Take care of NAXISN keywords
			  if grep(lit_key,'naxis')   then 
			    ss_dim = strsplit(lit_key,5); 
			    dim = msscanf(ss_dim(2),'%f'); 
				    if dim > 0 then   
				       info.Size(dim) = valu;
				    end   
			  end 		 // return;
   end   // select
endfunction    //knownprimarykeywords

function [info,endFound] = knowntablekeywords(info,keyword,value)
		endFound = 0;
		
		if keyword == "" | keyword =="COMMENT"  then
		     info = info;
		     return;
		end     

		valu = msscanf(value,'%f'); 
		lit_key = convstr(keyword,'l');				
		select lit_key, 
		 case 'end'     then		  endFound = 1;
		 case 'naxis1'  then		  info.RowSize = valu;
		 case 'naxis2'  then		  info.Rows = valu;
		 case 'tfields' then
		  info.NFields = valu;
		  info.Slope = ones(1,valu);
		  info.Intercept = zeros(1,valu);
		  info.MissingDataValue = zeros(1,valu); 
		 else
				  //Take care of indexed keywords
				  if grep(lit_key,'tform')   then       //No1
				    ss_dim=strsubst(keyword,'tform',''),    
				    idx = msscanf(ss_dim,'%s');
				    
				    select(msscanf(value,'%c'))
				     case 'A'    then				      formats = 'Char';
				     case 'I'    then				      formats = 'Integer';
				     case 'E'    then             formats = 'Single';
				     case 'F'    then				      formats = 'Single';
				     case 'D'    then				      formats = 'Double';
				     else           				      formats = 'Unknown'
				    end
				    
				    info.FieldFormat($+1) = value;
				    info.FieldPrecision($+1) = formats;
				    width = msscanf(value,' %*c%f');
				    info.FieldWidth(idx) = width;
				  else
				  if grep(convstr(keyword,'l'),'tbcol')   then     //No2
				    ss_dim=strsubst(keyword,'tbcol',''),  
				    idx = msscanf(ss_dim,'%s');
				    info.FieldPos(idx) = valu;
				  else
				  if grep(convstr(keyword,'l'), 'tscal')   then     //No3
				    ss_dim=strsubst(keyword,'tscal',''),  
				    tscale_idx = msscanf(ss_dim,'%i');
				    info.Slope(tscale_idx) = valu;
				  else
				  if grep(convstr(keyword,'l'), 'tzero')   then     //No4
				    ss_dim=strsubst(keyword,'tzero',''),  
				    tzero_idx = msscanf(ss_dim,'%i');
				    info.Intercept(tzero_idx) = valu;
				  else
				  if grep(convstr(keyword,'l'), 'tnull')    then     //No5
				    ss_dim=strsubst(keyword,'tnull',''),  
				    tnull_idx = msscanf(ss_dim,'%i');
				    info.MissingDataValue(tnull_idx) = valu;
				  end   //No5
				  end   //No4
				  end   //No3
				  end   //No2
				  end   //No1
		end    //select
endfunction    //knowntablekeywords

function [info,endFound] = knownbintablekeywords(info,keyword,value)
	  endFound = 0;
		
		if keyword == "" | keyword =="COMMENT"  then
		     info = info;
		     return;
		end     

		valu = msscanf(value,'%f'); 
		lit_key = convstr(keyword,'l');						
			select lit_key, 
			 case 'end'      then  endFound = 1;
			 case 'naxis1'   then  info.RowSize = valu;
			 case 'naxis2'   then  info.Rows = valu;
			 case 'tfields'  then
				  info.NFields = valu;
				  info.Slope = ones(1,valu);
				  info.Intercept = zeros(1,valu);
				  info.MissingDataValue = zeros(1,valu);
			 case 'pcount'   then  info.ExtensionSize = valu;
			 else
			   if grep(lit_key, 'tscal')  then   //No1
			    ss_dim=strsplit(lit_key,5);
			    tscale_idx = msscanf(ss_dim(2),'%i'); 
			    info.Slope(tscale_idx) = valu;
			   elseif grep(lit_key, 'tzero')   then   //No2
			    ss_dim=strsplit(lit_key,5);
			    tzero_idx = msscanf(ss_dim(2),'%i');
			    info.Intercept(tzero_idx) = valu;
			   elseif grep(lit_key, 'tnull')    then   //No3
			           ss_dim=strsplit(lit_key,5);
			           tnull_idx = msscanf(ss_dim(2),'%i');
			           info.MissingDataValue(tnull_idx) = valu;
			   elseif grep(lit_key, 'tform')     then   //No4
			           ss_dim=strsplit(lit_key,5);
			           idx = sscanf(ss_dim(2),'%i');
			           repeat = msscanf(1,value,'%i');
								    if isempty(repeat)  then
								      repeat = 1;
								    end
							   formats = msscanf(1,value,'%*i%c');
								    if isempty(formats)  then
								      formats = msscanf(1,value,' %c');
								    end
			
								   select formats,  
								     case 'L'    then      formats = 'char';
								     case 'X'    then
								       if repeat <> 0   then	
								          repeat = 1;     
								       end;
								     case 'B'    then      formats = 'uint8 ';
								     case 'I'    then      formats = 'int16 ';
								     case 'J'    then      formats = 'int32 ';
								     case 'A'    then      formats = 'char ';
								     case 'E'    then      formats = 'single ';
								     case 'D'    then      formats = 'double ';
								     case 'C'    then      formats = 'single complex ';
								     case 'M'    then      formats = 'double complex ';
								     case 'P'    then      formats = 'int32 ';
								      if repeat <> 0  then
									        repeat = 2;
								      end  //if
								    end  //select
				    info.FieldFormat($+1) = value;
				    info.FieldPrecision($+1) = formats;
				    info.FieldSize(idx) = repeat;								    
			  end    //No1
		 end  //select			  
endfunction    //knownbintablekeywords

function [info,endFound] = knownimagekeywords(info,keyword,value)
			endFound = 0;
		
		if keyword == "" | keyword =="COMMENT"  then
		     info = info;
		     return;
		end     

		valu = msscanf(value,'%f'); 
		lit_key = convstr(keyword,'l');					
			select lit_key,  
			 case 'end'     then  endFound = 1;
			 case 'bitpix'  then  info.DataType =  bitdepth2precision(valu);
			 case 'naxis'   then  info.Intercept = 0,  info.Slope = 1;
			 case 'blank'   then  info.MissingDataValue= valu;
			 case 'bzero'   then  info.Intercept = valu;
			 case 'bscale'  then  info.Slope = valu;
			 else
					  //NAXISN keywords
			  if grep(lit_key,'naxis')   then 
			    ss_dim = strsplit(keyword,5); 
			    dim = msscanf(ss_dim(2),'%f'); 
			    if dim > 0 then
			       info.Size(dim) = valu;
			    end   
			  end					     
			end  
endfunction    //knownimagekeywords

function [info,endFound] = knowngenerickeywords(info,keyword,value)
		endFound = 0;
		
		if keyword == "" | keyword =="COMMENT"  then
		     info = info;
		     return;
		end     

		valu = msscanf(value,'%f'); 
		lit_key = convstr(keyword,'l');		
		select lit_key,   
		 case 'bitpix'   then    info.DataType =  bitdepth2precision(valu),
		 case 'naxis'    then    info.Intercept = 0,  info.Slope = 1,
		 case 'extend'   then    extensions = 1,
		 case 'bzero'    then    info.Intercept = valu,
		 case 'bscale'   then    info.Slope = valu,
		 case 'blank'    then    info.MissingDataValue = valu,
		 case 'pcount'   then    info.PCOUNT = valu,
		 case 'gcount'   then    info.GCOUNT = valu,
		 case 'end'      then    endFound = 1,
		 else
			  //Take care of NAXISN keywords
			  if grep(lit_key,'naxis')   then 
			    ss_dim = strsplit(keyword,5); 
			    dim = msscanf(ss_dim(2),'%f'); 
			    if dim > 0 then
			       info.Size(dim) = valu;
			    end   
			  end	 	   
		end 
endfunction    //knowngenerickeywords

function precision = bitdepth2precision(value)
		select value,
		 case 8   then  precision = 'uint8',
		 case 16  then  precision = 'int16',
		 case 32  then  precision = 'int32',
		 case -32 then  precision = 'single',
		 case -64 then  precision = 'double',
		 else return,
		end
endfunction    //bitdepth2precision

function bitdepth = precision2bitdepth(precision)
	 select precision,
	  case 'uint8' then  bitdepth = 8,
	  case 'int16' then  bitdepth = 16,
	  case 'int32' then  bitdepth = 32,
	  case 'single' then  bitdepth = 32,
	  case 'double' then  bitdepth = 64;
	  else return,
	 end
endfunction    //precision2bitdepth