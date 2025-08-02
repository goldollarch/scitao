#ifndef SVD_H
#define SVD_H

#include "linear_algebra.h"

namespace Arroyo {

  ////////////////////////////////////
  // Template class to compute and store
  // the singular value decomposition of 
  // a matrix
  template<class precision>
    class svd :
  public pixel_array<precision> {
    
  private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("singular value decomposition"));};

  protected:

    precision *singular_values;
    precision *g_gtranspose_eigenmode_data;
    precision *gtranspose_g_eigenmode_data;

  public:
    
    static int verbose_level;

    svd(){};

    svd(const pixel_array<precision> & matrix){

      try{
	// sanity checks
	if(matrix.get_axes().size()!=2 ||
	   matrix.get_axes()[0] == 0 ||
	   matrix.get_axes()[1] == 0){
	  std::cerr << "svd::svd error - bad argument\n";
	  throw(std::string("svd::svd"));
	}
	
	singular_values = g_gtranspose_eigenmode_data = gtranspose_g_eigenmode_data = NULL;

	char jobu;
	char jobvt;
	int info;
	jobu = jobvt = 'A';
	
	this->pixel_array<precision>::operator=(matrix);
	int lwork = -1;

	precision optimal_workspace_size;
	
	int axes_0 = this->pixel_array<precision>::get_axes()[0];
	int axes_1 = this->pixel_array<precision>::get_axes()[1];

	if(svd::verbose_level)
	  std::cout << "svd::svd - computing optimal array size\n";
	
	singular_value_decomposition<precision>(&jobu, 
							&jobvt, 
							axes_0,
							axes_1, 
							(precision *)this->pixeldata, 
							axes_0, 
							(precision *)singular_values, 
							(precision *)g_gtranspose_eigenmode_data, 
							axes_0, 
							(precision *)gtranspose_g_eigenmode_data, 
							axes_1,
							&optimal_workspace_size, 
							lwork, 
							info);

	precision * workspace = new precision[(int)optimal_workspace_size];
	
	try{
	  this->singular_values = 
	    new precision[this->get_axes()[1]];

	    this->g_gtranspose_eigenmode_data = 
	      new precision[this->get_axes()[0]*this->get_axes()[0]];

	  this->gtranspose_g_eigenmode_data = 
	    new precision[this->get_axes()[1]*this->get_axes()[1]];
	  workspace = new precision[(int)(optimal_workspace_size)];
	} catch(...) {
	  cerr << "svd::svd error - "
	       << "unable to allocate memory to perform the singular value decomposition\n";
	  throw(string("svd::svd"));
	}
	
	// Warning - the svd appears to overwrite the input array.  Here we make the copy
	// to preserve the geometry matrix, if this has been requested
	pixel_array<precision> bkp_geometry_matrix(*this);
	
	int ispace = (int)optimal_workspace_size; 
	
	if(svd::verbose_level)
	  std::cout << "svd::svd - computing svd\n";

	singular_value_decomposition<precision>(&jobu, 
							&jobvt, 
							axes_0,
							axes_1, 
							(precision *)this->pixeldata, 
							axes_0, 
							(precision *)singular_values, 
							(precision *)g_gtranspose_eigenmode_data, 
							axes_0, 
							(precision *)gtranspose_g_eigenmode_data, 
							axes_1,
							(precision *)workspace, 
							ispace,
							info);

	
	/*
	  matrix_transpose(this->axes,
	  this->pixeldata);
	  // Transpose g_gtranspose required?
	  precision tmp;
	  for (int i = 0; i < this->axes[0]; i++) {
	  for (int j = i + 1; j < this->axes[0]; j++) {
	  tmp = this->g_gtranspose_eigenmode_data[i*this->axes[1]+j];
	  this->g_gtranspose_eigenmode_data[i*this->axes[0]+j] = 
	  this->g_gtranspose_eigenmode_data[j*this->axes[0]+i];
	  this->g_gtranspose_eigenmode_data[j*this->axes[0]+i] = tmp;
	  }
	  }
	*/


	if(svd::verbose_level)
	  std::cout << "svd::svd - cleaning up\n";

	delete [] workspace;
	
	if(svd::verbose_level)
	  std::cout << "svd::svd - cleaning up\n";

	this->pixel_array<precision>::operator=(bkp_geometry_matrix);
      } catch(...) {
	std::cerr << "svd::svd - error constructing svd\n";
	throw(std::string("svd::svd"));
      }
    };

    svd(const svd & s){
      this->operator=(s);
    };

    svd & operator=(const svd & s){
      if(this==&s)
	return(*this);

      
      int ggt_size = s.get_axes()[0]*s.get_axes()[1];
      
      int gtg_size = 
	s.get_axes()[1]*s.get_axes()[1];
      
      if(this->pixel_array<precision>::get_axes() !=
	 s.get_axes()){

	delete [] this->singular_values;
	delete [] this->g_gtranspose_eigenmode_data;
	delete [] this->gtranspose_g_eigenmode_data;

	if(s.get_axes()[0]!=0){
	  try{
	    this->singular_values = new precision[s.get_axes()[1]];
	    this->g_gtranspose_eigenmode_data = new precision[ggt_size];
	    this->gtranspose_g_eigenmode_data = new precision[gtg_size];
	  } catch(...) {
	    std::cerr << "svd::operator= - error allocating memory\n";
	    throw(std::string("svd::operator="));
	  }
	}
      }

      memcpy(s.singular_values, this->singular_values, s.get_axes()[1]*sizeof(precision));
      memcpy(s.singular_values, this->g_gtranspose_eigenmode_data, ggt_size*sizeof(precision));
      memcpy(s.singular_values, this->gtranspose_g_eigenmode_data, gtg_size*sizeof(precision));
      this->pixel_array<precision>::operator=(s);
      
      return(*this);
    };

    ~svd(){
      delete [] this->singular_values;
      delete [] this->g_gtranspose_eigenmode_data;
      delete [] this->gtranspose_g_eigenmode_data;
    };

    void read(const char * filename){
      iofits iof;
      try{iof.open(filename);}
      catch(...){
	cerr << "svd::read - "
	     << "error opening file " << filename << endl;
	throw(string("svd::read"));
      }
      try{this->read(iof);}
      catch(...){
	cerr << "svd::read - "
	     << "error reading "
	     << this->unique_name() << " from file " 
	     << filename << endl;
	throw(string("svd::read"));
      }
    };

    void read(const iofits & iof) {
      if(!iof.key_exists("TYPE")){
	cerr << "svd::read error - "
	     << "unrecognized type of file\n";
	throw(string("svd::read"));
      }
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "svd::read error - file of type " 
	     << type << " rather than of type "
	     << this->unique_name() << endl;
	throw(string("svd::read"));
      }
      
      try{
	this->pixel_array<precision>::read(iof);
      } catch(...) {
	std::cerr << "svd::read - error reading matrix from iofits\n";
	throw(std::string("svd::read"));
      }

      try{
	if(this->singular_values!=NULL)
	  delete [] this->singular_values;
	this->singular_values = new precision[this->get_axes()[1]];
	iof.read_image(0, this->get_axes()[1]-1, this->singular_values);
      } catch(...) {
	std::cerr << "svd::read - error reading eigenvalues from iofits\n";
	throw(std::string("svd::read"));
      }

      
      iof.movrel_hdu(1);
      
      int nelem = this->get_axes()[0]*this->get_axes()[0];
      
      try{
	if(this->g_gtranspose_eigenmode_data!=NULL)
	  delete [] this->g_gtranspose_eigenmode_data;
	this->g_gtranspose_eigenmode_data = 
	  new precision[nelem];
	iof.read_image(0, nelem-1, this->g_gtranspose_eigenmode_data);
      } catch(...) {
	std::cerr << "svd::read - error reading AAt from iofits\n";
	throw(std::string("svd::read"));
      }
      
      iof.movrel_hdu(1);
      
      nelem = this->get_axes()[1]*this->get_axes()[1];
      
      try{
	if(this->gtranspose_g_eigenmode_data!=NULL)
	  delete [] this->gtranspose_g_eigenmode_data;
	this->gtranspose_g_eigenmode_data = 
	  new precision[nelem];
	iof.read_image(0, nelem-1, this->gtranspose_g_eigenmode_data);
      } catch(...) {
	std::cerr << "svd::read - error reading AtA from iofits\n";
	throw(std::string("svd::read"));
      }

      // leave iof pointing to the next header, if it exists
      if(iof.get_hdu_num()!=iof.get_num_hdus())
	iof.movrel_hdu(1);
    };

    void write(const char * filename) const {

      iofits iof;
      try{iof.create(filename);}
      catch(...){
	cerr << "svd::write - "
	     << "error opening file " << filename << endl;
	throw(string("svd::write"));
      }
      
      try{this->write(iof);}
      catch(...){
	cerr << "svd::write - "
	     << "error writing "
	     << this->unique_name() << " to file "
	     << filename << endl;
	throw(string("svd::write"));
      }
    };

    void write(iofits & iof) const {

      fits_header_data<precision> fhd(this->axes);
      fhd.write(iof);
      string type = this->unique_name();
      string comment = "object type";
      iof.write_key("TYPE", type, comment);
      this->pixel_array<precision>::write(iof);

      // The eigenvalues
      fhd = fits_header_data<precision> (std::vector<long>(1,this->axes[1]));
      fhd.write(iof);
      iof.write_image(0, this->get_axes()[1]-1, this->singular_values);

      vector<long> ggt_axes(2,this->get_axes()[0]);
      int nelem = ggt_axes[0]*ggt_axes[1];
      
      // The eigenmodes of AAt
      fhd = fits_header_data<precision> (std::vector<long>(2,this->axes[0]));
      fhd.write(iof);
      iof.write_image(0, nelem-1, this->g_gtranspose_eigenmode_data);
      
      vector<long> gtg_axes(2,this->get_axes()[1]);
      nelem = gtg_axes[0]*gtg_axes[1];
      
      // The eigenmodes of AtA
      fhd = fits_header_data<precision> (std::vector<long>(2,this->axes[1]));
      fhd.write(iof);
      iof.write_image(0, nelem-1, this->gtranspose_g_eigenmode_data);
      
    };

    ///////////////////////////////////////////
    ///  Get the number of eigenmodes in the AAt matrix
    ///
    ///  If the SVD was performed without saving the
    ///  eigenmodes, this function throws an error.
    int get_n_AAt_eigenmodes() const {
      return this->axes[0];
    };

    ///////////////////////////////////////////
    ///  Get the number of eigenmodes in the AtA matrix
    ///
    ///  If the SVD was performed without saving the
    ///  eigenmodes, this function throws an error.
    int get_n_AtA_eigenmodes() const {
      return this->axes[1];
    };

    ///////////////////////////////////////////
    ///  Get singular value matrix
    pixel_array<precision> get_singular_value_matrix() const {
      
      int min_axis = this->axes[0]>this->axes[1] ? this->axes[1] : this->axes[0];

      pixel_array<precision> singular_value_matrix(std::vector<long>(2,min_axis));

      for(int i=0; i<min_axis; i++){
	//std::cout  << "min axis " << min_axis << "\ti " << i << "\t" << this->singular_values[i] << std::endl;
	singular_value_matrix.set_data(i*min_axis+i,this->singular_values[i]);
      }

      return singular_value_matrix;
    };

    ///////////////////////////////////////////
    ///  Get inverse singular value matrix
    pixel_array<precision> get_inverse_singular_value_matrix() const {
      
      int min_axis = this->axes[0]>this->axes[1] ? this->axes[1] : this->axes[0];

      pixel_array<precision> singular_value_matrix(std::vector<long>(2,min_axis));

      for(int i=0; i<min_axis; i++){
	//std::cout  << "min axis " << min_axis << "\ti " << i << "\t" << this->singular_values[i] << std::endl;
	singular_value_matrix.set_data(i*min_axis+i,1/this->singular_values[i]);
      }

      return singular_value_matrix;
    };

    ///////////////////////////////////////////
    ///  Get condition number
    double get_condition_number() const {
      int min_axis = this->axes[0]>this->axes[1] ? this->axes[1] : this->axes[0];
      return(this->singular_values[0]/this->singular_values[min_axis-1]);
    }; 


    ///////////////////////////////////////////
    ///  Retrieve the eigenvalue for a particular mode
    ///
    ///  If the SVD was performed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.  Note - number of eigenvalues is equal to the dimension
    ///  of AAt
    precision get_eigenvalue(int mode_number) const {
      if(mode_number>=this->axes[0]){
	std::cerr << "svd::get_eigenvalue error - mode number " 
		  << mode_number
		  << " out of range 0 - "
		  << this->axes[0]
		  << std::endl;
	throw(std::string("svd::get_eigenvalue"));
      }
      return this->singular_values[mode_number];
    };

    ///////////////////////////////////////////
    ///  Retrieve AAt
    pixel_array<precision> get_U() const {
      return(pixel_array<precision>(std::vector<long>(2,this->axes[0]),
					    this->g_gtranspose_eigenmode_data));
    };

    ///////////////////////////////////////////
    ///  Retrieve AAt range
    pixel_array<precision> get_range_U(double thresh = -1) const {
      
      int dimen = 0;
      for(int i=0; i<this->axes[1]; i++)
	if(this->singular_values[i]> thresh) dimen++;	  

      std::vector<long> axes(2,this->axes[0]);
      axes[1] = dimen;

      pixel_array<precision> range_U(axes);

      for(int i=0; i<dimen; i++)
	for(int j=0; j<this->axes[0]; j++)
	  range_U.set_data(j*dimen+i,
			   this->g_gtranspose_eigenmode_data[j*dimen+i]);

      return(range_U);
    };

    ///////////////////////////////////////////
    ///  Retrieve AAt nullspace
    pixel_array<precision> get_null_U(double thresh = -1) const {
      
      int dimen = this->axes[0];
      for(int i=0; i<this->axes[1]; i++)
	if(this->singular_values[i]> thresh) dimen--;	  

      std::vector<long> axes(2,this->axes[0]);
      axes[1] = dimen;

      pixel_array<precision> null_U(axes);

      int index;
      for(int i=(this->axes[0]-dimen); i<this->axes[0]; i++){

	index = i-(this->axes[0]-dimen);

	for(int j=0; j<this->axes[0]; j++)
	  null_U.set_data(index * this->axes[0] + j, 
			  this->g_gtranspose_eigenmode_data[i*this->axes[0] + j]);
      }
      return(null_U);
    };

    ///////////////////////////////////////////
    ///  Retrieve an AAt eigenmode.  This function returns the
    ///  eigenvalue of the mode.
    ///
    ///  If the SVD was performed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.
    pixel_array<precision> 
      get_AAt_eigenmode(int mode_number) const {

      pixel_array<precision> pixarr(std::vector<long>(1,this->axes[0]));

      if(mode_number>=this->axes[0]){
	std::cerr << "svd::get_AAt_eigenmode error - mode number " 
		  << mode_number
		  << " out of range 0 - "
		  << this->axes[0]
		  << std::endl;
	throw(std::string("svd::get_AAt_eigenmode"));
      }

      for(int i=0; i<this->axes[0]; i++)
	pixarr.set_data(i, this->g_gtranspose_eigenmode_data[mode_number*this->axes[0]+i]);
      return(pixarr);
    };

    ///////////////////////////////////////////
    ///  Retrieve AtA
    pixel_array<precision> get_Vt() const {
      return(pixel_array<precision>(std::vector<long>(2,this->axes[1]),
					    this->gtranspose_g_eigenmode_data));
    };

    ///////////////////////////////////////////
    ///  Retrieve an AtA eigenmode.  This function returns the
    ///  eigenvalue of the mode.
    ///
    ///  If the SVD was performed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.
    pixel_array<precision> get_AtA_eigenmode(int mode_number) const {

      pixel_array<precision> pixarr(std::vector<long>(1,this->axes[1]));

      if(mode_number<0 || mode_number>=this->axes[1]){
	std::cerr << "svd::get_AtA_eigenmode error - mode number " 
		  << mode_number
		  << " out of range 0 - "
		  << this->axes[1]
		  << std::endl;
	throw(std::string("svd::get_AtA_eigenmode"));
      }

      for(int i=0; i<this->axes[1]; i++)
	pixarr.set_data(i, this->gtranspose_g_eigenmode_data[mode_number*this->axes[1]+i]);
      return(pixarr);
    };

    ///////////////////////////////////////////
    ///  Retrieve the pseudoinverse
    ///
    ///  Eigenvalue threshold is a fraction of the 
    ///  largest eigenvalue
    pixel_array<precision> get_pseudoinverse(double eigenvalue_threshold=0) const {

      std::vector<long> pseudoinverse_axes(2,this->axes[1]);
      pseudoinverse_axes[1] = this->axes[0];

      pixel_array<precision> pseudoinverse(pseudoinverse_axes);

      // The eigenvalues are sorted in descending order.  We only process
      // entries for which the eigenvalues are above the threshold
      int limit = 0;
      double scaled_singular_value_threshold = 
	sqrt(eigenvalue_threshold*this->singular_values[0]*this->singular_values[0]);
      while(limit<this->axes[1] && 
	    this->singular_values[limit]>scaled_singular_value_threshold) 
	limit++;
      
      precision tmp;
      int indexa;
      for(int i=0; i<this->axes[1]; i++){
	for(int j=0; j<this->axes[0]; j++){
	  tmp = 0;

	  indexa = i*this->axes[1];
	  //indexa = i*this->axes[0];

	  for(int k=0; k<limit; k++)
	    tmp += this->gtranspose_g_eigenmode_data[indexa+k]*
	      this->g_gtranspose_eigenmode_data[j+k*this->axes[0]]/this->singular_values[k];
	  //pseudoinverse.set_data(i*this->axes[0]+j,tmp);
	  pseudoinverse.set_data(j*this->axes[1]+i,tmp);
	}  
      }
      return pseudoinverse;
    };
  };

  template<class T>
    int svd<T>::verbose_level = 0;

};

#endif
