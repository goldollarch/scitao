/*
Arroyo - software for the simulation of electromagnetic wave propagation
through turbulence and optics.

Copyright (c) 2000-2004 California Institute of Technology.  Written by
Dr. Matthew Britton.  For comments or questions about this software,
please contact the author at mbritton@astro.caltech.edu.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  In no
event shall California Institute of Technology be liable to any party
for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if the California Institute of Technology has been
advised of the possibility of such damage.   The California Institute of
Technology has no obligation to provide maintenance, support, updates,
enhancements or modifications.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef DETECTOR_H
#define DETECTOR_H

namespace Arroyo {

  using std::vector;

  ///  A class to represent a detector
  ///  the way this detector is intended
  ///  to be used is by detecting a single
  ///  wavefront from a point emitter repeatedly.
  ///  When you want to read out the detection,
  ///  first call convolve on the extended_emitter 
  ///  if you want to detect an extended source.  
  ///  Then call readout, which will add the pixel
  ///  gain fluctuations, read noise, and bias,
  ///  if these are present.

  class detector :
    public three_frame, 
    public electromagnetic_spectrum {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    detector(){};

    ///////////////////////////////////////////
    ///  Construct from iofits object
    detector(const iofits & iof);

    ///////////////////////////////////////////
    ///  Destructor
    virtual ~detector(){};

    ///////////////////////////////////////////
    ///  Operator = 
    detector & operator=(const detector & dtcr);

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from iofits object
    virtual void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Convolve with an emitter
    /// virtual void convolve(const extended_emitter & xtnd_emtr) = 0;

    ///////////////////////////////////////////
    ///  Detect a wavefront - adds wavefront
    ///  into the detector's pixels.
    virtual void detect(const wavefront & wf) = 0;

    ///////////////////////////////////////////
    ///  Read out detector
    ///  virtual simulated_AO_observation readout() const = 0;

    ///////////////////////////////////////////
    ///  Clear detector
    virtual void clear() = 0;

    ///////////////////////////////////////////
    ///  Factory constructor from file
    static detector * detector_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory constructor from iofits object
    static detector * detector_factory(const iofits & iof);

  };

  ///
  /// A class to represent an infrared detector
  ///
  /// According to Roger Smith (4/12/05):
  /// 
  /// Detector response may be modeled as
  ///
  /// (output) = C +
  ///            (bias) +
  ///            (dark)*t +
  ///            (gain)*(signal)*(time)
  ///
  /// where C is the detector read noise,
  /// (bias)+(dark)*t are detector contributions
  /// and (gain)*(signal)*(time) is the actual
  /// response due to the signal.  The latter three
  /// terms generate poisson noise contributions.
  ///
  /// The first term is composed of a number of 
  /// contributions:
  ///   kTC noise (multiplied by the gain)
  ///   Bias voltage noise coupling  (multiplied by the gain)
  ///   RF interference
  ///   ADC noise
  ///   Preamp noise  (multiplied by the gain)
  ///   quantization noise (roundoff in digitization)
  ///
  /// The noise process describing the first term is not
  /// well determined without an understanding of all these
  /// effects.  However, to first order it can be assumed 
  /// Gaussian for any single pixel in the array.  Over the
  /// whole array the distribution appears sort of Poisson-like.
  ///
  /// When detector engineers quote something like "7 electrons read
  /// noise" they are really saying that this is the residual noise
  /// after correlated double sampling.
  ///
  /// There is another issue here in that there is a conversion 
  /// from ADU to electrons, which Roger finds by comparing mean
  /// to variance in images with different signal levels.  Some
  /// of the above terms are quoted in electrons, and others in
  /// ADU.

  class infrared_detector :
    public detector {

    protected:
  
    /// The pixel dimensions of the detector
    vector<long> axes;

    /// The pixel scale
    double pixel_scale;

    /// The read noise, in electrons
    double read_noise;

    /// The detector bias, pixel by pixel
    float * pixel_bias;

    /// The detector dark current, pixel by pixel
    float * pixel_dark_current;

    /// The detector gains, pixel by pixel
    float * pixel_gains;

    /// The ADU's detected by each pixel
    float * pixel_counts;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    infrared_detector();

    ///////////////////////////////////////////
    ///  Copy constructor
    infrared_detector(const infrared_detector & dtcr);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    infrared_detector(iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    infrared_detector(const char * filename);

    ///////////////////////////////////////////
    ///  Destructor
    ~infrared_detector();

    ///////////////////////////////////////////
    ///  Operator = 
    infrared_detector & operator=(const infrared_detector & dctr);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Convolve with an emitter
    ///void convolve(const extended_emitter & xtnd_emtr);

    ///////////////////////////////////////////
    ///  Detect a wavefront
    void detect(const wavefront & wf);

    ///////////////////////////////////////////
    ///  Read out infrared_detector
    ///simulated_AO_observation readout() const;

    ///////////////////////////////////////////
    ///  Clear infrared_detector
    void clear();

    ///////////////////////////////////////////
    ///  Realize random pixel gains with 
    ///  mean gain <mean> and rms gain <rms>
    ///  resulting gains have a gaussian distribution
    void randomize_pixel_gains();
 
  };

}

#endif

