#ifndef MANTID_DATAHANDLING_LOADRAWBIN0_H
#define MANTID_DATAHANDLING_LOADRAWBIN0_H


//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/Algorithm.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataHandling/LoadRawHelper.h"
#include <climits>

//----------------------------------------------------------------------
// Forward declaration
//----------------------------------------------------------------------
class ISISRAW2;

namespace Mantid
{
  namespace DataHandling
  {
    /** @class LoadRaw3 LoadRaw3.h DataHandling/LoadRawBin0.h

    Loads bin zero for all spectra from ISIS RAW file and stores it in a 2D workspace
    (Workspace2D class). LoadRawBin0 is an algorithm and as such inherits
    from the Algorithm class and overrides the init() & exec() methods.
    
    Required Properties:
    <UL>
    <LI> Filename - The name of and path to the input RAW file </LI>
    <LI> OutputWorkspace - The name of the workspace in which to store the imported data
         (a multiperiod file will store higher periods in workspaces called OutputWorkspace_PeriodNo)</LI>
    </UL>

	Optional Properties: (note that these options are not available if reading a multiperiod file)
    <UL>
    <LI> spectrum_min  - The spectrum to start loading from</LI>
    <LI> spectrum_max  - The spectrum to load to</LI>
    <LI> spectrum_list - An ArrayProperty of spectra to load</LI>
    </UL>
    
    @author Sofia Antony,ISIS,RAL
    @date 12/04/2010

    Copyright &copy; 2010 STFC Rutherford Appleton Laboratory

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>.
    Code Documentation is available at: <http://doxygen.mantidproject.org>
    */
    class DLLExport LoadRawBin0: public LoadRawHelper
    {
    public:
      /// Default constructor
      LoadRawBin0();
      /// Destructor
      ~LoadRawBin0();
      /// Algorithm's name for identification overriding a virtual method
      virtual const std::string name() const { return "LoadRawBin0"; }
      /// Algorithm's version for identification overriding a virtual method
      virtual const int version() const { return 1; }
      /// Algorithm's category for identification overriding a virtual method
      virtual const std::string category() const { return "Diagnostics"; }

    private:
      /// Overwrites Algorithm method.
      void init();
      /// Overwrites Algorithm method
      void exec();

	  /// 
      void goManagedRaw();

	  ///
	  void setOptionalProperties();
       
      /// ISISRAW class instance which does raw file reading. Shared pointer to prevent memory leak when an exception is thrown.
      boost::shared_ptr<ISISRAW2> isisRaw;
      /// The name and path of the input file
      std::string m_filename;

      /// The number of spectra in the raw file
      int m_numberOfSpectra;
      /// The number of periods in the raw file
      int m_numberOfPeriods;
	  /// number of time regime
	  int m_noTimeRegimes;

	   /// Allowed values for the cache property
      std::vector<std::string> m_cache_options;
      /// A map for storing the time regime for each spectrum
      std::map<int,int> m_specTimeRegimes;
      /// The current value of the progress counter
      double m_prog;

      /// A flag int value to indicate that the value wasn't set by users
      static const int unSetInt = INT_MAX-15;
      
      // Read in the time bin boundaries
      int m_lengthIn;
      /// boolean for list spectra options
      bool m_bmspeclist;
      /// TimeSeriesProperty<int> containing data periods.
      boost::shared_ptr<Kernel::Property> m_perioids;

	  /// total number of specs
	  int m_total_specs;
	   /// time channel vector
	  std::vector<boost::shared_ptr<MantidVec> > m_timeChannelsVec;
	};
  }
}
#endif