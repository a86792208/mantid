#ifndef MANTID_ALGORITHMS_REBINBYPULSETIMES_H_
#define MANTID_ALGORITHMS_REBINBYPULSETIMES_H_

#include "MantidKernel/System.h"
#include "MantidAPI/Algorithm.h"

namespace Mantid
{
namespace Algorithms
{

  /** RebinByPulseTimes : Rebin an input EventWorkspace according to the pulse times of the events.
    
    Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

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

    File change history is stored at: <https://github.com/mantidproject/mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
  */
  class DLLExport RebinByPulseTimes  : public API::Algorithm
  {
  public:
    RebinByPulseTimes();
    virtual ~RebinByPulseTimes();
    
    virtual const std::string name() const;
    ///Summary of algorithms purpose
    virtual const std::string summary() const {return "Bins events according to pulse time. Binning parameters are specified relative to the start of the run.";}

    virtual int version() const;
    virtual const std::string category() const;

  private:

    void init();
    void exec();
  };


} // namespace Algorithms
} // namespace Mantid

#endif  /* MANTID_ALGORITHMS_REBINBYPULSETIMES_H_ */