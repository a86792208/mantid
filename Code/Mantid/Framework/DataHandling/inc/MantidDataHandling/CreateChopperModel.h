#ifndef MANTID_DATAHANDLING_CREATECHOPPERMODEL_H_
#define MANTID_DATAHANDLING_CREATECHOPPERMODEL_H_

#include "MantidAPI/Algorithm.h"

namespace Mantid
{
  namespace DataHandling
  {

    /**
      Defines an algorithm to set chopper properties on a given workspace. The
      properties are mostly specifed by string to allow flexibility for future
      models with different parameters

      Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge National Laboratory & European Spallation Source

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
    class DLLExport CreateChopperModel  : public API::Algorithm
    {
    public:

      virtual const std::string name() const;
    ///Summary of algorithms purpose
    virtual const std::string summary() const {return "Creates a chopper model for a given workspace";}

      virtual int version() const;
      virtual const std::string category() const;

    private:
      void init();
      void exec();
    };


  } // namespace DataHandling
} // namespace Mantid

#endif  /* MANTID_DATAHANDLING_CREATECHOPPERMODEL_H_ */
