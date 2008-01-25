#ifndef MANTIDAPI_TRIPLE_ITERATOR_H
#define MANTIDAPI_TRIPLE_ITERATOR_H

#include <vector>
#include "MantidKernel/System.h"
#include "TripleRef.h"

namespace Mantid
{
namespace API
{
/**
 triple_iterator iterates over a workspace providing values as TripleRefs

 \class triple_iterator
 \author S. Ansell
 \date November 2007
 \version 1.0 
 
 Copyright &copy; 2007-8 STFC Rutherford Appleton Laboratories

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
template<typename _Iterator, typename _Container>
class DLLExport triple_iterator
{
private:
  ///internal workspace pointer
  _Container * const m_workspace;
  /// pointer to a TripleRef of doubles
  TripleRef<double> m_CPoint;
//  _Iterator m_CPoint;
  ///The number of times this iterator should loop before ending
  int m_loopCount;
  ///The number of times this iterator should loop before ending
  const unsigned int m_loopOrientation;
  /// internal index of location within the workspace
  int m_index;
  ///Internal cache of the workspace size
  int m_wsSize;
  ///Internal cache of the workspace blocksize
  int m_blocksize;

  ///Internal cache of the current datablock index
  int m_dataBlockIndex;
  ///Internal cache of the current datablock index minimum value
  int m_blockMin;
  ///Internal cache of the current datablock index maximum value
  int m_blockMax;

  
  /// @cond
  template<typename T>
  struct internal_iterator_type {};
  
  template<typename T>
  struct internal_iterator_type<T*>
  {
    typedef std::vector<double>::iterator iterator_type;
  };
  
  template<typename T>
  struct internal_iterator_type<const T*>
  {
    typedef std::vector<double>::const_iterator iterator_type;    
  };
  /// @endcond
  
  /// Makes the underlying std::vector<double> iterator be const or not according to whether
  /// the triple_iterator is (or not).
  typedef typename internal_iterator_type<_Iterator*>::iterator_type iterator_type;
  
  ///Internal cache of X iterator for current datablock
  iterator_type it_dataX;
  ///Internal cache of Y iterator for current datablock
  iterator_type it_dataY;
  ///Internal cache of E iterator for current datablock
  iterator_type it_dataE;

  ///Validates the index and updates the current m_CPoint
  void validateIndex();

public:
  /// @cond
  typedef typename std::iterator_traits<_Iterator*>::iterator_category iterator_category;
  typedef typename std::iterator_traits<_Iterator*>::value_type        value_type;
  typedef typename std::iterator_traits<_Iterator*>::difference_type   difference_type;
  typedef typename std::iterator_traits<_Iterator*>::reference         reference;
  typedef typename std::iterator_traits<_Iterator*>::pointer           pointer;
  /// @endcond
  
  triple_iterator();
  triple_iterator(_Container&);
  triple_iterator(_Container&, int loopCount);
  triple_iterator(_Container&, int loopCount, const unsigned int loopOrientation);
  triple_iterator(const triple_iterator&);

  reference operator*() { return m_CPoint; }   ///< Base Accessor
  pointer operator->() { return &m_CPoint; }   ///< Base Pointer accessor
  /// Random accessor
  reference operator[](const difference_type& N)
  {
    m_index=N;
    validateIndex();
    return m_CPoint;
  }
  
  triple_iterator& operator++();
  triple_iterator operator++(int);
  triple_iterator& operator--();
  triple_iterator operator--(int);
  triple_iterator& operator+=(difference_type);
  triple_iterator& operator-=(difference_type);
  triple_iterator operator+(difference_type) const;
  triple_iterator operator-(difference_type) const;
  difference_type operator-(const triple_iterator&) const;

  /*!
   lessthan operator
   \param A :: Iterator to compare
   \return  status
   */
  bool operator<(const triple_iterator& A) const
  {
    if (!m_workspace)
      return 0;
    if (!A.m_workspace)
      return 1;
    return (m_index<A.m_index);
  }

  /*!
   Equality operator
   \param A :: Iterator to compare
   \return equality status
   */
  bool operator==(const triple_iterator& A) const
  {
    if (!m_workspace)
    {
      if (!A.m_workspace)
        return 1;
      return (A.m_wsSize==A.m_index) ? 1 : 0;
    }
    if (!A.m_workspace)
      return (m_wsSize==m_index) ? 1 : 0;

    return (m_index==A.m_index);
  }

  /*!
   InEquality operator
   \param A :: Iterator to compare
   \return equality status
   */
  bool operator!=(const triple_iterator& A) const
  {
    return (!(operator==(A)));
  }

  triple_iterator begin_const() const
  /*!
   Begin iterator (Effective copy+set zero)
   \return beginning iterator
   */
  {
    triple_iterator Out(*this);
    Out.m_index=0;
    Out.validateIndex();
    return Out;
  }

  /*!
   Begin iterator (Set zero)
   \return beginning iterator
   */
  triple_iterator& begin()
  {
    m_index=0;
    validateIndex();
    return *this;
  }

  /*!
   End iterator
   \return Null/end iterator
   */
  triple_iterator& end() const
  {
    static triple_iterator endIter;
    return endIter;
  }

};

/// Describes the orientation of the looping when using a looping iterator
struct LoopOrientation
{
  /// Enum giving the possible orientations
  enum
  { 
    Horizontal=0,    ///< Iterate repeatedly over a horizontal workspace
    Vertical=1       ///< Iterate repeatedly over a vertical workspace
  };
};


} // NAMESPACE API
} // NAMESPACE Mantid

#endif //MANTIDAPI_TRIPLE_ITERATOR_H
