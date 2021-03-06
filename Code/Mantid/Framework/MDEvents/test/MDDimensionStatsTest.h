#ifndef MANTID_MDEVENTS_MDDIMENSIONSTATSTEST_H_
#define MANTID_MDEVENTS_MDDIMENSIONSTATSTEST_H_

#include <cxxtest/TestSuite.h>
#include "MantidKernel/Timer.h"
#include "MantidKernel/System.h"
#include <iostream>
#include <iomanip>

#include "MantidMDEvents/MDDimensionStats.h"

using namespace Mantid::MDEvents;

class MDDimensionStatsTest : public CxxTest::TestSuite
{
public:

  void test_constructor()
  {
    MDDimensionStats s;
    TS_ASSERT_EQUALS(s.total, 0.0);
    TS_ASSERT_EQUALS(s.totalApproxVariance, 0.0);
    TS_ASSERT_EQUALS(s.numPoints, 0);
  }


  void test_addPoint()
  {
    MDDimensionStats s;
    s.addPoint(10.0);
    TS_ASSERT_DELTA( s.getMean(), 10.0, 1e-4);
    TS_ASSERT_DELTA( s.getApproxVariance(), 0.0, 1e-4);
    s.addPoint(0.0);
    TS_ASSERT_DELTA( s.getMean(), 5.0, 1e-4);
    TS_ASSERT_DELTA( s.getApproxVariance(), 25.0/2.0, 1e-4);
    s.addPoint(5.0);
    TS_ASSERT_DELTA( s.getMean(), 5.0, 1e-4);
    TS_ASSERT_DELTA( s.getApproxVariance(), 25.0/3.0, 1e-4);
  }


};


#endif /* MANTID_MDEVENTS_MDDIMENSIONSTATSTEST_H_ */

