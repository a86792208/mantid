#include "MantidMDEvents/ReflectometryMDTransform.h"

using namespace Mantid::API;

namespace Mantid
{
  namespace MDEvents
  {

    ReflectometryMDTransform::ReflectometryMDTransform() : m_nbinsx(10), m_nbinsz(10)
    {
    }

    ReflectometryMDTransform::~ReflectometryMDTransform()
    {
    }

    boost::shared_ptr<MDEventWorkspace2Lean> ReflectometryMDTransform::createMDWorkspace(
        Mantid::Geometry::IMDDimension_sptr a, Mantid::Geometry::IMDDimension_sptr b, BoxController_sptr boxController) const
    {
      auto ws = boost::make_shared<MDEventWorkspace2Lean>();

      ws->addDimension(a);
      ws->addDimension(b);

      BoxController_sptr wsbc = ws->getBoxController();// Get the box controller
      wsbc->setSplitInto(boxController->getSplitInto(0));
      wsbc->setMaxDepth(boxController->getMaxDepth());
      wsbc->setSplitThreshold(boxController->getSplitThreshold());

      // Initialize the workspace.
      ws->initialize();

      // Start with a MDGridBox.
      ws->splitBox();
      return ws;
    }

  }
}
