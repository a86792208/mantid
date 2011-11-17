#include "MantidQwtMatrixWorkspaceData.h"

/// Constructor
MantidQwtMatrixWorkspaceData::MantidQwtMatrixWorkspaceData(Mantid::API::MatrixWorkspace_const_sptr workspace,int specIndex, const bool logScale, bool distr)
: QObject(),
m_workspace(workspace),
m_spec(specIndex),
m_X(workspace->readX(specIndex)),
m_Y(workspace->readY(specIndex)),
m_E(workspace->readE(specIndex)),
m_isHistogram(workspace->isHistogramData()),
m_binCentres(false),
m_logScale(logScale),
m_minPositive(0),
m_isDistribution(distr)
{}

/// Copy constructor
MantidQwtMatrixWorkspaceData::MantidQwtMatrixWorkspaceData(const MantidQwtMatrixWorkspaceData& data)
: QObject(),
m_workspace(data.m_workspace),
m_spec(data.m_spec),
m_X(data.m_workspace->readX(data.m_spec)),
m_Y(data.m_workspace->readY(data.m_spec)),
m_E(data.m_workspace->readE(data.m_spec)),
m_isHistogram(m_workspace->isHistogramData()),
m_binCentres(data.m_binCentres),
m_logScale(data.m_logScale),
m_minPositive(0),
m_isDistribution(data.m_isDistribution)
{}

/** Size of the data set
 */
size_t MantidQwtMatrixWorkspaceData::size() const
{
  if (m_binCentres || !m_isHistogram)
  {
    return m_Y.size();
  }

  return m_X.size();
}

/**
Return the x value of data point i
@param i :: Index
@return x X value of data point i
*/
double MantidQwtMatrixWorkspaceData::x(size_t i) const
{
  return m_binCentres ? (m_X[i] + m_X[i+1])/2 : m_X[i];
}

/**
Return the y value of data point i
@param i :: Index
@return y Y value of data point i
*/
double MantidQwtMatrixWorkspaceData::y(size_t i) const
{
  double tmp = i < m_Y.size() ? m_Y[i] : m_Y[m_Y.size()-1];
  if (m_isDistribution)
  {
    tmp /= (m_X[i+1] - m_X[i]);
  }
  if (m_logScale && tmp <= 0.)
  {
    tmp = m_minPositive;
  }

  return tmp;
}

double MantidQwtMatrixWorkspaceData::ex(size_t i) const
{
  return m_isHistogram ? (m_X[i] + m_X[i+1])/2 : m_X[i];
}

double MantidQwtMatrixWorkspaceData::e(size_t i) const
{
  return m_E[i];
}

size_t MantidQwtMatrixWorkspaceData::esize() const
{
  return m_E.size();
}

bool MantidQwtMatrixWorkspaceData::sameWorkspace(boost::shared_ptr<const Mantid::API::MatrixWorkspace> workspace)const
{
  return workspace.get() == m_workspace.get();
}

void MantidQwtMatrixWorkspaceData::setLogScale(bool on)
{
  m_logScale = on;
}

void MantidQwtMatrixWorkspaceData::saveLowestPositiveValue(const double v)
{
  if (v > 0) m_minPositive = v;
}

void MantidQwtMatrixWorkspaceData::applyOffsets(const double xOffset, const double yOffset)
{
  std::transform(m_workspace->readX(m_spec).begin(),m_workspace->readX(m_spec).end(),m_X.begin(),std::bind2nd(std::plus<double>(),xOffset));
  std::transform(m_workspace->readY(m_spec).begin(),m_workspace->readY(m_spec).end(),m_Y.begin(),std::bind2nd(std::plus<double>(),yOffset));
}

bool MantidQwtMatrixWorkspaceData::setAsDistribution(bool on)
{
  m_isDistribution = on && m_isHistogram;
  return m_isDistribution;
}
