#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/WorkspaceGroup.h"

#include "MantidSINQ/PoldiFitPeaks1D.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidAPI/TableRow.h"

#include "MantidSINQ/PoldiUtilities/UncertainValue.h"
#include "MantidSINQ/PoldiUtilities/UncertainValueIO.h"

#include "MantidAPI/CompositeFunction.h"


namespace Mantid
{
namespace Poldi
{

using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::CurveFitting;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(PoldiFitPeaks1D)


PoldiFitPeaks1D::PoldiFitPeaks1D() :
    m_peaks(),
    m_profileTemplate(),
    m_backgroundTemplate(),
    m_profileTies(),
    m_fitCharacteristics(),
    m_peakResultOutput(),
    m_fwhmMultiples(1.0)
{

}

PoldiFitPeaks1D::~PoldiFitPeaks1D()
{
}


/// Algorithm's name for identification. @see Algorithm::name
const std::string PoldiFitPeaks1D::name() const { return "PoldiFitPeaks1D";}

/// Algorithm's version for identification. @see Algorithm::version
int PoldiFitPeaks1D::version() const { return 1;}

/// Algorithm's category for identification. @see Algorithm::category
const std::string PoldiFitPeaks1D::category() const { return "SINQ\\Poldi"; }

void PoldiFitPeaks1D::init()
{
    declareProperty(new WorkspaceProperty<Workspace2D>("InputWorkspace","",Direction::Input), "An input workspace containing a POLDI auto-correlation spectrum.");
    boost::shared_ptr<BoundedValidator<double> > minFwhmPerDirection = boost::make_shared<BoundedValidator<double> >();
    minFwhmPerDirection->setLower(2.0);
    declareProperty("FwhmMultiples", 2.0, minFwhmPerDirection, "Each peak will be fitted using x times FWHM data in each direction.", Direction::Input);

    std::vector<std::string> peakFunctions = FunctionFactory::Instance().getFunctionNames<IPeakFunction>();
    boost::shared_ptr<ListValidator<std::string> > peakFunctionNames(new ListValidator<std::string>(peakFunctions));
    declareProperty("PeakFunction", "Gaussian", peakFunctionNames, "Peak function that will be fitted to all peaks.", Direction::Input);

    declareProperty(new WorkspaceProperty<TableWorkspace>("PoldiPeakTable","",Direction::Input), "A table workspace containing POLDI peak data.");

    declareProperty(new WorkspaceProperty<TableWorkspace>("OutputWorkspace","RefinedPeakTable",Direction::Output), "Output workspace with refined peak data.");
    declareProperty(new WorkspaceProperty<TableWorkspace>("ResultTableWorkspace","ResultTable",Direction::Output), "Fit results.");
    declareProperty(new WorkspaceProperty<TableWorkspace>("FitCharacteristicsWorkspace","FitCharacteristics",Direction::Output), "Fit characteristics for each peak.");
    declareProperty(new WorkspaceProperty<Workspace>("FitPlotsWorkspace","FitPlots",Direction::Output), "Plots of all peak fits.");

    m_backgroundTemplate = FunctionFactory::Instance().createInitialized("name=UserFunction, Formula=A0 + A1*(x - x0)^2");
    m_profileTies = "f1.x0 = f0.PeakCentre";
}

void PoldiFitPeaks1D::setPeakFunction(std::string peakFunction)
{
    m_profileTemplate = peakFunction;
}

PoldiPeakCollection_sptr PoldiFitPeaks1D::getInitializedPeakCollection(TableWorkspace_sptr peakTable)
{
    return PoldiPeakCollection_sptr(new PoldiPeakCollection(peakTable));
}

IFunction_sptr PoldiFitPeaks1D::getPeakProfile(PoldiPeak_sptr poldiPeak) {
    IPeakFunction_sptr clonedProfile = boost::dynamic_pointer_cast<IPeakFunction>(FunctionFactory::Instance().createFunction(m_profileTemplate));
    clonedProfile->setCentre(poldiPeak->q());
    clonedProfile->setFwhm(poldiPeak->fwhm(PoldiPeak::AbsoluteQ));
    clonedProfile->setHeight(poldiPeak->intensity());

    IFunction_sptr clonedBackground = m_backgroundTemplate->clone();

    CompositeFunction_sptr totalProfile(new CompositeFunction);
    totalProfile->initialize();
    totalProfile->addFunction(clonedProfile);
    totalProfile->addFunction(clonedBackground);

    if(!m_profileTies.empty()) {
        totalProfile->addTies(m_profileTies);
    }

    return totalProfile;
}

void PoldiFitPeaks1D::setValuesFromProfileFunction(PoldiPeak_sptr poldiPeak, IFunction_sptr fittedFunction)
{
    CompositeFunction_sptr totalFunction = boost::dynamic_pointer_cast<CompositeFunction>(fittedFunction);

    if(totalFunction) {
        IPeakFunction_sptr peakFunction = boost::dynamic_pointer_cast<IPeakFunction>(totalFunction->getFunction(0));

        if(peakFunction) {
            poldiPeak->setIntensity(UncertainValue(peakFunction->height(), peakFunction->getError(0)));
            poldiPeak->setQ(UncertainValue(peakFunction->centre(), peakFunction->getError(1)));
            poldiPeak->setFwhm(UncertainValue(peakFunction->fwhm(), getFwhmWidthRelation(peakFunction) * peakFunction->getError(2)));
        }
    }
}

double PoldiFitPeaks1D::getFwhmWidthRelation(IPeakFunction_sptr peakFunction)
{
    return peakFunction->fwhm() / peakFunction->getParameter(2);
}

void PoldiFitPeaks1D::exec()
{
    setPeakFunction(getProperty("PeakFunction"));

    // Number of points around the peak center to use for the fit
    m_fwhmMultiples = getProperty("FwhmMultiples");

    // try to construct PoldiPeakCollection from provided TableWorkspace
    TableWorkspace_sptr poldiPeakTable = getProperty("PoldiPeakTable");    
    m_peaks = getInitializedPeakCollection(poldiPeakTable);

    g_log.information() << "Peaks to fit: " << m_peaks->peakCount() << std::endl;

    Workspace2D_sptr dataWorkspace = getProperty("InputWorkspace");

    m_fitCharacteristics = boost::dynamic_pointer_cast<TableWorkspace>(WorkspaceFactory::Instance().createTable());
    WorkspaceGroup_sptr fitPlotGroup(new WorkspaceGroup);

    for(size_t i = 0; i < m_peaks->peakCount(); ++i) {
        PoldiPeak_sptr currentPeak = m_peaks->peak(i);
        IFunction_sptr currentProfile = getPeakProfile(currentPeak);

        IAlgorithm_sptr fit = getFitAlgorithm(dataWorkspace, currentPeak, currentProfile);

        bool fitSuccess = fit->execute();

        if(fitSuccess) {
            setValuesFromProfileFunction(currentPeak, fit->getProperty("Function"));
            addPeakFitCharacteristics(fit->getProperty("OutputParameters"));

            MatrixWorkspace_sptr fpg = fit->getProperty("OutputWorkspace");
            fitPlotGroup->addWorkspace(fpg);
        }
    }

    m_peakResultOutput = generateResultTable(m_peaks);

    setProperty("OutputWorkspace", m_peaks->asTableWorkspace());
    setProperty("FitCharacteristicsWorkspace", m_fitCharacteristics);
    setProperty("ResultTableWorkspace", m_peakResultOutput);
    setProperty("FitPlotsWorkspace", fitPlotGroup);
}

IAlgorithm_sptr PoldiFitPeaks1D::getFitAlgorithm(Workspace2D_sptr dataWorkspace, PoldiPeak_sptr peak, IFunction_sptr profile)
{
    double width = peak->fwhm();
    double extent = std::min(0.05, std::max(0.002, width)) * m_fwhmMultiples;

    std::pair<double, double> xBorders(peak->q() - extent, peak->q() + extent);

    IAlgorithm_sptr fitAlgorithm = createChildAlgorithm("Fit", -1, -1, false);
    fitAlgorithm->setProperty("CreateOutput", true);
    fitAlgorithm->setProperty("Output", "FitPeaks1D");
    fitAlgorithm->setProperty("CalcErrors", true);
    fitAlgorithm->setProperty("Function", profile);
    fitAlgorithm->setProperty("InputWorkspace", dataWorkspace);
    fitAlgorithm->setProperty("WorkspaceIndex", 0);
    fitAlgorithm->setProperty("StartX", xBorders.first);
    fitAlgorithm->setProperty("EndX", xBorders.second);

    return fitAlgorithm;
}

void PoldiFitPeaks1D::addPeakFitCharacteristics(ITableWorkspace_sptr fitResult)
{
    if(m_fitCharacteristics->columnCount() == 0) {
        initializeFitResultWorkspace(fitResult);
    }

    TableRow newRow = m_fitCharacteristics->appendRow();

    for(size_t i = 0; i < fitResult->rowCount(); ++i) {
        TableRow currentRow = fitResult->getRow(i);

        newRow << UncertainValueIO::toString(UncertainValue(currentRow.Double(1), currentRow.Double(2)));
    }
}

void PoldiFitPeaks1D::initializeFitResultWorkspace(ITableWorkspace_sptr fitResult)
{
    for(size_t i = 0; i < fitResult->rowCount(); ++i) {
        TableRow currentRow = fitResult->getRow(i);
        m_fitCharacteristics->addColumn("str", currentRow.cell<std::string>(0));
    }
}

void PoldiFitPeaks1D::initializePeakResultWorkspace(TableWorkspace_sptr peakResultWorkspace)
{
    peakResultWorkspace->addColumn("str", "Q");
    peakResultWorkspace->addColumn("str", "d");
    peakResultWorkspace->addColumn("double", "deltaD/d *10^3");
    peakResultWorkspace->addColumn("str", "FWHM rel. *10^3");
    peakResultWorkspace->addColumn("str", "Intensity");
}

void PoldiFitPeaks1D::storePeakResult(TableRow tableRow, PoldiPeak_sptr peak)
{
    UncertainValue q = peak->q();
    UncertainValue d = peak->d();

    tableRow << UncertainValueIO::toString(q)
             << UncertainValueIO::toString(d)
             << d.error() / d.value() * 1e3
             << UncertainValueIO::toString(peak->fwhm(PoldiPeak::Relative) * 1e3)
             << UncertainValueIO::toString(peak->intensity());
}

TableWorkspace_sptr PoldiFitPeaks1D::generateResultTable(PoldiPeakCollection_sptr peaks)
{
    TableWorkspace_sptr outputTable = boost::dynamic_pointer_cast<TableWorkspace>(WorkspaceFactory::Instance().createTable());
    initializePeakResultWorkspace(outputTable);

    for(size_t i = 0; i < peaks->peakCount(); ++i) {
        storePeakResult(outputTable->appendRow(), peaks->peak(i));
    }

    return outputTable;
}


} // namespace Poldi
} // namespace Mantid
