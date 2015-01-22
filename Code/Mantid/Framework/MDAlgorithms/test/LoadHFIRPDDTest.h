#ifndef MANTID_DATAHANDLING_LOADHFIRPDDTEST_H_
#define MANTID_DATAHANDLING_LOADHFIRPDDTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidMDAlgorithms/LoadHFIRPDD.h"
// #include "MantidDataHandling/LoadNexusProcessed.h"
#include "MantidDataHandling/LoadSpiceAscii.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidDataHandling/LoadInstrument.h"
#include "MantidAPI/IMDIterator.h"
#include "MantidMDAlgorithms/SaveMD.h"
#include "MantidGeometry/IComponent.h"

// #include "Mantid"

using Mantid::DataHandling::LoadHFIRPDD;
// using Mantid::DataHandling::LoadNexusProcessed;
using Mantid::DataHandling::LoadInstrument;
using Mantid::MDAlgorithms::SaveMD;
using Mantid::DataHandling::LoadSpiceAscii;

using namespace Mantid::API;
using namespace Mantid::DataObjects;

class LoadHFIRPDDTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static LoadHFIRPDDTest *createSuite() { return new LoadHFIRPDDTest(); }
  static void destroySuite(LoadHFIRPDDTest *suite) { delete suite; }

  void test_Init() {
    LoadHFIRPDD loader;
    loader.initialize();
    TS_ASSERT(loader.isInitialized());
  }

  /** Test loading HB2A's IDF file
   * @brief test_HB2BIDF
   */
  void test_HB2AIDF() {
    MatrixWorkspace_sptr dataws =
        WorkspaceFactory::Instance().create("Workspace2D", 44, 2, 1);
    AnalysisDataService::Instance().addOrReplace("EmptyWS", dataws);

    LoadInstrument loader;
    loader.initialize();

    loader.setProperty("InstrumentName", "HB2A");
    loader.setProperty("Workspace", dataws);

    loader.execute();
    TS_ASSERT(loader.isExecuted());

    MatrixWorkspace_sptr outws = boost::dynamic_pointer_cast<MatrixWorkspace>(
        AnalysisDataService::Instance().retrieve("EmptyWS"));
    TS_ASSERT(outws);

    Mantid::Geometry::Instrument_const_sptr hb2a = outws->getInstrument();
    TS_ASSERT(hb2a);
    if (!hb2a)
      return;

    TS_ASSERT_EQUALS(hb2a->getName(), "HB2A");

    Mantid::Geometry::IComponent_const_sptr source = hb2a->getSource();
    TS_ASSERT(source);
    Mantid::Geometry::IComponent_const_sptr sample = hb2a->getSample();
    TS_ASSERT(sample);
    const Mantid::Kernel::V3D &samplepos = sample->getPos();
    TS_ASSERT_DELTA(samplepos.X(), 0.0, 0.00001);
    std::vector<Mantid::detid_t> detids = hb2a->getDetectorIDs();
    TS_ASSERT_EQUALS(detids.size(), 44);
  }

  void test_LoadHB2AData() {
    // This is the solution of stage 1. Stage 2 will be different
#if 0
    LoadNexusProcessed nxsloader;
    nxsloader.initialize();
    nxsloader.setProperty("Filename", "HB2A_exp0231_tabledata.nxs");
    nxsloader.setProperty("OutputWorkspace", "DataTable");
    nxsloader.execute();
    ITableWorkspace_sptr datatablews =
        boost::dynamic_pointer_cast<ITableWorkspace>(
            AnalysisDataService::Instance().retrieve("DataTable"));
    TS_ASSERT(datatablews);

    nxsloader.setProperty("Filename", "HB2A_exp0231_log.nxs");
    nxsloader.setProperty("OutputWorkspace", "LogParentWS");
    nxsloader.execute();
    MatrixWorkspace_sptr parentlogws =
        boost::dynamic_pointer_cast<MatrixWorkspace>(
            AnalysisDataService::Instance().retrieve("LogParentWS"));
    TS_ASSERT(parentlogws);
#endif

    LoadSpiceAscii spcloader;
    spcloader.initialize();

    TS_ASSERT_THROWS_NOTHING(
        spcloader.setProperty("Filename", "HB2A_exp0231_scan0001.dat"));
    TS_ASSERT_THROWS_NOTHING(
        spcloader.setProperty("OutputWorkspace", "DataTable"));
    TS_ASSERT_THROWS_NOTHING(
        spcloader.setProperty("RunInfoWorkspace", "LogParentWS"));
    TS_ASSERT_THROWS_NOTHING(spcloader.setPropertyValue(
        "DateAndTimeLog", "date,MM/DD/YYYY,time,HH:MM:SS AM"));
    TS_ASSERT_THROWS_NOTHING(
        spcloader.setProperty("IgnoreUnlistedLogs", false));

    spcloader.execute();

    ITableWorkspace_sptr datatablews =
        boost::dynamic_pointer_cast<ITableWorkspace>(
            AnalysisDataService::Instance().retrieve("DataTable"));
    TS_ASSERT(datatablews);

    MatrixWorkspace_sptr parentlogws =
        boost::dynamic_pointer_cast<MatrixWorkspace>(
            AnalysisDataService::Instance().retrieve("LogParentWS"));
    TS_ASSERT(parentlogws);

    LoadHFIRPDD loader;
    loader.initialize();

    // loader.setPropertyValue("Filename", "HB2A_exp0231_scan0001.dat");
    loader.setProperty("InputWorkspace", datatablews);
    loader.setProperty("ParentWorkspace", parentlogws);
    loader.setProperty("Instrument", "HB2A");
    // loader.setProperty("RunStart", "2012-08-13T11:57:33");
    loader.setPropertyValue("OutputWorkspace", "HB2A_MD");
    loader.setPropertyValue("OutputMonitorWorkspace", "MonitorMDW");
    loader.setPropertyValue("OutputReducedWorkspace", "ReducedWS");

    loader.execute();
    TS_ASSERT(loader.isExecuted());

    // Get IMDEventWorkspace
    IMDEventWorkspace_sptr mdws =
        boost::dynamic_pointer_cast<IMDEventWorkspace>(
            AnalysisDataService::Instance().retrieve("HB2A_MD"));
    TS_ASSERT(mdws);

    // Check the IMDEvent workspace generated
    size_t numevents = mdws->getNEvents();
    TS_ASSERT_EQUALS(numevents, 44 * 61);

    IMDIterator *mditer = mdws->createIterator();
    TS_ASSERT_EQUALS(mditer->getNumEvents(), 44 * 61);

    Mantid::coord_t lastx = mditer->getInnerPosition(44 * 61 - 1, 0);
    TS_ASSERT_DELTA(lastx, 1.57956, 0.0001);

    ExperimentInfo_const_sptr expinfo = mdws->getExperimentInfo(0);
    TS_ASSERT(expinfo);

    SaveMD saver;
    saver.initialize();
    saver.setProperty("InputWorkspace", mdws);
    saver.setProperty("Filename", "/tmp/hb2amd.nxs");
    saver.execute();
  }
};

#endif /* MANTID_DATAHANDLING_LOADHFIRPDDTEST_H_ */
