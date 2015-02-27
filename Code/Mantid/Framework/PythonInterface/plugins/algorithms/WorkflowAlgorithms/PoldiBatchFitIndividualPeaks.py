from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
import os

from pyparsing import *


class PoldiSample:
    def __init__(self, initString):

        self._name, self._runList = self.parseInitString(initString)


    def name(self):
        return self._name


    def runList(self):
        return self._runList


    def runNames(self):
        return [x[0] for x in self._runList]


    def parseInitString(self, initString):
        parts = initString.split(':')

        if not len(parts) == 2:
            raise RuntimeError("Unable to parse init string.")

        sampleName = parts[0]
        runDict = self.parseRunDictString(parts[1])

        return (sampleName, runDict)


    def parseRunDictString(self, runDictString):
        runRanges = filter(bool, runDictString.strip().split(';'))

        newRunList = []
        for runRange in runRanges:
            runPoints = self.parseRunRangeString(runRange)
            newRunList += runPoints

        return newRunList


    def parseRunRangeString(self, runRangeString):
        runRangeParameters = [int(x) for x in runRangeString.split(',')]

        if not len(runRangeParameters) == 3:
            raise RuntimeError("Unable to parse run range string.")

        begin, end, step = runRangeParameters

        runList = []
        for i in range(begin, end + 1, step):
            runList.append((i + step - 1, [str(x) for x in range(i, i + step)]))

        return runList


class PoldiProject:
    def __init__(self, projectName, sampleListFile):
        self._name = projectName
        self._sampleList = self.parseSampleListFile(sampleListFile)

    def sampleList(self):
        return self._sampleList

    def name(self):
        return self._name

    def parseSampleListFile(self, sampleListFile):
        fh = open(sampleListFile, 'r')

        sampleList = []
        for line in fh:
            if len(line.strip()) > 0 and not line.strip()[0] == '#':
                sampleList.append(PoldiSample(line))

        return sampleList


class PoldiCompound:
    def __init__(self, name, content, tolerance, elements):
        self._name = name
        self._content = content
        self._tolerance = tolerance

        self.assign(elements)

    def assign(self, elements):
        for c in elements:
            if c[0] == "atoms":
                self._atomString = c[1][0]
            elif c[0] == "lattice":
                pNames = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
                self._latticeDict = dict(zip(pNames, c[1:]))
            elif c[0] == "spacegroup":
                self._spacegroup = c[1]

    def getAtoms(self):
        return self._atomString

    def getLatticeDict(self):
        return self._latticeDict

    def getSpaceGroup(self):
        return self._spacegroup

    def getContent(self):
        return self._content

    def getTolerance(self):
        return self._tolerance

    def getName(self):
        return self._name


class PoldiCrystalFile:
    elementSymbol = Word(alphas, exact=2)
    integerNumber = Word(nums)
    decimalSeparator = Literal('.')
    floatNumber = Combine(
        integerNumber +
        Optional(decimalSeparator + Optional(integerNumber))
    )

    atomLine = Combine(
        elementSymbol + Suppress(White()) +
        delimitedList(floatNumber, delim=White()),
        joinString=' '
    )

    keyValueSeparator = Suppress(Literal(":"))

    atomsGroup = Group(CaselessLiteral("atoms") + keyValueSeparator + nestedExpr(
        opener="{", closer="}",
        content=Combine(
            delimitedList(atomLine, delim='\n'),
            joinString=";")))

    unitCell = Group(CaselessLiteral("lattice") + keyValueSeparator + delimitedList(
        floatNumber, delim=White()))

    spaceGroup = Group(CaselessLiteral("spacegroup") + keyValueSeparator + Word(
        alphanums + "-" + ' '))

    compoundContent = Each([atomsGroup, unitCell, spaceGroup])

    compoundName = Word(alphanums)

    compound = Group(Suppress(CaselessLiteral("compound")) + Suppress(White()) + \
                     compoundName + Suppress(White()) + floatNumber + \
                     Suppress(White()) + floatNumber + \
                     nestedExpr(opener='{', closer='}', content=compoundContent))

    comment = Suppress(Literal('#') + restOfLine)

    compounds = OneOrMore(compound).ignore(comment)

    def __init__(self, filename):
        parsedContent = self.compounds.parseFile(filename)

        self._compounds = []
        for comp in parsedContent:
            self._compounds.append(PoldiCompound(comp[0], comp[1], comp[2], comp[3]))

        self._contributions = [x.getContent() for x in self._compounds]
        self._tolerances = [x.getTolerance() for x in self._compounds]
        self._names = [x.getName() for x in self._compounds]

    def getCompounds(self):
        return self._compounds

    def getContributions(self):
        return self._contributions

    def getTolerances(self):
        return self._tolerances

    def getNames(self):
        return self._names


class PoldiBatchFitIndividualPeaks(PythonAlgorithm):
    progress = None
    outputBase = ""
    crystalFile = ""

    prefixes = {'raw': 'raw_',
                'data': 'data_'}

    suffixes = {'correlation': '_correlation',
                'raw_peaks': '_peaks_raw',
                'fit_peaks_1d': '_peaks_fit_1d',
                'fit_peaks_2d': '_peaks_fit_2d',
                'calculated_2d': '_spectrum_fitted_2d',
                'calculated_1d': '_spectrum_fitted_1d',
                'residuals': '_residuals',
                'total': '_total'}

    def category(self):
        return "Workflow\\Poldi"

    def name(self):
        return "PoldiBatchFitIndividualPeaks"

    def summary(self):
        return "Load a file with instructions for batch processing POLDI fits."

    def PyInit(self):
        self.declareProperty("ProjectName", "", "Project name to use for labeling")
        self.declareProperty(FileProperty(name="SampleListFile", defaultValue="",
                                          action=FileAction.Load))
        self.declareProperty(FileProperty(name="OutputDirectory", defaultValue="",
                                          action=FileAction.Directory))
        self.declareProperty(FileProperty(name="CrystalStructures",
                                          defaultValue="",
                                          action=FileAction.Load))
        self.declareProperty("Year", 2014, "The year in which the files ")
        self.declareProperty("PeakNumber", 10, "Number of peaks to fit.")

    def PyExec(self):
        poldiProject = PoldiProject(self.getProperty("ProjectName").valueAsStr,
                                    self.getProperty("SampleListFile").valueAsStr)

        self.crystalFile = PoldiCrystalFile(self.getProperty("CrystalStructures").valueAsStr)

        self.outputBase = self.getProperty("OutputDirectory").valueAsStr + "/" + poldiProject.name()

        if not os.path.isdir(self.outputBase):
            os.mkdir(self.outputBase)

        for compound in self.crystalFile.getCompounds():
            self.generateCompoundPeaks(compound)

        pointCount = sum([len(x.runList()) for x in poldiProject.sampleList()])
        self.progress = Progress(self, start=0.0, end=1.0, nreports=pointCount)

        for sample in poldiProject.sampleList():
            self.processSample(sample)

    def generateCompoundPeaks(self, compound):
        latticeParameters = compound.getLatticeDict()

        PoldiCreatePeaksFromCell(
            SpaceGroup=compound.getSpaceGroup(),
            a=latticeParameters['a'], b=latticeParameters['b'], c=latticeParameters['c'],
            alpha=latticeParameters['alpha'], beta=latticeParameters['beta'], gamma=latticeParameters['gamma'],
            Atoms=compound.getAtoms(),
            OutputWorkspace=compound.getName()
        )

    def processSample(self, poldiSample):
        sampleList = poldiSample.runList()

        points = []
        for dataPoint in sampleList:
            points.append(dataPoint[0])
            try:
                self.getLogger().error("Processing %s" % (dataPoint[0]))
                self.processDataPoint(dataPoint)
            except RuntimeError as error:
                self.getLogger().error("An error occured when processing %s: '%s'" % (dataPoint[0], error.message))
            self.progress.report()

        sampleDir = self.outputBase + "/" + poldiSample.name()

        if not os.path.isdir(sampleDir):
            os.mkdir(sampleDir)

        for phase in self.crystalFile.getNames():
            peaks = []
            for p in points:
                peaks.append(self.extractPeakInformation("Indexed_%s_data_%s_peaks_fit_2d" % (phase, str(p))))

            hkls = peaks[0].keys()

            for h in hkls:
                hklFile = open(sampleDir + "/%s_%s_%s.out" % (poldiSample.name(), phase, h), 'w')
                hklFile.write("# Run Q delta_Q d delta_d FWHMrel delta_FWHMrel Intensity delta_Intensity\n")

                for i, d in enumerate(peaks):
                    hklFile.write(str(points[i]) + ' ')
                    if h in d:
                        current = [str(x[0]) + ' ' + str(x[1]) for x in d[h]]
                        hklFile.write(' '.join(current))
                    hklFile.write('\n')

                hklFile.close()


    def extractPeakInformation(self, wsName):
        peakDict = {}
        try:
            tableWorkspace = mtd[wsName]

            keys = ['Q', 'd', 'FWHM (rel.)', 'Intensity']

            for i in range(tableWorkspace.rowCount()):
                current = tableWorkspace.row(i)

                da = []
                for k in keys:
                    data = current[k].split()
                    da.append((float(data[0]), float(data[-1])))

                peakDict[current['HKL']] = da
        except:
            pass

        return peakDict


    def processDataPoint(self, dataPointTuple):
        self.loadDataPoint(dataPointTuple)
        self.runPoldiAnalysis(dataPointTuple[0])

    def loadDataPoint(self, dataPointTuple):
        dataWorkspaceName = self.prefixes['data'] + str(dataPointTuple[0])

        rawWorkspaceNames = [
            self.prefixes['raw'] + str(x) for x in dataPointTuple[1]]

        self.getLogger().error("Loading " + str(dataPointTuple[1]))

        year = int(self.getProperty("Year").valueAsStr)

        for numor, rawWsName in zip(dataPointTuple[1], rawWorkspaceNames):
            LoadSINQ(Instrument="POLDI", Year=year,
                     Numor=numor, OutputWorkspace=rawWsName)
            LoadInstrument(rawWsName, InstrumentName="POLDI")
            PoldiTruncateData(InputWorkspace=rawWsName, OutputWorkspace=rawWsName)


        if len(dataPointTuple[1]) > 1:
            PoldiMerge(rawWorkspaceNames, OutputWorkspace=dataWorkspaceName)

            for ws in rawWorkspaceNames:
                DeleteWorkspace(ws)

        else:
            RenameWorkspace(rawWorkspaceNames[0], OutputWorkspace=dataWorkspaceName)

        return dataWorkspaceName

    def runPoldiAnalysis(self, dataPointName):
        dataWorkspaceName = self.prefixes['data'] + str(dataPointName)

        correlationWsName = dataWorkspaceName + self.suffixes['correlation']
        PoldiAutoCorrelation(dataWorkspaceName, wlenmin=1.1, wlenmax=5.0,
                             OutputWorkspace=correlationWsName)

        peakSearchWsName = dataWorkspaceName + self.suffixes['raw_peaks']
        PoldiPeakSearch(correlationWsName, OutputWorkspace=peakSearchWsName)

        DeleteTableRows(peakSearchWsName,
                        self.getProperty("PeakNumber").valueAsStr + '-30')

        peakFit1DWsName = dataWorkspaceName + self.suffixes['fit_peaks_1d']
        PoldiFitPeaks1D(correlationWsName, PoldiPeakTable=peakSearchWsName,
                        OutputWorkspace=peakFit1DWsName)

        peakFit2DWsName = dataWorkspaceName + self.suffixes['fit_peaks_2d']
        spectrum2DCalc = dataWorkspaceName + self.suffixes['calculated_2d']
        spectrum1DCalc = dataWorkspaceName + self.suffixes['calculated_1d']
        PoldiFitPeaks2D(dataWorkspaceName, PoldiPeakWorkspace=peakFit1DWsName,
                        MaximumIterations=30, OutputWorkspace=spectrum2DCalc,
                        RefinedPoldiPeakWorkspace=peakFit2DWsName,
                        Calculated1DSpectrum=spectrum1DCalc)

        residualsWsName = dataWorkspaceName + self.suffixes['residuals']
        PoldiAnalyseResiduals(dataWorkspaceName, spectrum2DCalc,
                              OutputWorkspace=residualsWsName,
                              MaxIterations=5)

        totalWsName = dataWorkspaceName + self.suffixes['total']
        Plus(LHSWorkspace=spectrum1DCalc, RHSWorkspace=residualsWsName, OutputWorkspace=totalWsName)

        indexedPeaksWs = dataWorkspaceName + "_indexed"
        PoldiIndexKnownCompounds(peakFit2DWsName,
                                 CompoundWorkspaces=','.join(self.crystalFile.getNames()),
                                 Tolerances=','.join(self.crystalFile.getTolerances()),
                                 ScatteringContributions=','.join(self.crystalFile.getContributions()),
                                 OutputWorkspace=indexedPeaksWs)

        return indexedPeaksWs


AlgorithmFactory.subscribe(PoldiBatchFitIndividualPeaks)