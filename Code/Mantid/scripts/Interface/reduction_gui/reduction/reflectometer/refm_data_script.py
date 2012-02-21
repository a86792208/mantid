"""
    Reduction script for REFM
"""
import xml.dom.minidom
import os
import time
from reduction_gui.reduction.scripter import BaseScriptElement

class DataSets(BaseScriptElement):

    DataPeakPixels = [215, 225]
    DataBackgroundFlag = False
    DataBackgroundRoi = [115, 137,123, 137]
    DataTofRange = [10700., 24500.]
    crop_TOF_range = True
    
    data_x_range_flag = True
    data_x_range = [115,210]
    
    norm_x_range_flag = True
    norm_x_range = [90, 160]

    NormFlag = True
    NormPeakPixels = [120, 130]
    NormBackgroundFlag = False
    NormBackgroundRoi = [115, 137]

    # Data files
    data_files = [0]
    norm_file = 0
    
    # Q range
    q_min = 0.0025
    q_step = -0.01
    
    # scattering angle
    theta = 0.0
    center_pixel = 220
    use_center_pixel = False

    def __init__(self):
        super(DataSets, self).__init__()
        self.reset()

    def to_script(self, for_automated_reduction=False):
        """
            Generate reduction script
            @param execute: if true, the script will be executed
        """
        if for_automated_reduction:
            return self._automated_reduction()
        
        script =  "RefMReduction(RunNumbers=%s,\n" % ','.join([str(i) for i in self.data_files])
        script += "              NormalizationRunNumber=%d,\n" % self.norm_file
        script += "              SignalPeakPixelRange=%s,\n" % str(self.DataPeakPixels)
        script += "              SubtractSignalBackground=%s,\n" % str(self.DataBackgroundFlag)
        script += "              SignalBackgroundPixelRange=%s,\n" % str(self.DataBackgroundRoi[:2])
        script += "              PerformNormalization=%s,\n" % str(self.NormFlag)
        script += "              NormPeakPixelRange=%s,\n" % str(self.NormPeakPixels)
        script += "              NormBackgroundPixelRange=%s,\n" % str(self.NormBackgroundRoi)
        script += "              SubtractNormBackground=%s,\n" % str(self.NormBackgroundFlag)
                        
        script += "              CropLowResDataAxis=%s,\n" % str(self.data_x_range_flag)
        if self.data_x_range_flag:
            script += "              LowResDataAxisPixelRange=%s,\n" % str(self.data_x_range)
            
        script += "              CropLowResNormAxis=%s,\n" % str(self.norm_x_range_flag)
        if self.norm_x_range_flag:
            script += "              LowResNormAxisPixelRange=%s,\n" % str(self.norm_x_range)
            
        if self.crop_TOF_range:
            script += "              TOFRange=%s,\n" % str(self.DataTofRange)
        script += "              QMin=%s,\n" % str(self.q_min)
        script += "              QStep=%s,\n" % str(self.q_step)
        
        # Angle offset
        if self.use_center_pixel:
            script += "              CenterPixel=%s,\n" % str(self.center_pixel)
        else:
            script += "              Theta=%s,\n" % str(self.theta)
            
        # The output should be slightly different if we are generating
        # a script for the automated reduction
        script += "              OutputWorkspace='reflectivity_Off_Off_%s')" % str(self.data_files[0])
        script += "\n"

        return script

    def _automated_reduction(self):
        script =  "# REF_M automated reduction\n"
        
        script += "estimates = RefEstimates(RunNumber=runNumber)\n"
        script += "peak_min = estimates.getProperty('PeakMin').value\n"
        script += "peak_max = estimates.getProperty('PeakMax').value\n"
        script += "ref_pixel = (peak_max+peak_min)/2.0\n\n"
        
        script += "estimates = RefEstimates(RunNumber=%s)\n" % str(self.norm_file)
        script += "peak_min_norm = estimates.getProperty('PeakMin').value\n"
        script += "peak_max_norm = estimates.getProperty('PeakMax').value\n\n"
        
        script += "RefMReduction(RunNumbers=%s,\n" % ','.join([str(i) for i in self.data_files])
        script += "              NormalizationRunNumber=%d,\n" % self.norm_file
        script += "              SignalPeakPixelRange=[peak_min, peak_max],\n"
        script += "              SubtractSignalBackground=False,\n"
        script += "              PerformNormalization=%s,\n" % str(self.NormFlag)
        script += "              NormPeakPixelRange=[peak_min_norm, peak_max_norm],\n"
        script += "              SubtractNormBackground=False,\n"
                        
        script += "              CropLowResDataAxis=%s,\n" % str(self.data_x_range_flag)
        if self.data_x_range_flag:
            script += "              LowResDataAxisPixelRange=%s,\n" % str(self.data_x_range)
            
        script += "              CropLowResNormAxis=%s,\n" % str(self.norm_x_range_flag)
        if self.norm_x_range_flag:
            script += "              LowResNormAxisPixelRange=%s,\n" % str(self.norm_x_range)
            
        script += "              QMin=%s,\n" % str(self.q_min)
        script += "              QStep=%s,\n" % str(self.q_step)
        script += "              Theta=ref_pixel,\n"
            
        # The output should be slightly different if we are generating
        # a script for the automated reduction
        script += "              OutputWorkspace='reflectivity_Off_Off_'+%s)\n" % str(self.data_files[0])

        return script

    def update(self):
        """
            Update transmission from reduction output
        """
        pass

    def to_xml(self):
        """
            Create XML from the current data.
        """
        xml  = "<RefMData>\n"
        xml += "<from_peak_pixels>%s</from_peak_pixels>\n" % str(self.DataPeakPixels[0])
        xml += "<to_peak_pixels>%s</to_peak_pixels>\n" % str(self.DataPeakPixels[1])
        xml += "<background_flag>%s</background_flag>\n" % str(self.DataBackgroundFlag)
        xml += "<back_roi1_from>%s</back_roi1_from>\n" % str(self.DataBackgroundRoi[0])
        xml += "<back_roi1_to>%s</back_roi1_to>\n" % str(self.DataBackgroundRoi[1])
        xml += "<back_roi2_from>%s</back_roi2_from>\n" % str(self.DataBackgroundRoi[2])
        xml += "<back_roi2_to>%s</back_roi2_to>\n" % str(self.DataBackgroundRoi[3])
        xml += "<crop_tof>%s</crop_tof>\n" % str(self.crop_TOF_range)
        xml += "<from_tof_range>%s</from_tof_range>\n" % str(self.DataTofRange[0])
        xml += "<to_tof_range>%s</to_tof_range>\n" % str(self.DataTofRange[1])
        xml += "<data_sets>%s</data_sets>\n" % ','.join([str(i) for i in self.data_files])
        xml += "<x_min_pixel>%s</x_min_pixel>\n" % str(self.data_x_range[0])
        xml += "<x_max_pixel>%s</x_max_pixel>\n" % str(self.data_x_range[1])
        xml += "<x_range_flag>%s</x_range_flag>\n" % str(self.data_x_range_flag)

        xml += "<norm_flag>%s</norm_flag>\n" % str(self.NormFlag)
        xml += "<norm_x_range_flag>%s</norm_x_range_flag>\n" % str(self.norm_x_range_flag)
        xml += "<norm_x_max>%s</norm_x_max>\n" % str(self.norm_x_range[1])
        xml += "<norm_x_min>%s</norm_x_min>\n" % str(self.norm_x_range[0])
        
        xml += "<norm_from_peak_pixels>%s</norm_from_peak_pixels>\n" % str(self.NormPeakPixels[0])
        xml += "<norm_to_peak_pixels>%s</norm_to_peak_pixels>\n" % str(self.NormPeakPixels[1])
        xml += "<norm_background_flag>%s</norm_background_flag>\n" % str(self.NormBackgroundFlag)
        xml += "<norm_from_back_pixels>%s</norm_from_back_pixels>\n" % str(self.NormBackgroundRoi[0])
        xml += "<norm_to_back_pixels>%s</norm_to_back_pixels>\n" % str(self.NormBackgroundRoi[1])
        xml += "<norm_dataset>%s</norm_dataset>\n" % str(self.norm_file)
        
        # Q cut
        xml += "<q_min>%s</q_min>\n" % str(self.q_min)
        xml += "<q_step>%s</q_step>\n" % str(self.q_step)
        
        # Scattering angle
        xml += "<theta>%s</theta>\n" % str(self.theta)
        xml += "<center_pixel>%s</center_pixel>\n" % str(self.center_pixel)
        xml += "<use_center_pixel>%s</use_center_pixel>\n" % str(self.use_center_pixel)
        
        xml += "</RefMData>\n"

        return xml

    def from_xml(self, xml_str):
        self.reset()    
        dom = xml.dom.minidom.parseString(xml_str)
        self.from_xml_element(dom)
        element_list = dom.getElementsByTagName("RefMData")
        if len(element_list)>0:
            instrument_dom = element_list[0]

    def from_xml_element(self, instrument_dom):
        """
            Read in data from XML
            @param xml_str: text to read the data from
        """   
        #Peak from/to pixels
        self.DataPeakPixels = [BaseScriptElement.getIntElement(instrument_dom, "from_peak_pixels"),
                               BaseScriptElement.getIntElement(instrument_dom, "to_peak_pixels")]
        
        
        #low resolution range
        self.data_x_range_flag = BaseScriptElement.getBoolElement(instrument_dom, "data_x_range_flag",
                                                                  default=DataSets.data_x_range_flag)
        
        self.data_x_range = [BaseScriptElement.getIntElement(instrument_dom, "x_min_pixel"),
                             BaseScriptElement.getIntElement(instrument_dom, "x_max_pixel")]
        
        self.norm_x_range_flag = BaseScriptElement.getBoolElement(instrument_dom, "norm_x_range_flag",
                                                                  default=DataSets.norm_x_range_flag)

        self.norm_x_range = [BaseScriptElement.getIntElement(instrument_dom, "norm_x_min"),
                             BaseScriptElement.getIntElement(instrument_dom, "norm_x_max")]
        
        #discrete selection string
        self.DataPeakDiscreteSelection = BaseScriptElement.getStringElement(instrument_dom, "peak_discrete_selection")
        
        #background flag
        self.DataBackgroundFlag = BaseScriptElement.getBoolElement(instrument_dom,
                                                                   "background_flag",
                                                                   default=DataSets.DataBackgroundFlag)

        #background from/to pixels
        self.DataBackgroundRoi = [BaseScriptElement.getIntElement(instrument_dom, "back_roi1_from"),
                                  BaseScriptElement.getIntElement(instrument_dom, "back_roi1_to"),
                                  BaseScriptElement.getIntElement(instrument_dom, "back_roi2_from"),
                                  BaseScriptElement.getIntElement(instrument_dom, "back_roi2_to")]

        #from TOF and to TOF
        self.crop_TOF_range = BaseScriptElement.getBoolElement(instrument_dom, "crop_tof",
                                                               default=DataSets.crop_TOF_range)
        self.DataTofRange = [BaseScriptElement.getFloatElement(instrument_dom, "from_tof_range"),
                             BaseScriptElement.getFloatElement(instrument_dom, "to_tof_range")]

        self.data_files = BaseScriptElement.getIntList(instrument_dom, "data_sets")
            
        #with or without norm 
        self.NormFlag = BaseScriptElement.getBoolElement(instrument_dom, "norm_flag",
                                                         default=DataSets.NormFlag)
        
        #Peak from/to pixels
        self.NormPeakPixels = [BaseScriptElement.getIntElement(instrument_dom, "norm_from_peak_pixels"),
                               BaseScriptElement.getIntElement(instrument_dom, "norm_to_peak_pixels")]

        #background flag
        self.NormBackgroundFlag = BaseScriptElement.getBoolElement(instrument_dom, 
                                                                   "norm_background_flag", 
                                                                   default=DataSets.NormBackgroundFlag)
        
        #background from/to pixels
        self.NormBackgroundRoi = [BaseScriptElement.getIntElement(instrument_dom, "norm_from_back_pixels"),
                                  BaseScriptElement.getIntElement(instrument_dom, "norm_to_back_pixels")]
        
        self.norm_file = BaseScriptElement.getIntElement(instrument_dom, "norm_dataset")
    
        # Q cut
        self.q_min = BaseScriptElement.getFloatElement(instrument_dom, "q_min", default=DataSets.q_min)    
        self.q_step = BaseScriptElement.getFloatElement(instrument_dom, "q_step", default=DataSets.q_step)
    
        # scattering angle
        self.theta = BaseScriptElement.getFloatElement(instrument_dom, "theta", default=DataSets.theta)
        self.center_pixel = BaseScriptElement.getFloatElement(instrument_dom, "center_pixel", default=DataSets.center_pixel)
        self.use_center_pixel = BaseScriptElement.getBoolElement(instrument_dom, 
                                                                 "use_center_pixel", 
                                                                 default=DataSets.use_center_pixel)
        
    def reset(self):
        """
            Reset state
        """
        self.DataBackgroundFlag = DataSets.DataBackgroundFlag
        self.DataBackgroundRoi = DataSets.DataBackgroundRoi
        self.DataPeakPixels = DataSets.DataPeakPixels
        self.DataTofRange = DataSets.DataTofRange
        self.crop_TOF_range = DataSets.crop_TOF_range
        self.data_files = DataSets.data_files
        
        self.NormFlag = DataSets.NormFlag
        self.NormBackgroundFlag = DataSets.NormBackgroundFlag
        self.NormBackgroundRoi = DataSets.NormBackgroundRoi
        self.NormPeakPixels = DataSets.NormPeakPixels
        self.norm_file = DataSets.norm_file
        self.data_x_range_flag = DataSets.data_x_range_flag
        self.data_x_range = DataSets.data_x_range
        self.norm_x_range_flag = DataSets.norm_x_range_flag
        self.norm_x_range = DataSets.norm_x_range
        
        # Q range
        self.q_min = DataSets.q_min
        self.q_step = DataSets.q_step
        
        # Scattering angle
        self.theta = DataSets.theta
        self.center_pixel = DataSets.center_pixel
        self.use_center_pixel = DataSets.use_center_pixel
