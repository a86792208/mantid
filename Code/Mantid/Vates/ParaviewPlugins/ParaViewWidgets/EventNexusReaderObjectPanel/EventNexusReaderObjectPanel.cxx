#include "EventNexusReaderObjectPanel.h"
#include "pqObjectPanelInterface.h"
#include "pqPropertyManager.h"
#include "pqNamedWidgets.h"
#include "vtkSMStringVectorProperty.h"
#include "vtkSMDoubleVectorProperty.h"
#include <QLayout>
#include <QMessageBox>
#include <QCheckBox>
#include <QGridLayout>
#include "GeometryWidget.h"
#include "ThresholdRangeWidget.h"
#include "MantidGeometry/MDGeometry/MDGeometryXMLParser.h"
#include "MantidVatesAPI/SynchronisingGeometryPresenter.h"

using namespace Mantid::VATES;
using namespace Mantid::Geometry;

EventNexusReaderObjectPanel::EventNexusReaderObjectPanel(pqProxy* pxy, QWidget* p) :
pqAutoGeneratedObjectPanel(pxy, p), m_cachedMinThreshold(0), m_cachedMaxThreshold(0), m_geometryXMLString(""), m_geometryWidget(NULL), m_thresholdWidget(NULL)
{
  //Auto generated widgets are replaced by Custom Widgets. Autogenerated ones need to be removed.
  removeAutoGeneratedWidgets();
}

/// Event handler for framework event.
void EventNexusReaderObjectPanel::updateInformationAndDomains()
{
  this->proxy()->UpdatePropertyInformation();
  QGridLayout* gLayout = dynamic_cast<QGridLayout*>(this->layout());
  try
  {
    this->constructThresholdRanges(gLayout);
    this->constructGeometry(gLayout);
  }
  catch(std::exception& ex)
  {
    UNUSED_ARG(ex);
    QMessageBox::information(NULL, "Setup Not possible.",
      "Could not interpret metadata. Are you using a rebinning source? Check field data.");
  }
}

void EventNexusReaderObjectPanel::constructGeometry(QGridLayout* gLayout)
{
  vtkSMStringVectorProperty* inputGeometryProperty = vtkSMStringVectorProperty::SafeDownCast(
    this->proxy()->GetProperty("InputGeometryXML"));

  std::string geometryXMLString = inputGeometryProperty->GetElement(0);

  if(m_geometryXMLString != geometryXMLString) //Only attempt to reconstruct the geometry widget if the xml has changed.
  {
    MDGeometryXMLParser xmlParser(geometryXMLString);
    xmlParser.execute();

    //Empty geometry widget added to layout.
    if(m_geometryWidget != NULL)
    {
      this->layout()->removeWidget(m_geometryWidget);
      delete m_geometryWidget;
    }

    // Construct custom widget instance.
    m_geometryWidget = new GeometryWidget(new SynchronisingGeometryPresenter(xmlParser), false);
    gLayout->addWidget(m_geometryWidget, gLayout->rowCount() + 1, 0, Qt::AlignLeft);

    connect(m_geometryWidget, SIGNAL(ignoreBinChanges()), this, SLOT(ignoreBinChangesListner()));

    // Property used as setter.
    vtkSMProperty * appliedGeometryXML = this->proxy()->GetProperty("AppliedGeometryXML");

    //Hook up geometry change event to listener on filter.
    this->propertyManager()->registerLink(m_geometryWidget, "GeometryXML",
      SIGNAL(valueChanged()), this->proxy(), appliedGeometryXML);

    m_geometryXMLString = geometryXMLString;
  }
}

void EventNexusReaderObjectPanel::constructThresholdRanges(QGridLayout* gLayout)
{
  // Access getter property to extract original input max threshold value.
  vtkSMDoubleVectorProperty* inputMaxThresholdProperty = vtkSMDoubleVectorProperty::SafeDownCast(
    this->proxy()->GetProperty("InputMaxThreshold"));
  double inputMaxThreshold = inputMaxThresholdProperty->GetElement(0);

  // Access getter property to extract original input min threshold value.
  vtkSMDoubleVectorProperty* inputMinThresholdProperty = vtkSMDoubleVectorProperty::SafeDownCast(
    this->proxy()->GetProperty("InputMinThreshold"));
  double inputMinThreshold = inputMinThresholdProperty->GetElement(0);

  if(inputMaxThreshold != m_cachedMaxThreshold || inputMinThreshold != m_cachedMinThreshold)
  {

    if(m_thresholdWidget == NULL)
    {
      m_thresholdWidget = new ThresholdRangeWidget(inputMinThreshold, inputMaxThreshold);
      gLayout->addWidget(m_thresholdWidget, gLayout->rowCount() + 1, 0, Qt::AlignCenter);
    }
    else
    {
      m_thresholdWidget->setMaximum(inputMaxThreshold);
      m_thresholdWidget->setMinimum(inputMinThreshold);
    }

    // Property used as setter
    vtkSMProperty * minThreshold = this->proxy()->GetProperty("MinThreshold");

    // Property used as setter
    vtkSMProperty * maxThreshold = this->proxy()->GetProperty("MaxThreshold");

    // Property used as setter
    vtkSMProperty * rangeStrategy = this->proxy()->GetProperty("ThresholdRangeStrategyIndex");

    // Hook-up events to PV properties.
    this->propertyManager()->registerLink(m_thresholdWidget, "MinSignal",
      SIGNAL(minChanged()), this->proxy(), minThreshold);

    // Hook-up events to PV properties.
    this->propertyManager()->registerLink(m_thresholdWidget, "MaxSignal",
      SIGNAL(maxChanged()), this->proxy(), maxThreshold);

    // Hook-up events to PV properties.
    this->propertyManager()->registerLink(m_thresholdWidget, "ChosenStrategy",
      SIGNAL(chosenStrategyChanged()), this->proxy(), rangeStrategy);

    m_cachedMaxThreshold = inputMaxThreshold;
    m_cachedMinThreshold = inputMinThreshold;
  }
}

/// Direct removal of autogenerated widgets.
void EventNexusReaderObjectPanel::removeAutoGeneratedWidgets()
{
  popWidget(); // Autogenerated Geometry QLineEdit
  popWidget(); // Autogenerated Geometry QLabel
  popWidget(); // Autogenerated Max threshold QLineEdit
  popWidget(); // Autogenerated Max threshold QLabel
  popWidget(); // Autogenerated Min threshold QLineEdit
  popWidget(); // Autogenerated Min threshold QLabel
  popWidget(); // Autogenerated User defined ComboBox.
  popWidget(); // Autogenerated User defined ComboBox.
}

/// Pop widgets off the layout and hide them.
void EventNexusReaderObjectPanel::popWidget()
{
  unsigned int size = layout()->count();
  if(size >= 1)
  {
    //Pop the last widget off the layout and hide it.
    QLayoutItem* pLayoutItem = layout()->itemAt(size - 1);
    QWidget* pWidget = pLayoutItem->widget();
    if (NULL == pWidget)
    {
      throw std::domain_error(
        "Error ::popWidget(). Attempting to pop a non-widget object off the layout!");
    }
    else
    {
      pWidget->setHidden(true);
      this->layout()->removeItem(pLayoutItem);
    }
  }
}

/// Listener for ignore bin changes.
void EventNexusReaderObjectPanel::ignoreBinChangesListner()
{
  QLayoutItem *child;
  unsigned int size = layout()->count();
  for(unsigned int i = 0; i < size; i++)
  {
    child = layout()->itemAt(i);
    QWidget* pWidget = child->widget();
    if (NULL != pWidget) // capability query of layout item.
    {
      QCheckBox* checkBox;
      if((checkBox = dynamic_cast<QCheckBox*>(pWidget)) != NULL) //capability query of widget.
      {
        // Apply Clip check box set to unchecked.
        checkBox->setChecked(false);
        break;
      }
    }
  }
}
