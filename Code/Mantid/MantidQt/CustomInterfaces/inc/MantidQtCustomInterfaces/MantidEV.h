#ifndef MANTIDQTCUSTOMINTERFACES_MANTID_EV_H_
#define MANTIDQTCUSTOMINTERFACES_MANTID_EV_H_

#include <QtCore/QtCore>
#include <QtGui/QWidget>
#include <QActionGroup>
#include <QRunnable>
#include <MantidKernel/System.h>

#include "ui_MantidEV.h"
#include "MantidEVWorker.h"
#include "MantidQtAPI/UserSubWindow.h"

namespace MantidQt
{
namespace CustomInterfaces
{

/// Local class to load file and convert to MD in a Non-Qt thread.
class RunLoadAndConvertToMD : public QRunnable
{
  public:

  /// Constructor just saves the info needed by the run() method
  RunLoadAndConvertToMD(       MantidEVWorker * worker, 
                         const std::string    & file_name,
                         const std::string    & ev_ws_name,
                         const std::string    & md_ws_name );

  /// Calls worker->loadAndConvertToMD from a separate thread
  void run();

  private:
    MantidEVWorker * worker;
    std::string      file_name;
    std::string      ev_ws_name;
    std::string      md_ws_name;
};


/// Local class to run FindPeaks in a Non-Qt thread.
class RunFindPeaks : public QRunnable
{
  public:

  /// Constructor just saves the info needed by the run() method
  RunFindPeaks(       MantidEVWorker * worker,
                const std::string    & md_ws_name,
                const std::string    & peaks_ws_name,
                      double           max_abc,
                      size_t           num_to_find,
                      double           min_intensity );

  /// Calls worker->findPeaks from a separate thread
  void run();

  private:
    MantidEVWorker * worker;
    std::string      md_ws_name;
    std::string      peaks_ws_name;
    double           max_abc;
    size_t           num_to_find;
    double           min_intensity;
};


/// Local class to run IntegratePeaksMD in a Non-Qt thread.
class RunSphereIntegrate : public QRunnable
{
  public:

  /// Constructor just saves the info needed by the run() method
  RunSphereIntegrate(       MantidEVWorker * worker,
                      const std::string    & peaks_ws_name,
                      const std::string    & event_ws_name,
                            double           peak_radius,
                            double           inner_radius,
                            double           outer_radius,
                            bool             integrate_edge );

  /// Calls worker->sphereIntegrate from a separate thread
  void run();

  private:
    MantidEVWorker * worker;
    std::string      peaks_ws_name;
    std::string      event_ws_name;
    double           peak_radius; 
    double           inner_radius; 
    double           outer_radius; 
    bool             integrate_edge; 
};


/// Local class to run PeakIntegration in a Non-Qt thread.
class RunFitIntegrate : public QRunnable
{
  public:

  /// Constructor just saves the info needed by the run() method
  RunFitIntegrate(       MantidEVWorker * worker,
                   const std::string    & peaks_ws_name,
                   const std::string    & event_ws_name,
                   const std::string    & rebin_params,
                         size_t           n_bad_edge_pix,
                         bool             use_ikeda_carpenter );

  /// Calls worker->fitIntegrate from a separate thread
  void run();

  private:
    MantidEVWorker * worker;
    std::string      peaks_ws_name;
    std::string      event_ws_name;
    std::string      rebin_params;
    size_t           n_bad_edge_pix;
    bool             use_ikeda_carpenter;
};


/// Local class to run ellipsoidIntegrate in a Non-Qt thread.
class RunEllipsoidIntegrate : public QRunnable
{
  public:

  /// Constructor just saves the info needed by the run() method
  RunEllipsoidIntegrate(       MantidEVWorker * worker,
                         const std::string    & peaks_ws_name,
                         const std::string    & event_ws_name,
                               double           region_radius,
                               bool             specify_size,
                               double           peak_size, 
                               double           inner_size, 
                               double           outer_size );

  /// Calls worker->ellipsoidIntegrate from a separate thread
  void run();

  private:
    MantidEVWorker * worker;
    std::string      peaks_ws_name;
    std::string      event_ws_name;
    double           region_radius;
    bool             specify_size;
    double           peak_size;
    double           inner_size;
    double           outer_size;
};


class MantidEV : public API::UserSubWindow
{
  Q_OBJECT

public:

  MantidEV(QWidget *parent = 0);
  ~MantidEV();

  /// The name of the interface as registered into the factory
  static std::string name() { return "SCD Event Data Reduction"; }


private slots:

  /// Slot for the select workspace tab's Apply button 
  void selectWorkspace_slot();

  /// Slot for the Browse button for loading an event file 
  void loadEventFile_slot();

  /// Slot for the find peaks tab's Apply button 
  void findPeaks_slot();

  /// Slot for the find UB tab's Apply button 
  void findUB_slot();

  /// Slot for the Browse button for choosing a matrix file 
  void loadUB_slot();

  /// Slot for the choose cell tab's Apply button 
  void chooseCell_slot();

  /// Slot for the change HKL tab's Apply button 
  void changeHKL_slot();

  /// Slot for the integrate tab's Apply button 
  void integratePeaks_slot();

  /// Slot to enable/disable the Load Event File controls
  void setEnabledLoadEventFileParams_slot( bool on );

  /// Slot to enable/disable the find peaks controls
  void setEnabledFindPeaksParams_slot( bool on );

  /// Slot to enable/disable the find UB using FFT controls
  void setEnabledFindUBFFTParams_slot( bool on );

  /// Slot to enable/disable the load UB controls
  void setEnabledLoadUBParams_slot( bool on );

  /// Slot to enable/disable the optimize goniometer angles controls
  void setEnabledMaxOptimizeDegrees_slot();

  /// Slot to enable/disable the index peaks controls
  void setEnabledIndexParams_slot( bool on );

  /// Slot to enable/disable the show conventional cell controls
  void setEnabledShowCellsParams_slot( bool on );

  /// Slot to enable/disable the select cell based on cell type controls
  void setEnabledSetCellTypeParams_slot(bool on );

  /// Slot to enable/disable the select cell based on cell number controls
  void setEnabledSetCellFormParams_slot(bool on );

  /// Slot to enable/disable the sphere integration controls
  void setEnabledSphereIntParams_slot( bool on );

  /// Slot to enable/disable the 2D fitting integration controls
  void setEnabledFitIntParams_slot( bool on );

  /// Slot to enable/disable the ellipsoidal integration controls
  void setEnabledEllipseIntParams_slot( bool on );

  /// Slot to enable/disable the ellipse size options controls
  void setEnabledEllipseSizeOptions_slot();


private:
  /// super class pure virtual method we MUST implement
  virtual void initLayout();

  /// Utility method to display an error message
  void errorMessage( const std::string message );

  /// Utility method to parse a double value in a string
  bool getDouble( std::string str, double & value );

  /// Utility method to get a double value from a QtLineEdit widget 
  bool getDouble( QLineEdit *ledt, double &value );

  /// Utility method to get a positive double from a QtLineEdit widget 
  bool getPositiveDouble( QLineEdit *ledt, double &value );

  /// Utility method to get a positive integer value from a QtLineEdit widget 
  bool getPositiveInt( QLineEdit *ledt, size_t &value );


  Ui::MantidEV   m_uiForm;     ///< The form generated by Qt Designer

  MantidEVWorker *worker;      /// class that uses Mantid algorithms
                               /// to do the actual work

  std::string  last_UB_file;   /// filename of last UB file that was
                               /// loaded, if any.

  std::string  last_event_file;/// filename of last event file that
                               /// was loaded, if any.

  QThreadPool  *m_thread_pool; /// thread pool with only one thread, to 
                               /// allow running precisely one operation 
                               /// at a time in a separate thread.
};

} // namespace CustomInterfaces
} // namespace MantidQt

#endif //MANTIDQTCUSTOMINTERFACES_MANTID_EV_H_