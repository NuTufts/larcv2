/**
 * \file UBLArFlowStitcher.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBLArFlowStitcher
 *
 * @author twongjirad
 *
 * We use crops from UBSplitDetector to crop out
 * LArFlow images. This does the work we have to do
 * to change the pixel flow variables from one crop into another.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBLARFLOWSTITCHER_H__
#define __UBLARFLOWSTITCHER_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBLArFlowStitcher ... these comments are used to generate
     doxygen documentation!
  */
  class UBLArFlowStitcher : public ProcessBase {

  public:

    /// Default constructor
    UBLArFlowStitcher(const std::string name = "UBLArFlowStitcher");

    /// Default destructor
    ~UBLArFlowStitcher() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    void setupEvent( const std::vector<larcv::Image2D>& src_adc_v );
    void insertFlowSubimage( const larcv::Image2D& flow_predict, const larcv::ImageMeta& flow_target );
    

  protected:

    bool _output_initialized;
    int  _verbosity;
    std::string _output_flo_y2u_producer;
    std::string _output_vis_y2u_producer;
    
    std::vector< larcv::Image2D > _output_y2u; // [0] flow, [1] visi
    const std::vector<larcv::Image2D>* _psrc_adc_v; // pointer to source ADC image

    void setInputImage( const std::vector<larcv::Image2D>& src_adc_v );    
    void initializeOutput( const std::vector<larcv::Image2D>& img_v ); // set the output image size
    void clearOutput();
    
  };


  /**
     \class larcv::UBLArFlowStitcherFactory
     \brief A concrete factory class for larcv::UBLArFlowStitcher
  */
  class UBLArFlowStitcherProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBLArFlowStitcherProcessFactory() { ProcessFactory::get().add_factory("UBLArFlowStitcher", this); }
    /// dtor
    ~UBLArFlowStitcherProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBLArFlowStitcher(instance_name); }
  };

}

#endif
