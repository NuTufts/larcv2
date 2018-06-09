#ifndef __UBLARFLOWSTITCHER_CXX__
#define __UBLARFLOWSTITCHER_CXX__

#include <sstream>

#include "UBLArFlowStitcher.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventMeta.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/ROOTUtil/ROOTUtils.h"

//larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

namespace larcv {

  
  UBLArFlowStitcher::UBLArFlowStitcher(const std::string name)
    : ProcessBase(name)
  {}

  void UBLArFlowStitcher::configure(const PSet& cfg)
  {

    // _verbosity_             = cfg.get<int>("Verbosity");
    // _input_bbox_producer    = cfg.get<std::string>("InputBBoxProducer");
    // _input_adc_producer     = cfg.get<std::string>("InputADCProducer");
    // _input_cropped_producer = cfg.get<std::string>("InputCroppedADCProducer");
    // _input_vis_producer     = cfg.get<std::string>("InputVisiProducer");
    // _input_flo_producer     = cfg.get<std::string>("InputFlowProducer");
    // _output_adc_producer    = cfg.get<std::string>("OutputCroppedADCProducer");
    // _output_vis_producer    = cfg.get<std::string>("OutputCroppedVisiProducer");
    // _output_flo_producer    = cfg.get<std::string>("OutputCroppedFlowProducer");
    // _output_meta_producer   = cfg.get<std::string>("OutputCroppedMetaProducer");    
    // _output_filename        = cfg.get<std::string>("OutputFilename");

    // _max_images             = cfg.get<int>("MaxImages",-1);
    // _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) );

    // // Max pooling. Shrink image by some downsampling factor. must be factor of image size.
    // _do_maxpool             = cfg.get<bool>("DoMaxPool",false);
    // if (_do_maxpool) {
    //   _row_downsample_factor  = cfg.get<int>("RowDownsampleFactor");
    //   _col_downsample_factor  = cfg.get<int>("ColDownsampleFactor");
    // }
    // else {
    //   _row_downsample_factor = -1;
    //   _col_downsample_factor = -1;
    // }

    // // sparsity requirement: prevent cropping over the same regions
    // _limit_overlap        = cfg.get<bool>("LimitOverlap",false);
    // _max_overlap_fraction = cfg.get<float>("MaxOverlapFraction", 0.5 );

    // // debug options
    // _check_flow             = cfg.get<bool>("CheckFlow",false);      // output to screen, checks of cropped images
    // _make_check_image       = cfg.get<bool>("MakeCheckImage",false); // dump png of image checks
    // if ( _make_check_image )
    //   gStyle->SetOptStat(0);
          
  }

  void UBLArFlowStitcher::initialize()
  {}

  bool UBLArFlowStitcher::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data

    // // input ADC
    // auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_adc_producer));
    // if (!ev_in_adc) {
    //   LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
    //   throw larbys();
    // }
    // const std::vector< larcv::Image2D >& img_v = ev_in_adc->image2d_array();

    // // input visibility/matchability
    // auto ev_in_vis  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_vis_producer));
    // if (!ev_in_vis) {
    //   LARCV_CRITICAL() << "No Input VIS Image2D found with a name: " << _input_vis_producer << std::endl;
    //   throw larbys();
    // }
    // const std::vector< larcv::Image2D >& vis_v = ev_in_vis->image2d_array();

    // // input flo
    // auto ev_in_flo  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_flo_producer));
    // if (!ev_in_flo) {
    //   LARCV_CRITICAL() << "No Input flo Image2D found with a name: " << _input_flo_producer << std::endl;
    //   throw larbys();
    // }
    // const std::vector< larcv::Image2D >& flo_v = ev_in_flo->image2d_array();

    // // output BBox
    // auto ev_in_bbox  = (larcv::EventImage2D*)(mgr.get_data("bbox2d", _input_bbox_producer));
    // if (!ev_in_bbox) {
    //   LARCV_CRITICAL() << "No Input BBox2D found with a name: " << _input_bbox_producer << std::endl;
    //   throw larbys();
    // }
    
    return true;
  }
  
  
  void UBLArFlowStitcher::finalize()
  {
  }

  
}
#endif
