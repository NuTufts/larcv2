#ifndef __UBCROPINFILL_CXX__
#define __UBCROPINFILL_CXX__

#include <sstream>

#include "UBCropInfill.h"
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

  static UBCropInfillProcessFactory __global_UBCropInfillProcessFactory__;
  int   UBCropInfill::_check_img_counter = 0;
  const float UBCropInfill::_NO_FLOW_VALUE_ = -4000;

  UBCropInfill::UBCropInfill(const std::string name)
    : ProcessBase(name)
  {}

  void UBCropInfill::configure(const PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");
    _input_bbox_producer    = cfg.get<std::string>("InputBBoxProducer");
    _input_adc_producer     = cfg.get<std::string>("InputADCProducer");
    _input_wire_producer     = cfg.get<std::string>("InputWireProducer");
    _input_labels_producer     = cfg.get<std::string>("InputLabelsProducer");
    _output_adc_producer    = cfg.get<std::string>("OutputCroppedADCProducer");
    _output_wire_producer    = cfg.get<std::string>("OutputCroppedWireProducer");
    _output_labels_producer    = cfg.get<std::string>("OutputCroppedLabelsProducer");
    _output_weights_producer   = cfg.get<std::string>("OutputCroppedWeightsProducer");
    _output_filename        = cfg.get<std::string>("OutputFilename");

    _max_images             = cfg.get<int>("MaxImages",-1);
    _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) );


    // sparsity requirement: prevent cropping over the same regions
    _limit_overlap        = cfg.get<bool>("LimitOverlap",false);
    _max_overlap_fraction = cfg.get<float>("MaxOverlapFraction", 0.5 );

    // debug options
    _check_flow             = cfg.get<bool>("CheckFlow",false);      // output to screen, checks of cropped images
    _make_check_image       = cfg.get<bool>("MakeCheckImage",false); // dump png of image checks
    if ( _make_check_image )
      gStyle->SetOptStat(0);

    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();

  }

  void UBCropInfill::initialize()
  {}

  bool UBCropInfill::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data

    // input ADC
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_adc_v = ev_in_adc->image2d_array();

    // input wire
    auto ev_in_wire  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_wire_producer));
    if (!ev_in_wire) {
      LARCV_CRITICAL() << "No Input wire Image2D found with a name: " << _input_wire_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_wire_v = ev_in_wire->image2d_array();

    // input labels
    auto ev_in_labels  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_labels_producer));
    if (!ev_in_labels) {
      LARCV_CRITICAL() << "No Input labels Image2D found with a name: " << _input_labels_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_labels_v = ev_in_labels->image2d_array();

    // input BBox
    auto ev_in_bbox  = (larcv::EventBBox2D *)(mgr.get_data("bbox2d", _input_bbox_producer));
    if (!ev_in_bbox) {
      LARCV_CRITICAL() << "No Input BBox2D found with a name: " << _input_bbox_producer << std::endl;
      throw larbys();
    }

    // ----------------------------------------------------------------

    // Output containers
    larcv::EventImage2D* ev_out_adc  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_adc_producer);
    larcv::EventImage2D* ev_out_wire  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_wire_producer);
    larcv::EventImage2D* ev_out_labels  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_labels_producer);
    larcv::EventImage2D* ev_out_weights  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_weights_producer);

    ev_out_adc->clear();
    ev_out_wire->clear();
    ev_out_labels->clear();
    ev_out_weights->clear();

    // ----------------------------------------------------------------

    // Run, subrun, event
    int run    = ev_in_adc->run();
    int subrun = ev_in_adc->subrun();
    int event  = ev_in_adc->event();
    int nboxes = ev_in_bbox->size()/3;
    int n=0;
    for (int nimages=0;nimages < nboxes ;nimages++){
      std::vector<larcv::Image2D> cropped_wire;
      std::vector<larcv::Image2D> cropped_adc;
      std::vector<larcv::Image2D> cropped_labels;

      cropUsingBBox2D(
                *ev_in_bbox,
                img_wire_v,
                cropped_wire,
                n );
      cropUsingBBox2D(
              *ev_in_bbox,
               img_adc_v,
               cropped_adc,
                n );
      cropUsingBBox2D(
               *ev_in_bbox,
               img_labels_v,
               cropped_labels,
                n );
      n++;
      // check the quality of the crop
      bool passes_check_filter = true;
      passes_check_filter = CheckImages(cropped_adc,cropped_labels);
      if ( passes_check_filter ) {

        	ev_out_wire->emplace( std::move(cropped_wire) );
        	ev_out_adc->emplace( std::move(cropped_adc) );
        	ev_out_labels->emplace( std::move(cropped_labels) );

          //----------------------------------------------------------
          //creating weighted outputs
          const std::vector< larcv::Image2D >& labels2_v = ev_out_labels->image2d_array();
          const std::vector< larcv::Image2D >& ADC2_v = ev_out_adc->image2d_array();

          std::vector< larcv::Image2D > weights_v;
          for (auto const& img : labels2_v){
            Image2D copy = img;
            weights_v.emplace_back(std::move(copy));
          }

          float ncharge = 0;
          float nempty = 0;
          size_t planes = labels2_v.size();
          for ( unsigned int plane = 0 ; plane < planes; plane++){
            size_t cols = labels2_v[plane].meta().cols();
            size_t rows = labels2_v[plane].meta().rows();
              for ( unsigned int col = 0; col < cols; col++){
                for ( unsigned int row = 0; row < rows; row++){
                   if (labels2_v[plane].pixel(row,col) == 1){
                      if (ADC2_v[plane].pixel(row,col) >0 ) ncharge++;
                      else nempty++;
                   }
                }
             }
             for ( unsigned int col = 0; col <cols; col++){
                for ( unsigned int row = 0; row < rows; row++){
                   if (labels2_v[plane].pixel(row,col) ==1 ){
                      if (ADC2_v[plane].pixel(row,col) >0){
                         weights_v[plane].set_pixel(row,col, 1.0/ncharge);
                      }
                      else weights_v[plane].set_pixel(row,col, 1.0/nempty);

                    }
                 }
             }
             ncharge = 0;
             nempty = 0;
          }
          ev_out_weights->emplace(std::move(weights_v));

          //----------------------------------------------------------

        	foutIO->set_id( run, subrun, 100*event+nimages );
        	foutIO->save_entry();

      }

      if ( _max_images>0 && nimages>_max_images) break;

   }
    return true;
  }

  void UBCropInfill::cropUsingBBox2D(
           larcv::EventBBox2D& bbox_vec,
					 const std::vector<larcv::Image2D>& img_v,
					 std::vector<larcv::Image2D>& output_imgs,
            int nimages ) {
    // inputs
    // ------
    // bbox_v, vector of bounding boxes for (u,v,y)
    // img_v, source adc images
    // nimages, number of the bbox to use
    // outputs
    // --------
    // output_imgs, cropped output image2d instances filled into eventimage2d container


    // get bounding boxes
    const larcv::BBox2D& bbox_u = bbox_vec[nimages*3+0];
    const larcv::BBox2D& bbox_v = bbox_vec[nimages*3+1];
    const larcv::BBox2D& bbox_y = bbox_vec[nimages*3+2];
    //std::cout << bbox_vec.size() << std::endl;
    //crop y
    larcv::Image2D crop_yp = img_v[2].crop( bbox_y );
    output_imgs.emplace_back( std::move(crop_yp) );

    larcv::Image2D crop_up = img_v[0].crop( bbox_u );
    output_imgs.emplace_back( std::move(crop_up) );

    larcv::Image2D crop_vp = img_v[1].crop( bbox_v );
    output_imgs.emplace_back( std::move(crop_vp) );

    return;
  }

 bool UBCropInfill::CheckImages(std::vector<larcv::Image2D>& cropped_adc,
    std::vector<larcv::Image2D>& cropped_labels){
      float occupied = 0.;
      bool saveimg = true;
      int minpixfrac=10;
      if ( minpixfrac>=0 ) {
        for (int plane=0; plane<3; plane++){
        occupied=0;
        const larcv::Image2D& ADC = cropped_adc[plane];
        const larcv::Image2D& labels = cropped_labels[plane];
        for (int col=0; col<(int)ADC.meta().cols(); col++) {
          for (int row=0; row<(int)ADC.meta().rows(); row++) {
              if ( (ADC.pixel(row,col)*labels.pixel(row,col)) > 0 )
              occupied+=1.0;
           }
         }
        if ( occupied <= 0 ) saveimg = false;
        if ( !saveimg )
          {return false;}
        else
            {std::cout << "Pixels occupied in dead regions of accepted image: " << occupied << std::endl; }
        }
      }
      return true;
      }


  void UBCropInfill::finalize()
  {
    foutIO->finalize();
  }

}
#endif
