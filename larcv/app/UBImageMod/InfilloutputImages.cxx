//#ifndef __INFILLOUTPUTIMAGES_CXX__
#define __INFILLOUTPUTIMAGES_CXX__

#include <sstream>

#include "InfilloutputImages.h"
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

  static InfilloutputImagesProcessFactory __global_InfilloutputImagesProcessFactory__;

  InfilloutputImages::InfilloutputImages(const std::string name)
    : ProcessBase(name)
  {}

  void InfilloutputImages::configure(const PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");
    _input_wire_producer     = cfg.get<std::string>("InputWireProducer");
    _input_nofill_producer     = cfg.get<std::string>("InputNoFillProducer");
    _input_fill_producer     = cfg.get<std::string>("InputFillProducer");
    _input_weights_producer     = cfg.get<std::string>("InputWeightsProducer");
    _output_infill_producer    = cfg.get<std::string>("OutputInfillProducer");
    _output_overlay_producer    = cfg.get<std::string>("OutputOverlayProducer");
    _output_filename        = cfg.get<std::string>("OutputFilename");


    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();

  }

  void InfilloutputImages::initialize()
  {}

  bool InfilloutputImages::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data

    // input wire
    auto ev_in_wire  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_wire_producer));
    if (!ev_in_wire) {
      LARCV_CRITICAL() << "No Input wire Image2D found with a name: " << _input_wire_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_wire_v = ev_in_wire->image2d_array();

    // input weights
    auto ev_in_weights  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_weights_producer));
    if (!ev_in_weights) {
      LARCV_CRITICAL() << "No Input weights Image2D found with a name: " << _input_weights_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_weights_v = ev_in_weights->image2d_array();

    //input nofill
    auto ev_in_nofill  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_nofill_producer));
    if (!ev_in_nofill) {
      LARCV_CRITICAL() << "No Input wire Image2D found with a name: " << _input_nofill_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_nofill_v = ev_in_nofill->image2d_array();

    //input fill
    auto ev_in_fill  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_fill_producer));
    if (!ev_in_fill) {
      LARCV_CRITICAL() << "No Input wire Image2D found with a name: " << _input_fill_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_fill_v = ev_in_fill->image2d_array();

    // ----------------------------------------------------------------

    // Output containers
    larcv::EventImage2D* ev_out_overlay  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_overlay_producer);
    larcv::EventImage2D* ev_out_infill  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_infill_producer);

    ev_out_overlay->clear();
    ev_out_infill->clear();
    //std::cout<<"made outputs" <<std::endl;
    // ----------------------------------------------------------------

   std::vector<larcv::Image2D> temp_overlay;
   std::vector<larcv::Image2D> temp_infill;

   const std::vector<larcv::Image2D>& wire = ev_in_wire->image2d_array();

   for (auto const& img1 : wire){
      Image2D copy = img1;
      temp_overlay.emplace_back(std::move(copy));
      }
   for (auto const& img2 : wire){
      Image2D copy2 = img2;
      temp_infill.emplace_back(std::move(copy2));
      }

  const std::vector<larcv::Image2D>& nofill = ev_in_nofill->image2d_array();
  const std::vector<larcv::Image2D>& fill = ev_in_fill->image2d_array();
  const std::vector<larcv::Image2D>& weights = ev_in_weights->image2d_array();
  //loop through each pixel
  size_t columns = wire[0].meta().cols();
  size_t rows = wire[0].meta().rows();
  for (unsigned int col=0; col<columns; col++) {
      for (unsigned int row=0; row<rows; row++) {
            float weight_value = (weights[0].pixel(row,col));
            float wire_value = (wire[0].pixel(row,col));
               if (weight_value == 0){
                  temp_infill[0].set_pixel(row,col,wire_value);
                  temp_overlay[0].set_pixel(row,col,0);
               }
               else{
                  float fill_value = (fill[0].pixel(row,col));
                  float nofill_value = (nofill[0].pixel(row,col));
                  if (fill_value > nofill_value){
                    temp_infill[0].set_pixel(row,col,1);
                    temp_overlay[0].set_pixel(row,col,1);
                  }
                  else{
                    temp_infill[0].set_pixel(row,col,0);
                    temp_overlay[0].set_pixel(row,col,0);
                  }
               }
        }
   }

    ev_out_overlay->emplace (std::move(temp_overlay) );
    ev_out_infill->emplace( std::move(temp_infill) );
  	foutIO->save_entry();
    return true;
  }

 } 
