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
  {
    _output_initialized = false;
    _psrc_adc_v = NULL;
    _verbosity = 2;
    _output_flo_y2u_producer = "larflow_y2u";
    _output_vis_y2u_producer = "larvisi_y2u";
  }

  void UBLArFlowStitcher::configure(const PSet& cfg)
  {
    _verbosity               = cfg.get<int>("Verbosity",2);
    _output_flo_y2u_producer = cfg.get<std::string>("OutputY2UFlowProducer","larflow_y2u");
    _output_vis_y2u_producer = cfg.get<std::string>("OutputY2UVisiProducer","larvisi_y2u");
  }

  void UBLArFlowStitcher::initialize()
  {}

  bool UBLArFlowStitcher::process(IOManager& mgr)
  {
    // we simply save the output image

    // output flo
    auto ev_out_flo_y2u  = (larcv::EventImage2D*)(mgr.get_data("image2d", _output_flo_y2u_producer));
    auto ev_out_vis_y2u  = (larcv::EventImage2D*)(mgr.get_data("image2d", _output_vis_y2u_producer));    

    ev_out_flo_y2u->emplace( std::move(_output_y2u[2]) );    // follow first in u-plane pos
    ev_out_flo_y2u->emplace( std::move(_output_y2u[0]) );    // flow prediction in y-plane pos
    
    ev_out_vis_y2u->emplace( std::move(_output_y2u[1]) );
    
    _output_initialized = false; // reset
    
    return true;
  }

  void UBLArFlowStitcher::initializeOutput( const std::vector<larcv::Image2D>& src_adc_v ) {

    // clear output
    if ( !_output_y2u.empty() ) {
      _output_y2u.clear();
    }

    // Y2U flow
    larcv::Image2D y2u_flow( src_adc_v.at(2).meta() );    // flow prediction, back in whole-image space
    larcv::Image2D y2u_visi( src_adc_v.at(2).meta() );    // visibility, in whole-image space
    larcv::Image2D y2u_follow( src_adc_v.at(0).meta() );  // follow the flow (ideally creates image similar to u-plane adc image)
    larcv::Image2D y2u_trusted( src_adc_v.at(0).meta() ); // marked if filled by trusted region (center of yimage where u/v overlap by construction)
    _output_y2u.emplace_back( std::move(y2u_flow) );
    _output_y2u.emplace_back( std::move(y2u_visi) );
    _output_y2u.emplace_back( std::move(y2u_follow) );
    _output_y2u.emplace_back( std::move(y2u_trusted) );    

    _output_initialized = true;
  }

  void UBLArFlowStitcher::setInputImage( const std::vector<larcv::Image2D>& src_adc_v ) {

    // save pointer to source adc
    _psrc_adc_v = &src_adc_v;

  };
    
  void UBLArFlowStitcher::clearOutput() {
    for ( auto& img : _output_y2u )
      img.paint(0.0);
  }

  void UBLArFlowStitcher::setupEvent( const std::vector<larcv::Image2D>& src_adc_v ) {
    if ( !_output_initialized )
      initializeOutput( src_adc_v );

    clearOutput();
    setInputImage( src_adc_v );
  }

  void UBLArFlowStitcher::insertFlowSubimage( const larcv::Image2D& flow_predict, const larcv::ImageMeta& tar_meta ) {

    // flow_predict: larflow scores, in the y-image. a subimage of output_y2u
    // flow_target:  the u-plane adc image.

    if ( !_output_initialized ) {
      throw std::runtime_error( "output not yet initialized. call setupEvent( const std::vector<larcv::Image2D>& ) for each event first." );
    }
    // add check for event number

    larcv::Image2D& out_flow = _output_y2u[0];
    larcv::Image2D& out_visi = _output_y2u[1];
    larcv::Image2D& out_follow = _output_y2u[2];
    larcv::Image2D& out_trust  = _output_y2u[3];        
    const larcv::ImageMeta& out_meta = out_flow.meta();
    
    int src_min_c = out_meta.col( flow_predict.meta().min_x() );
    int src_min_r = out_meta.row( flow_predict.meta().min_y() );
    
    for ( int r=0; r<(int)flow_predict.meta().rows(); r++ ) {
      for ( int c=0; c<(int)flow_predict.meta().cols(); c++ ) {

	int src_c = src_min_c+c;
	int src_r = src_min_r+r;
	float src_adc = _psrc_adc_v->at(2).pixel(src_r,src_c);
	
	if ( src_adc<10.0 )
	  continue;

	if ( out_trust.pixel(src_r,src_c)>0 )
	  continue; // already filled by trusted region
	
	int flo = std::round(flow_predict.pixel(r,c)); // prediction
	// std::cout << "src(r,c)=(" << src_r << "," << src_c << ") adc=" << _psrc_adc_v->at(2).pixel(src_r,src_c) << " "
	// 	  << "flowsub(r,c)=(" << r << "," << c << ") "
	// 	  << "flo=" << flo
	// 	  << std::endl;
	
	int tar_col = c+flo; // target column
	if ( tar_col>(int)tar_meta.cols() || tar_col<0 ) {
	  //std::cout << "  flo out of bounds. target col=" << tar_col << std::endl;
	  continue;
	}
	
	float tar_wire = tar_meta.pos_x( tar_col ); // target wire
	int out_col = out_meta.col( tar_wire ); // column in output image
	
	float final_flow = out_col-src_c; // output flow
	//std::cout << "  ok value. target_wire=" << tar_wire << " final_flow=" << final_flow << std::endl;
	out_flow.set_pixel( src_r, src_c, final_flow );
	out_follow.set_pixel( src_r, out_col, src_adc );
	if ( c>=261 || c<=571 )
	  out_trust.set_pixel( src_r, src_c, 1.0 );	
      }
    }
    
  }
  
  void UBLArFlowStitcher::finalize()
  {
  }

  
}
#endif
