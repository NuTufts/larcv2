#ifndef __BatchFillerBBox2D_CXX__
#define __BatchFillerBBox2D_CXX__

#include "BatchFillerBBox2D.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include <random>
#include <array>
#include <sstream>
#include <exception>

namespace larcv {

  static BatchFillerBBox2DProcessFactory __global_BatchFillerBBox2DProcessFactory__;

  BatchFillerBBox2D::BatchFillerBBox2D(const std::string name)
    : BatchFillerTemplate<float>(name)
  {}

  void BatchFillerBBox2D::configure(const PSet& cfg)
  {
    _bbox2d_producer = cfg.get<std::string>("BBox2DProducer");
    _projid = (ProjectionID_t)cfg.get<int>("ProjectionID");
    _convert_xy_to_pix = cfg.get<bool>("ConvertXYtoPixel",false);
    _imageproducer = cfg.get<std::string>("ImageProducerForMeta","");
    _maxboxes = cfg.get<int>("MaxNumBoxes");
    _include_class = cfg.get<bool>("IncludeClass",true); // includes class index as last entry in array
    std::string format = cfg.get<std::string>("Format","MinMax");

    if (format=="MinMax")
      _format = kMinMax;
    else if (format=="Center")
      _format = kCenter;
    else {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << "::" << __FUNCTION__ << ":: 'Format' parameter must be either 'MinMax' or 'Center'" << std::endl;
      throw std::runtime_error( ss.str() );      
    }
      
  }

  void BatchFillerBBox2D::initialize()
  {}

  void BatchFillerBBox2D::_batch_begin_()
  {
    
    
    if (!batch_data().dim().empty() && (int)(batch_size()) != batch_data().dim().front()) {
      auto dim = batch_data().dim();
      dim[0] = batch_size();
      this->set_dim(dim);
    }

    if ( _include_class )
      _entry_data.resize(batch_size()*5*_maxboxes,0);
    else
      _entry_data.resize(batch_size()*4*_maxboxes,0);
  }

  void BatchFillerBBox2D::_batch_end_()
  {
  }

  void BatchFillerBBox2D::finalize()
  {}

  bool BatchFillerBBox2D::process(IOManager & mgr)
  {

    LARCV_DEBUG() << "start" << std::endl;

    // get the bboxes for the entry
    auto const& event_bbox2d = mgr.get_data<larcv::EventBBox2D>(_bbox2d_producer);

    // if we need to convert coordinate systems
    //  from (x,y) to (row,col)
    //  we need the image2d event container as well
    larcv::EventImage2D* ev_img2d = NULL;
    std::vector< const larcv::ImageMeta* > meta_v;
    if ( _convert_xy_to_pix ) {
      ev_img2d = (larcv::EventImage2D*)mgr.get_data("image2d",_imageproducer);

      // we need the metas for the projections
      // whats't he max id?
      int maxid = 0 ;
      for ( auto const& img : ev_img2d->image2d_array() ) {
	if ( maxid<img.meta().id() )
	  maxid = img.meta().id();
      }
      meta_v.resize( maxid+1, 0 );
      // store pointers to the metas
      for ( auto const& img : ev_img2d->image2d_array() ) {
	if ( !meta_v[img.meta().id()] )
	  meta_v[img.meta().id()] = &(img.meta());
      }
    }

    // if slow, one thing to try is to allocate this only once
    // by keeping as member data
    std::vector< std::array<float,5> > tempstore;
    tempstore.reserve( event_bbox2d.size() );

    for ( auto const& bbox2d : event_bbox2d ) {
      if ( bbox2d.id()!=_projid )
	continue;
      
      const ImageMeta* pmeta = meta_v[bbox2d.id()];
      
      std::array<float,5> xybh;
      xybh[0] = bbox2d.center_x();
      xybh[1] = bbox2d.center_y();
      xybh[2] = bbox2d.width();
      xybh[3] = bbox2d.height();
      xybh[4] = 0;      
      if ( _include_class )
	xybh[4] = 1;

      if  (_convert_xy_to_pix ) {
	xybh[0] = pmeta->col( xybh[0] );
	xybh[1] = pmeta->row( xybh[1] );
	xybh[2] /= pmeta->pixel_width();
	xybh[3] /= pmeta->pixel_height();
      }

      if ( _format==kMinMax ) {
	float x_min = xybh[0]-0.5*xybh[2];
	float y_min = xybh[1]-0.5*xybh[3];
	float x_max = xybh[0]+0.5*xybh[2];
	float y_max = xybh[1]+0.5*xybh[3];
	xybh[0] = x_min;
	xybh[1] = y_min;
	xybh[2] = x_max;
	xybh[3] = y_max;
      }
      
      tempstore.push_back( xybh );
    }
      
    int nboxes = tempstore.size();
    LARCV_DEBUG() << "Number of boxes in event: " << nboxes << std::endl;    
    if ( nboxes>_maxboxes )
      nboxes = _maxboxes;
    
    if (batch_data().dim().empty()) {
      LARCV_DEBUG() << "Setting batch size" << std::endl;
      std::vector<int> dim(4);
      dim[0] = batch_size();
      dim[1] = _maxboxes;
      dim[2] = 1;
      if ( _include_class )
	dim[3] = 5;
      else
	dim[3] = 4;
      
      this->set_dim(dim);
    }
    int stride = 5;
    if ( _include_class )  {
      _entry_data.resize(5*_maxboxes,0);
      stride = 5;
    }
    else {
      _entry_data.resize(4*_maxboxes,0);
      stride = 4;
    }
    for (auto& v : _entry_data) v = 0.;
	  
    for ( int ibox=0; ibox<nboxes; ibox++ ) {
      auto& bbox2d = tempstore[ibox];
      for (int i=0; i<4; i++) 
	_entry_data[stride*ibox+i] = bbox2d[i];
      if (_include_class)
	_entry_data[stride*ibox+4] = bbox2d[4];
    }
    this->set_entry_data(_entry_data);

    LARCV_DEBUG() << "end" << std::endl;    
    return true;
  }

}
#endif
