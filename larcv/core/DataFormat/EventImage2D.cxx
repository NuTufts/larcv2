#ifndef EVENTIMAGE2D_CXX
#define EVENTIMAGE2D_CXX

#include "EventImage2D.h"
#include "larcv/core/Base/larbys.h"

namespace larcv {

  /// Global larcv::SBClusterFactory to register ClusterAlgoFactory
  static EventImage2DFactory __global_EventImage2DFactory__;

  void EventImage2D::clear()
  {
    EventBase::clear();
    _image_v.clear();
  }

  const Image2D& EventImage2D::at(ImageIndex_t id) const
  {
    if( id >= _image_v.size() ) throw larbys("Invalid request (ImageIndex_t out-o-range)!");
    return _image_v[id];
  }

  void EventImage2D::append(const Image2D& img)
  {
    _image_v.push_back(img);
    _image_v.back().index((ImageIndex_t)(_image_v.size()-1));
  }

  void EventImage2D::emplace(Image2D&& img)
  {
    _image_v.emplace_back(std::move(img));
    _image_v.back().index((ImageIndex_t)(_image_v.size()-1));
  }

  void EventImage2D::emplace(std::vector<larcv::Image2D>&& image_v)
  {
    _image_v = std::move(image_v);
    for(size_t i=0; i<_image_v.size(); ++i) _image_v[i].index((ImageIndex_t)i);
  }

  void EventImage2D::paintall( float value ) {
    // calls paint on all stored images in container
    // intention is to be able to zero out all images
    for ( auto &img : _image_v ) {
      img.paint(value);
    }
  }

  Image2D& EventImage2D::modimgat(ImageIndex_t id) {
    if( id >= _image_v.size() ) throw larbys("Invalid request (ImageIndex_t out-o-range)!");
    return _image_v[id];
  }

  void EventImage2D::reserve(size_t s) {
    _image_v.reserve(s);
  }
  
}

#endif
