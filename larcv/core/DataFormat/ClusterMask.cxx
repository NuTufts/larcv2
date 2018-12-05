#ifndef LARCV_CLUSTERMASK_CXX
#define LARCV_CLUSTERMASK_CXX

#include "ClusterMask.h"

namespace larcv {

  ClusterMask::ClusterMask()
  : box(0,0,0,0,kINVALID_PROJECTIONID) , meta(ImageMeta()), points_v(0,Point2D(0,0))
  {
    type = 0;
    _box = {(float) meta.row(box.min_x()), (float) meta.col(box.min_y()), (float) meta.row(box.max_x()), (float) meta.col(box.max_y()), (float) type};
  }


  ClusterMask::ClusterMask(BBox2D box_in, ImageMeta meta_in, std::vector<Point2D> pts_in, InteractionID_t type_in)
  :  box(box_in), meta(meta_in), points_v(pts_in), type(type_in)
  {
    _box = {(float) meta_in.col(box.min_x()), (float) meta_in.row(box.min_y()), (float) meta_in.col(box.max_x()), (float) meta_in.row(box.max_y()), (float) type};
    _mask = std::vector<float>( (box.height()/meta_in.pixel_height()+1) * (box.width()/meta_in.pixel_width()+1), 0.0);
    for (Point2D pt : points_v){
      _mask[(int)(pt.x - meta_in.col(box.min_x())) * box.height()/meta_in.pixel_height()+1 + pt.y-meta_in.row(box.min_y())] = 1.0;
    }
  }


  bool ClusterMask::check_containment() {
    for (Point2D pt : points_v) {
      if ((pt.x >= 0 && pt.y >= 0 && pt.x <= box.width() && pt.y <= box.height())==false) {return false;}
    }
    return true;
  }
}
#endif
