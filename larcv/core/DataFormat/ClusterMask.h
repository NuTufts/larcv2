/**
 * \file ClusterMask.h
 *
 * \ingroup DataFormat
 *
 * \brief Class def header for a class ClusterMask
 *
 * @author JoshuaMills
 */

/** \addtogroup DataFormat

    @{*/
#ifndef LARCV_CLUSTERMASK_H
#define LARCV_CLUSTERMASK_H

#include <iostream>
#include <cmath>
#include "Point.h"
#include "BBox.h"
#include "ImageMeta.h"
namespace larcv {

  /**
     \class ClusterMask
     Simple 2D point struct (unit of "x" and "y" are not defined here and app specific)
  */
  typedef int InteractionID_t;


  class ClusterMask {
  public:
    ClusterMask();
    ClusterMask(BBox2D box_in, ImageMeta meta, std::vector<Point2D> pts_in, InteractionID_t type_in);

    ~ClusterMask() {}


    BBox2D box; //Placeholder for bbox
    ImageMeta meta;
    std::vector<Point2D> points_v;
    InteractionID_t type;
    std::vector<float> _box;
    std::vector<float> _mask;

    double get_fraction_clustered() {return points_v.size() / (box.area());}
    bool check_containment() ;
    const std::vector<float>& as_vector_mask() const {return _mask;}
    const std::vector<float>& as_vector_box() const {return _box;}



  };
  //End of Class

}
#endif
/** @} */ // end of doxygen group
