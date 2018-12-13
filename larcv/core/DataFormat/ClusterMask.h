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
     Points also stored along with bounding box. Intended use is to store ground-truth for
       instance segmentation. This means collecting the following into one object
       1) bounding box for an object
       2) points inside the object
       3) the class of that object (typedef int InteractionID_t)
     
     One application we keep in mind, is the need to have to generate this truth information
       for a crop, or sub-region, of an image. This means, the bounding box and mask have to 
       be modified when the crop only contains a portion of the object.
    
  */
  typedef int InteractionID_t;


  class ClusterMask {
  public:
    ClusterMask(); 
    ClusterMask(BBox2D box_in, ImageMeta meta, std::vector<Point2D> pts_in, InteractionID_t type_in);
    ClusterMask(ImageMeta meta_in, std::vector<Point2D> pts_in, InteractionID_t type_in);
    
    ~ClusterMask() {}


    BBox2D box; //Placeholder for bbox
    ImageMeta meta;
    std::vector<Point2D> points_v;
    InteractionID_t type;
    //std::vector<float> _box; // what's the purpose of this variable?
    //std::vector<float> _mask; // only generated if asked for it, we don't store this dense object unless we have to

    double get_fraction_clustered() { return points_v.size() / (box.area());}
    bool check_containment();
    //const std::vector<float>& as_vector_mask() const {return _mask;}
    std::vector<float> as_vector_mask( const larcv::ImageMeta& meta ) const;          // generate mask, given the meta
    std::vector<float> as_vector_mask() const { return as_vector_mask(meta); }; // generate mask, using original meta

    std::vector<float> as_vector_box( const larcv::ImageMeta& meta ) const;
    std::vector<float> as_vector_box() const { return as_vector_box(this->meta); };



  };
  //End of Class

}
#endif
/** @} */ // end of doxygen group
