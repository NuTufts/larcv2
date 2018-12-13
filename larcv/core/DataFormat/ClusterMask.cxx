#ifndef LARCV_CLUSTERMASK_CXX
#define LARCV_CLUSTERMASK_CXX

#include "ClusterMask.h"

namespace larcv {

  ClusterMask::ClusterMask()
  : box(0,0,0,0,kINVALID_PROJECTIONID) , meta(ImageMeta()), points_v(0,Point2D(0,0))
  {
    type = 0;
    //_box = {(float) meta.row(box.min_x()), (float) meta.col(box.min_y()), (float) meta.row(box.max_x()), (float) meta.col(box.max_y()), (float) type};
  }


  ClusterMask::ClusterMask(BBox2D box_in, ImageMeta meta_in, std::vector<Point2D> pts_in, InteractionID_t type_in)
  :  box(box_in), meta(meta_in), points_v(pts_in), type(type_in)
  {
    // constructor
    //
    // meta_in: the meta which defines how a points (x,y) values relate to the (col,row) position for an image.
    //          represents the meta the points originally embedded in
    // box_in: bounding box in units of (x,y) [not necessarily (col,row) unless meta defined as such]
    // pts_in: (x,y) points
    // type_id: some integer denoting the class
    
    //_box = {(float) box.min_x(), (float) box.min_y(), (float) box.max_x(), (float) box.max_y(), (float) type};
    // note from taritree -- seems to be segfaulting on occassion
    //_mask.reserve( points_v.size() );
    // _mask = std::vector<float>( (box.height()/meta_in.pixel_height()+1) * (box.width()/meta_in.pixel_width()+1), 0.0);
    // for (Point2D pt : points_v){
    //   _mask[(int)(pt.x - meta_in.col(box.min_x())) * box.height()/meta_in.pixel_height()+1 + pt.y-meta_in.row(box.min_y())] = 1.0;
    // }
  }

  ClusterMask::ClusterMask(ImageMeta meta_in, std::vector<Point2D> pts_in, InteractionID_t type_in)
  :  meta(meta_in), points_v(pts_in), type(type_in)
  {
    // constructor (we autogenerate bounding box to be minimal one around pts)
    //
    // meta_in: the meta which defines how a points (x,y) values relate to the (col,row) position for an image.
    //          represents the meta the points originally embedded in (doesn't seem necessary)
    // pts_in: (x,y) points
    // type_id: some integer denoting the class
    //
    // fills in 'box' data member

    if ( pts_in.size()==0 ) {
      throw larbys( "no points provided for ClusterMask" );
    }

    // see box with first point
    double ptmin[2] = { pts_in.front().x, pts_in.front().y };
    double ptmax[2] = { pts_in.front().x, pts_in.front().y };
    
    for ( auto const& pt : pts_in ) {
      if ( pt.x < ptmin[0] ) ptmin[0] = pt.x;
      if ( pt.y < ptmin[1] ) ptmin[1] = pt.y;      
      if ( pt.x > ptmax[0] ) ptmax[0] = pt.x;
      if ( pt.y > ptmax[1] ) ptmax[1] = pt.y;      
    }

    box = BBox2D( ptmin[0], ptmin[1], ptmax[0], ptmax[1], meta.id() );
    
  }
  
  
  std::vector<float> ClusterMask::as_vector_mask( const larcv::ImageMeta& meta ) const {
    std::vector< float > mask( meta.rows()*meta.cols(), 0 );
    for ( auto const& pt : points_v ) {
      if ( meta.contains( pt.y, pt.x ) ) {
	int col = meta.col(pt.x);
	int row = meta.row(pt.y);
	mask[ meta.index( row,col ) ] = 1;
      }
    }
    return mask;
  }
  

  std::vector<float> ClusterMask::as_vector_box( const larcv::ImageMeta& meta ) const {
    // purpose: output adjusted bounding box coordinates, including region overlapping with provided meta
    // inputs
    // ------
    //  meta: defines the image we want to embedd this box into
    // 
    // outputs
    // -------
    //  5-element vector<float>: boundingbox corners + type
    larcv::BBox2D overlap = box.overlap(meta);
    if ( overlap.area()<=0 )
      return std::vector<float>(5,0);

    // overlap probably not exactly divisible by pixel height nor width.
    // we adjust the box for this (output box snaps to pixel grid defined by embedding meta)
    // we accomplish that by first finding (col,row) of given corner in embedding image
    //   then ask for the value of the edges for those col and row values thus
    //   providing the snapped box

    larcv::BBox2D snapped( meta.pos_x( meta.col(overlap.min_x()) ),
			   meta.pos_y( meta.row(overlap.min_y()) ),
			   meta.pos_x( meta.col(overlap.max_x()) ),
			   meta.pos_y( meta.row(overlap.max_y()) ) );

    std::vector<float> bbox_v = { (float) snapped.min_x(),
				  (float) snapped.min_y(),
				  (float) snapped.max_x(),
				  (float) snapped.max_y(),
				  (float) type };
    return bbox_v;
  }

  
  bool ClusterMask::check_containment() {
    // returns truth if all points within the bounding box
    for (Point2D pt : points_v) {
      if ( !box.contains( pt ) )
	return false;
      //if ((pt.x >= 0 && pt.y >= 0 && pt.x <= box.width() && pt.y <= box.height())==false) {return false;}
    }
    return true;
  }
}
#endif
