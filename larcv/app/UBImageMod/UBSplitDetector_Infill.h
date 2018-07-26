/**
 * \file UBSplitDetector_Infill.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBSplitDetector_Infill
 *
 * @author kmason
 *
 * We carve the detector into 3D regions and define vectors of bounding boxes to
 * then produce cropped images.
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBSPLITDETECTOR_INFILL_H__
#define __UBSPLITDETECTOR_INFILL_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "larcv/core/DataFormat/BBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBSplitDetector_Infill ... these comments are used to generate
     doxygen documentation!
  */
  class UBSplitDetector_Infill : public ProcessBase {

  public:

    /// Default constructor
    UBSplitDetector_Infill(const std::string name = "UBSplitDetector_Infill");

    /// Default destructor
    ~UBSplitDetector_Infill() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // algo functions
    // static functions are defined to allow them to be reused in a stand-alone manner
    static std::vector<larcv::BBox2D> defineBoundingBoxFromCropCoords( const std::vector<larcv::Image2D>& img_v,
								       const int box_pixel_width, const int box_pixel_height,
								       const int t1, const int t2,
								       const int u1, const int u2,
								       const int v1, const int v2,
								       const int y1, const int y2);

    static bool cropUsingBBox2D( const std::vector<larcv::BBox2D>& bbox_vec,
				 const std::vector<larcv::Image2D>& img_v,
                                 const std::vector<larcv::Image2D>& labels_v,
                                 const std::vector<larcv::Image2D>& ADC_v,
				 const int y1, const int y2, bool fill_y_image,
				 const float minpixfrac,
                                 larcv::EventImage2D& out_ADC_imgs,
                                 larcv::EventImage2D& out_labels_imgs,
				 larcv::EventImage2D& output_imgs );

    static std::vector<int> defineImageBoundsFromPosZT( const float zwire, const float tmid, const float zwidth, const float dtick,
							const int box_pixel_width, const int box_pixel_height,
							const std::vector<larcv::Image2D>& img_v );




  private:

    // config parameters
    std::string _input_producer;
    std::string _labels_input;
    std::string _ADC_input;
    std::string _output_bbox_producer;
    std::string _output_img_producer;
    std::string _output_labels_producer;
    std::string _output_ADC_producer;
    std::string _output_weights_producer;
    bool _enable_img_crop;
    int _box_pixel_height;
    int _box_pixel_width;
    int _covered_z_width;
    bool _complete_y_crop;
    bool _debug_img;
    int _max_images;
    // randomize cropping parameters
    bool _randomize_crops;
    int  _randomize_attempts;
    float _randomize_minfracpix;
  };

  /**
     \class larcv::UBSplitDetectorFactory
     \brief A concrete factory class for larcv::UBSplitDetector
  */
  class UBSplitDetector_InfillProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBSplitDetector_InfillProcessFactory() { ProcessFactory::get().add_factory("UBSplitDetector_Infill", this); }
    /// dtor
    ~UBSplitDetector_InfillProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBSplitDetector_Infill(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
