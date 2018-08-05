/**
 * \file UBCropSegment.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropSegment
 *
 * @author jmills
 *
 * We use crops from UBSplitDetector to crop out
 * Segment Labeled images.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPSEGMENT_H__
#define __UBCROPSEGMENT_H__

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventBBox.h"

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBCropSegment ... these comments are used to generate
     doxygen documentation!
  */
  class UBCropSegment : public ProcessBase {

  public:

    /// Default constructor
    UBCropSegment(const std::string name = "UBCropSegment");

    /// Default destructor
    ~UBCropSegment() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions
    void cropUsingBBox2D( const larcv::EventBBox2D& bbox_vec,
         const std::vector<larcv::Image2D>& img_v,
         int num_calls,
         larcv::EventImage2D& output_imgs );



    // ----------------------------------------------------------------------
    // we save data ourselves

  public:

    const larcv::IOManager& getOutputIOMan() const { return *foutIO; };

  protected:

    larcv::IOManager* foutIO;

    // ----------------------------------------------------------------------

  private:

    // config parameters
    std::string _input_bbox_producer;
    std::string _input_adc_producer;
    std::string _input_labels_producer;
    std::string _input_cropped_producer;

    std::string _output_adc_producer;
    std::string _output_labels_producer;
    std::string _output_weights_producer;

    std::string _output_meta_producer;
    std::string _output_filename;
    int _max_images;
    std::vector<float> _thresholds_v;
    bool _check_flow;
    bool _make_check_image;


    bool _limit_overlap;
    float _max_overlap_fraction;
    int _verbosity_;


    static int _check_img_counter;
    static const float _NO_FLOW_VALUE_;

    int _background_offset;
    int _shower_offset;
    int _end_offset;
    int _track_offset;
  };

  /**
     \class larcv::UBCropSegmentFactory
     \brief A concrete factory class for larcv::UBCropSegment
  */
  class UBCropSegmentProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBCropSegmentProcessFactory() { ProcessFactory::get().add_factory("UBCropSegment", this); }
    /// dtor
    ~UBCropSegmentProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBCropSegment(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
