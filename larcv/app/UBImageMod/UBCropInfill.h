/**
 * \file UBCropInfill.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropInfill
 *
 * @author twongjirad
 *
 * We use crops from UBSplitDetector to crop out
 * Infill images. This does the work we have to do
 * to change the pixel flow variables from one crop into another.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPINFILL_H__
#define __UBCROPINFILL_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventBBox.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class UBCropInfill ... these comments are used to generate
     doxygen documentation!
  */
  class UBCropInfill : public ProcessBase {

  public:

    /// Default constructor
    UBCropInfill(const std::string name = "UBCropInfill");

    /// Default destructor
    ~UBCropInfill() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions
    static  void cropUsingBBox2D(
               larcv::EventBBox2D& bbox_vec,
    					 const std::vector<larcv::Image2D>& img_v,
    					 std::vector<larcv::Image2D>& output_imgs,
               int nimages );
    static bool CheckImages(
              std::vector<larcv::Image2D>& cropped_adc,
              std::vector<larcv::Image2D>& cropped_labels );

    /*static void make_cropped_images( const int src_plane,
    						const larcv::ImageMeta& srcmeta,
    						std::vector<larcv::Image2D>& croppedwire_v,
                const std::vector<larcv::Image2D>& img_adc_v,
                const std::vector<larcv::Image2D>& img_labels_v,
    						const std::vector<float>& thresholds,
    						std::vector<larcv::Image2D>& cropped_adc,
                std::vector<larcv::Image2D>& cropped_labels);*/

    /*static std::vector<float> check_cropped_images( const int src_plane,
						    const std::vector<larcv::Image2D>& cropped_adc_v,
						    const std::vector<float>& thresholds,
						    const std::vector<larcv::Image2D>& cropped_flow,
						    const std::vector<larcv::Image2D>& cropped_visi,
						    const bool visualize_flow,
						    const larcv::logger* log=NULL, const int verbosity=2 );
   */

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
    std::string _input_wire_producer;
    std::string _input_labels_producer;
    std::string _input_adc_producer;
    //std::string _input_cropped_producer;
    std::string _output_wire_producer;
    std::string _output_labels_producer;
    std::string _output_adc_producer;
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
  };

  /**
     \class larcv::UBCropInfillFactory
     \brief A concrete factory class for larcv::UBCropInfill
  */
  class UBCropInfillProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBCropInfillProcessFactory() { ProcessFactory::get().add_factory("UBCropInfill", this); }
    /// dtor
    ~UBCropInfillProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBCropInfill(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
