/**
 * \file UBCropMask.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropMask
 *
 * @author jmills
 *
 * We use crops from UBSplitDetector to adjust and carry forward ClusterMasks
 * for the different ancestors.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPMASK_H__
#define __UBCROPMASK_H__

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventClusterMask.h"
#include "larcv/core/DataFormat/EventBBox.h"

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBCropMask ... these comments are used to generate
     doxygen documentation!
  */
  class UBCropMask : public ProcessBase {

  public:

    /// Default constructor
    UBCropMask(const std::string name = "UBCropMask");

    /// Default destructor
    ~UBCropMask() {}

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

   void AdjustMasksUsingBBox2D( const larcv::EventBBox2D& bbox_vec,
            const std::vector<larcv::ImageMeta>& meta_orig_v,
            const std::vector<larcv::ImageMeta>& meta_crop_v,
            const std::vector<std::vector<larcv::ClusterMask>>& masks_vv,
            int num_calls,
            larcv::EventClusterMask& output_masks
            );

  int get_count_cosmic() {return nCosmic;}
  int get_count_electron() {return nElectron;}
  int get_count_proton() {return nProton;}
  int get_count_neutron() {return nNeutron;}
  int get_count_neutrino() {return nNeutrino;}
  int get_count_other() {return nOther;}




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
    std::string _input_masks_producer;
    std::string _input_cropped_producer;

    std::string _output_adc_producer;
    std::string _output_masks_producer;
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

    int nCosmic;
    int nElectron;
    int nProton;
    int nOther;
    int nNeutrino;
    int nNeutron;
  };

  /**
     \class larcv::UBCropMaskFactory
     \brief A concrete factory class for larcv::UBCropMask
  */
  class UBCropMaskProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBCropMaskProcessFactory() { ProcessFactory::get().add_factory("UBCropMask", this); }
    /// dtor
    ~UBCropMaskProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBCropMask(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
