/**
* \file InfilloutputImages.h
* \ingroup UBImageMod
* @author kmason 
*
* This function is for use after using run_infill_precropped.py
* It takes the output from the Infill Network and returns 2 sets of images
*   The network prediction of just the dead channel
*   The above prediction overlayed on the original wire image
**/

/** \addtogroup UBImageMod

    @{*/
#ifndef __INFILLOUTPUTIMAGES_H__
#define __INFILLOUTPUTIMAGES_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventBBox.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class Infill_outputImages ... these comments are used to generate
     doxygen documentation!
  */
  class InfilloutputImages : public ProcessBase {

  public:

    /// Default constructor
    InfilloutputImages(const std::string name = "InfilloutputImages");

    /// Default destructor
    ~InfilloutputImages() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize() {};
    // ----------------------------------------------------------------------
    // we save data ourselves

    public:

      const larcv::IOManager& getOutputIOMan() const { return *foutIO; };

    protected:

      larcv::IOManager* foutIO;

    // ----------------------------------------------------------------------

  private:

    // config parameters
    std::string _input_wire_producer;
    std::string _input_fill_producer;
    std::string _input_nofill_producer;
    std::string _input_weights_producer;
    std::string _output_infill_producer;
    std::string _output_overlay_producer;
    std::string _output_filename;
    int _verbosity_;

  };

  /**
     \class larcv::InfilloutputImagesFactory
     \brief A concrete factory class for larcv::InfilloutputImages
  */
  class InfilloutputImagesProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    InfilloutputImagesProcessFactory() { ProcessFactory::get().add_factory("InfilloutputImages", this); }
    /// dtor
    ~InfilloutputImagesProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new InfilloutputImages(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
