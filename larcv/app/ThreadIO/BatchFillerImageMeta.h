/**
 * \file BatchFillerImageMeta.h
 *
 * \ingroup ThreadIO
 * 
 * \brief Class def header for a class BatchFillerImageMeta
 *
 * @author taritree
 */

/** \addtogroup ThreadIO

    @{*/
#ifndef __BATCHFILLERIMAGEMETA_H__
#define __BATCHFILLERIMAGEMETA_H__

#include "larcv/core/Processor/ProcessFactory.h"
#include "BatchFillerTemplate.h"
#include "RandomCropper.h"
#include "larcv/core/DataFormat/EventImage2D.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class BatchFillerImageMeta ... these comments are used to generate
     doxygen documentation!
  */
  class BatchFillerImageMeta : public BatchFillerTemplate<float> {

  public:
    
    /// Default constructor
    BatchFillerImageMeta(const std::string name="BatchFillerImageMeta");
    
    /// Default destructor
    ~BatchFillerImageMeta(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();


  protected:

    size_t set_image_size(const EventImage2D* image_data);
    
    void _batch_begin_();
    void _batch_end_();

  private:

    std::string _image_producer;
    size_t _num_channels;
    std::vector<size_t> _slice_v;
    size_t _max_ch;
    std::vector<float>  _entry_data;
    
  };

  /**
     \class larcv::BatchFillerImageMetaFactory
     \brief A concrete factory class for larcv::BatchFillerImageMeta
  */
  class BatchFillerImageMetaProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    BatchFillerImageMetaProcessFactory() { ProcessFactory::get().add_factory("BatchFillerImageMeta",this); }
    /// dtor
    ~BatchFillerImageMetaProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) {
      return new BatchFillerImageMeta(instance_name);
    }
  };

}

#endif
/** @} */ // end of doxygen group 

