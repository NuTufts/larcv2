/**
 * \file BatchFillerBBox2D.h
 *
 * \ingroup ThreadIO
 * 
 * \brief Class def header for a class BatchFillerBBox2D
 *
 * @author twongjirad
 */

/** \addtogroup ThreadIO

    @{*/
#ifndef __BATCHFILLERBBOX2D_H__
#define __BATCHFILLERBBOX2D_H__

#include "larcv/core/Processor/ProcessFactory.h"
#include "BatchFillerTemplate.h"

namespace larcv {

  /**
     \class BatchFillerBBox2D
     From bbox2d trees, fills bbounding boxes.
  */
  class BatchFillerBBox2D : public BatchFillerTemplate<float> {

  public:
    
    /// Default constructor
    BatchFillerBBox2D(const std::string name="BatchFillerBBox2D");
    
    /// Default destructor
    ~BatchFillerBBox2D(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void _batch_begin_();

    void _batch_end_();

    void finalize();

  private:

    // we only load boxes from a certain projection (view).
    // create multiple instances of batchfillerbox2d to load
    // bounding boxes from different views
    ProjectionID_t _projid;

    // the name of the ROOT tree storing bbox2d instances
    std::string _bbox2d_producer;

    // entry_data is where we store values before passing it into the
    // data members of the parent BatchFillerTemplate class
    // we have to provide unrolled values.
    // so values ordered like: center-x,center-y,width,height,center-x,...    
    std::vector< float > _entry_data; 

  };

  /**
     \class larcv::BatchFillerBBox2DFactory
     \brief A concrete factory class for larcv::BatchFillerBBox2D
  */
  class BatchFillerBBox2DProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    BatchFillerBBox2DProcessFactory() { ProcessFactory::get().add_factory("BatchFillerBBox2D",this); }
    /// dtor
    ~BatchFillerBBox2DProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new BatchFillerBBox2D(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group 
