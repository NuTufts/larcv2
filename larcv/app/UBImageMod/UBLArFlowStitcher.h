/**
 * \file UBLArFlowStitcher.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBLArFlowStitcher
 *
 * @author twongjirad
 *
 * We use crops from UBSplitDetector to crop out
 * LArFlow images. This does the work we have to do
 * to change the pixel flow variables from one crop into another.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPLARFLOW_H__
#define __UBCROPLARFLOW_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBLArFlowStitcher ... these comments are used to generate
     doxygen documentation!
  */
  class UBLArFlowStitcher : public ProcessBase {

  public:

    /// Default constructor
    UBLArFlowStitcher(const std::string name = "UBLArFlowStitcher");

    /// Default destructor
    ~UBLArFlowStitcher() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();
    
  }


  /**
     \class larcv::UBLArFlowStitcherFactory
     \brief A concrete factory class for larcv::UBLArFlowStitcher
  */
  class UBLArFlowStitcherProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBLArFlowStitcherProcessFactory() { ProcessFactory::get().add_factory("UBLArFlowStitcher", this); }
    /// dtor
    ~UBLArFlowStitcherProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBLArFlowStitcher(instance_name); }
  };

}

#endif
