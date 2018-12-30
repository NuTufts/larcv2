// /**
//  * \file BatchFillerClusterMask.h
//  *
//  * \ingroup ThreadIO
//  *
//  * \brief Class def header for a class BatchFillerClusterMask
//  *
//  * @author kazuhiro
//  */
//
// /** \addtogroup ThreadIO
//
//     @{*/
// #ifndef __BATCHFILLERCLUSTERMASK_H__
// #define __BATCHFILLERCLUSTERMASK_H__
//
// #include "larcv/core/Processor/ProcessFactory.h"
// #include "BatchFillerTemplate.h"
// #include "RandomCropper.h"
// #include "larcv/core/DataFormat/EventClusterMask.h"
// namespace larcv {
//
//   /**
//      \class ProcessBase
//      User defined class BatchFillerClusterMask ... these comments are used to generate
//      doxygen documentation!
//   */
//   class BatchFillerClusterMask : public BatchFillerTemplate<float> {
//
//   public:
//
//     /// Default constructor
//     BatchFillerClusterMask(const std::string name="BatchFillerClusterMask");
//
//     /// Default destructor
//     ~BatchFillerClusterMask(){}
//
//     void configure(const PSet&);
//
//     void initialize();
//
//     bool process(IOManager& mgr);
//
//     void finalize();
//
//     const std::vector<bool>& mirrored() const { return _mirrored; }
//
//   protected:
//
//     void _batch_begin_();
//     void _batch_end_();
//
//   private:
//
//     size_t set_clustermask_size(const EventClusterMask* clustermask_data);
//     void assert_dimension(const EventClusterMask* clustermask_data) const;
//
//     bool _caffe_mode;
//     std::string _clustermask_producer;
//     size_t _rows;
//     size_t _cols;
//     size_t _num_channels;
//     std::vector<size_t> _slice_v;
//     size_t _max_ch;
//     std::vector<size_t> _caffe_idx_to_img_idx;
//     std::vector<size_t> _mirror_caffe_idx_to_img_idx;
//     std::vector<bool>   _mirrored;
//     std::vector<float>  _entry_data;
//     bool _mirror_image;
//     bool _crop_image;
//
//     RandomCropper _cropper;
//
//   };
//
//   /**
//      \class larcv::BatchFillerClusterMaskFactory
//      \brief A concrete factory class for larcv::BatchFillerClusterMask
//   */
//   class BatchFillerClusterMaskProcessFactory : public ProcessFactoryBase {
//   public:
//     /// ctor
//     BatchFillerClusterMaskProcessFactory() { ProcessFactory::get().add_factory("BatchFillerClusterMask",this); }
//     /// dtor
//     ~BatchFillerClusterMaskProcessFactory() {}
//     /// creation method
//     ProcessBase* create(const std::string instance_name) {
//       return new BatchFillerClusterMask(instance_name);
//     }
//   };
//
// }
//
// #endif
// /** @} */ // end of doxygen group
