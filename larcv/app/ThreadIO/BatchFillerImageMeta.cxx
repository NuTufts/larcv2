#ifndef __BatchFillerImageMeta_CXX__
#define __BatchFillerImageMeta_CXX__

#include "BatchFillerImageMeta.h"

#include <random>

namespace larcv {

  static BatchFillerImageMetaProcessFactory __global_BatchFillerImageMetaProcessFactory__;

  BatchFillerImageMeta::BatchFillerImageMeta(const std::string name)
    : BatchFillerTemplate<float>(name)
    , _slice_v()
    , _max_ch(0)
  {}

  void BatchFillerImageMeta::configure(const PSet& cfg)
  {
    LARCV_DEBUG() << "start" << std::endl;
    _image_producer = cfg.get<std::string>("ImageProducer");
    LARCV_DEBUG() << "done" << std::endl;
  }

  void BatchFillerImageMeta::initialize()
  {}

  void BatchFillerImageMeta::_batch_begin_()
  {
    if (!batch_data().dim().empty() && (int)(batch_size()) != batch_data().dim().front()) {
      auto dim = batch_data().dim();
      dim[0] = batch_size();
      this->set_dim(dim);
    }
  }

  void BatchFillerImageMeta::_batch_end_()
  {
    if (logger().level() <= msg::kINFO) {
      LARCV_INFO() << "Total data size: " << batch_data().data_size() << std::endl;
    }
  }

  void BatchFillerImageMeta::finalize()
  { _entry_data.clear(); }

  size_t BatchFillerImageMeta::set_image_size(const EventImage2D* image_data)
  {
    auto const& image_v = image_data->image2d_array();
    if (image_v.empty()) {
      LARCV_CRITICAL() << "Input image is empty!" << std::endl;
      throw larbys();
    }
    if (_slice_v.empty()) {
      // store: xmin,xmax,ymin,ymax
      _slice_v.resize(image_v.size());
      for (size_t i = 0; i < _slice_v.size(); ++i) _slice_v.at(i) = i;
    }

    _num_channels = _slice_v.size();
    _max_ch = 0;
    for (auto const& v : _slice_v) if (_max_ch < v) _max_ch = v;

    if (image_v.size() <= _max_ch) {
      LARCV_CRITICAL() << "Requested slice max channel (" << _max_ch
                       << ") exceeds available # of channels in the input image" << std::endl;
      throw larbys();
    }
    return (4 * _num_channels);
  }

  bool BatchFillerImageMeta::process(IOManager& mgr)
  {
    LARCV_DEBUG() << "start" << std::endl;
    auto image_data = (EventImage2D*)(mgr.get_data("image2d", _image_producer));
    if (!image_data) {
      LARCV_CRITICAL() << "Could not locate image data w/ producer name " << _image_producer << std::endl;
      throw larbys();
    }
    if (image_data->image2d_array().empty()) {
      LARCV_CRITICAL() << "Image data w/ producer name " << _image_producer << " is empty!" << std::endl;
      throw larbys();
    }
    // one time operation: get image dimension
    if (batch_data().dim().empty()) {
      this->set_image_size(image_data);

      std::vector<int> dim;
      dim.resize(4);
      dim[0] = batch_size();
      dim[1] = _num_channels;
      dim[2] = 1;
      dim[3] = 4;
      this->set_dim(dim);
    }

    if (_entry_data.size() != batch_data().entry_data_size())
      _entry_data.resize(batch_data().entry_data_size(), 0.);

    for (auto& v : _entry_data) v = 0.;


    auto const& image_v = ((EventImage2D*)image_data)->image2d_array();

    for (size_t ch = 0; ch < _num_channels; ++ch) {

      size_t input_ch = _slice_v.at(ch);

      auto const& input_meta  = image_v[input_ch].meta();

      _entry_data.at( ch*4 + 0 ) = input_meta.min_x();
      _entry_data.at( ch*4 + 1 ) = input_meta.min_y();
      _entry_data.at( ch*4 + 2 ) = input_meta.max_x();
      _entry_data.at( ch*4 + 3 ) = input_meta.max_y();
      
    }

    // record the entry data
    LARCV_INFO() << "Inserting entry data of size " << _entry_data.size() << std::endl;
    set_entry_data(_entry_data);

    return true;
  }

}
#endif
