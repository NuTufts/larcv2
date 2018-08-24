/**
 * \file UBCropLArFlow.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropLArFlow
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
     User defined class UBCropLArFlow ... these comments are used to generate
     doxygen documentation!
  */
  
  class UBCropLArFlow : public ProcessBase {

  public:

    /// Default constructor
    UBCropLArFlow(const std::string name = "UBCropLArFlow");

    /// Default destructor
    ~UBCropLArFlow() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions
    static void make_cropped_flow_images( const int src_plane,
					  const larcv::ImageMeta& srcmeta,
					  const std::vector<const larcv::Image2D* >& croppedadc_v,
					  const std::vector<larcv::Image2D>& srcflow,
					  const std::vector<larcv::Image2D>& srcvisi,
					  const std::vector<float>& thresholds,
					  const bool maxpool, const int row_ds_factor, const int col_ds_factor,					  
					  std::vector<larcv::Image2D>& cropped_flow,
					  std::vector<larcv::Image2D>& cropped_visi,
					  const larcv::logger* log=NULL );
    
    static std::vector<float> check_cropped_images( const int src_plane,
						    const std::vector<const larcv::Image2D*>& croppedadc_v,
						    const std::vector<float>& thresholds,
						    const std::vector<larcv::Image2D>& cropped_flow,
						    const std::vector<larcv::Image2D>& cropped_visi,
						    const bool visualize_flow,				      
						    const larcv::logger* log=NULL, const int verbosity=2 );
    
    static void downsample_crops( const std::vector<larcv::Image2D*>& cropped_adc_v,
				  const std::vector<larcv::Image2D>& cropped_flow_v,
				  const std::vector<larcv::Image2D>& cropped_visi_v,
				  std::vector<larcv::Image2D>& downsampled_adc_v,
				  std::vector<larcv::Image2D>& downsampled_flow_v,
				  std::vector<larcv::Image2D>& downsampled_visi_v );

    void maxPool( const int row_downsample_factor, const int col_downsample_factor,
		  const larcv::Image2D& src_adc, const larcv::Image2D& target_adc,
		  const larcv::Image2D& flow, const larcv::Image2D& visi,
		  const std::vector<float>& thresholds,
		  larcv::Image2D& ds_src_adc, larcv::Image2D& ds_target_adc,
		  larcv::Image2D& ds_flow, larcv::Image2D& ds_visi );
    
    
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
    std::string _input_vis_producer;
    std::string _input_flo_producer;
    std::string _input_cropped_producer;
    std::string _output_adc_producer;
    std::string _output_vis_producer;
    std::string _output_flo_producer;
    std::string _output_meta_producer;    
    std::string _output_filename;
    std::vector<float> _thresholds_v;
    int   _max_images;    
    bool  _check_flow;
    bool  _make_check_image;
    bool  _do_maxpool;
    int   _row_downsample_factor;
    int   _col_downsample_factor;
    bool  _limit_overlap;
    bool  _require_min_goodpixels;
    float _max_overlap_fraction;
    int   _verbosity_;
    bool  _save_output;
    bool  _is_mc;

    // algorithms
    struct FlowOffset {
      FlowOffset( float offset ) { _offset = offset; };
      float _offset;
      float operator()( const float& flo_value ) const { return flo_value+_offset; };      
    };
    struct ModVisibility {
      ModVisibility( float target_xmin, float target_xmax, float col ) { _xmin = target_xmin; _xmax = target_xmax; _col = col; };
      float _xmin;
      float _xmax;
      float _col;
      float operator()( const float& flo_value, const float& vis_value ) {
	// too much branching?
	// could break up intro separate pieces
	if ( flo_value<=-4000 ) return 0.0; // no-flow value in that pixel
	if ( vis_value<1.0 ) return 0.0; // don't check visi not 1.0, this means there was dead region or below thresh
	float targetwire = _col+flo_value;
	if ( targetwire < _xmin || targetwire >=_xmax ) return 0.0;
	return 1.0;
      };
    };
    struct MaskBelowThreshold {
      MaskBelowThreshold( float threshold, float maskvalue ) { _threshold = threshold; _maskvalue = maskvalue; };
      float _threshold;
      float _maskvalue;
      float operator()( const float& adc_value, const float& pixvalue ) {
	if ( adc_value<_threshold )
	  return _maskvalue;
	else
	  return pixvalue;
      };
    };
    
    static int _check_img_counter;
    static const float _NO_FLOW_VALUE_;

    // color palette for visualization
    static void setBWRPalette();
    static void setRainbowPalette();
    static int* _colors;
  };

  /**
     \class larcv::UBCropLArFlowFactory
     \brief A concrete factory class for larcv::UBCropLArFlow
  */
  class UBCropLArFlowProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBCropLArFlowProcessFactory() { ProcessFactory::get().add_factory("UBCropLArFlow", this); }
    /// dtor
    ~UBCropLArFlowProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBCropLArFlow(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group

