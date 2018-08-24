#ifndef __UBCROPLARFLOW_CXX__
#define __UBCROPLARFLOW_CXX__

#include <sstream>
#include <algorithm>

#include "UBCropLArFlow.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventMeta.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/ROOTUtil/ROOTUtils.h"

//larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

namespace larcv {

  static UBCropLArFlowProcessFactory __global_UBCropLArFlowProcessFactory__;
  int   UBCropLArFlow::_check_img_counter = 0;
  bool  UBCropLArFlow::_fusevector = false;
  const float UBCropLArFlow::_NO_FLOW_VALUE_ = -4000;
  int* UBCropLArFlow::_colors = NULL;
  
  UBCropLArFlow::UBCropLArFlow(const std::string name)
    : ProcessBase(name)
  {}

  void UBCropLArFlow::configure(const PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");
    _input_adc_producer     = cfg.get<std::string>("InputADCProducer"); // whole image
    _input_bbox_producer    = cfg.get<std::string>("InputBBoxProducer");  // from UBSplitDetector
    _input_cropped_producer = cfg.get<std::string>("InputCroppedADCProducer"); // from UBSplitDetector
    _input_vis_producer     = cfg.get<std::string>("InputVisiProducer"); // whole image with visi truth
    _input_flo_producer     = cfg.get<std::string>("InputFlowProducer"); // whole image with flow truth
    _output_adc_producer    = cfg.get<std::string>("OutputCroppedADCProducer");  
    _output_vis_producer    = cfg.get<std::string>("OutputCroppedVisiProducer");
    _output_flo_producer    = cfg.get<std::string>("OutputCroppedFlowProducer");
    _output_meta_producer   = cfg.get<std::string>("OutputCroppedMetaProducer");    
    _output_filename        = cfg.get<std::string>("OutputFilename"); // creates new larcv file to store cropped images, if SaveOuput=true
    _is_mc                  = cfg.get<bool>("IsMC"); // if we expect to have truth visi and flow
    UBCropLArFlow::_fusevector = cfg.get<bool>("UseVectorizedCode",false); // optimization (experimental)

    _max_images             = cfg.get<int>("MaxImages",-1); // maximum number of images to save
    _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) ); // ADC thresholds for each plane

    // Max pooling. Shrink image by some downsampling factor. must be factor of image size.
    _do_maxpool             = cfg.get<bool>("DoMaxPool",false);
    if (_do_maxpool) {
      _row_downsample_factor  = cfg.get<int>("RowDownsampleFactor");
      _col_downsample_factor  = cfg.get<int>("ColDownsampleFactor");
    }
    else {
      _row_downsample_factor = -1;
      _col_downsample_factor = -1;
    }

    // minimum pixel requirements
    _require_min_goodpixels = cfg.get<bool>("RequireMinGoodPixels"); // minimum number of pixels in Y-plane that need to be above threshold
    
    // sparsity requirement: prevent cropping over the same regions
    _limit_overlap        = cfg.get<bool>("LimitOverlap",false); // prevents random crops from overlapping by some fraction
    _max_overlap_fraction = cfg.get<float>("MaxOverlapFraction", 0.5 ); // fraction random crops can overlap

    // debug options
    _check_flow             = cfg.get<bool>("CheckFlow",false);      // output to screen, checks of cropped images
    _make_check_image       = cfg.get<bool>("MakeCheckImage",false); // dump png of image checks
    if ( _make_check_image )
      gStyle->SetOptStat(0);

    // save output
    _save_output = cfg.get<bool>("SaveOutput"); // save cropped images to _output_filename
    
    // output file
    if ( _save_output ) {
      foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
      foutIO->set_out_file( _output_filename );
      foutIO->initialize();
    }
    else {
      foutIO = NULL;
    }
    
  }

  void UBCropLArFlow::initialize()
  {}

  bool UBCropLArFlow::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data

    // input ADC: Whole image
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_v = ev_in_adc->as_vector();

    // input visibility/matchability (whole image)
    larcv::EventImage2D* ev_in_vis = NULL;
    const std::vector< larcv::Image2D >* vis_v = NULL;    
    if ( _is_mc ) {
      ev_in_vis  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_vis_producer));
      if (!ev_in_vis) {
	LARCV_CRITICAL() << "No Input VIS Image2D found with a name: " << _input_vis_producer << std::endl;
	throw larbys();
      }
      vis_v = &(ev_in_vis->image2d_array());
    }
    
    // input flo (whole image)
    larcv::EventImage2D* ev_in_flo = NULL;
    const std::vector< larcv::Image2D >* flo_v = NULL;
    if ( _is_mc ) {
      ev_in_flo  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_flo_producer));
      if (!ev_in_flo) {
	LARCV_CRITICAL() << "No Input flo Image2D found with a name: " << _input_flo_producer << std::endl;
	throw larbys();
      }
      flo_v = &(ev_in_flo->image2d_array());
    }

    // input BBox (defines crops. from UBSplitDetector)
    auto ev_in_bbox  = (larcv::EventImage2D*)(mgr.get_data("bbox2d", _input_bbox_producer));
    if (!ev_in_bbox) {
      LARCV_CRITICAL() << "No Input BBox2D found with a name: " << _input_bbox_producer << ". make using UBSplitDetector." << std::endl;
      throw larbys();
    }

    // cropped input ADC (from UBSplitDetector)
    auto ev_in_cropped  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_cropped_producer));
    if (!ev_in_cropped) {
      LARCV_CRITICAL() << "No Input Cropped ADC Image2D found with a name: " << _input_cropped_producer << ". make using UBSplitDetector" << std::endl;
      throw larbys();
    }
    std::vector< larcv::Image2D >& cropped_v = ev_in_cropped->as_mod_vector(); // we use mutable, because we might pass objects to output container
    

    // ----------------------------------------------------------------

    // Output ADC containers
    larcv::EventImage2D* ev_out_adc  = NULL;
    larcv::EventImage2D* ev_vis_adc  = NULL;
    larcv::EventImage2D* ev_flo_adc  = NULL;

    if ( _save_output ) {
      ev_out_adc = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_adc_producer);
      if ( _is_mc ) {
	ev_vis_adc = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_vis_producer);
	ev_flo_adc = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_flo_producer);
      }
    }
    else {
      ev_out_adc = (larcv::EventImage2D*)mgr.get_data("image2d",_output_adc_producer);
      if ( _is_mc ) {
	ev_vis_adc = (larcv::EventImage2D*)mgr.get_data("image2d",_output_vis_producer);
	ev_flo_adc = (larcv::EventImage2D*)mgr.get_data("image2d",_output_flo_producer);
      }
    }
    ev_out_adc->clear();
    if ( _is_mc ) {
      ev_vis_adc->clear();
      ev_flo_adc->clear();
    }
    
    // Output Meta containers: meta of each crop
    larcv::EventMeta* ev_meta = NULL;
    if ( _save_output ) {
      ev_meta = (larcv::EventMeta*)foutIO->get_data("meta",_output_meta_producer);
    }
    else {
      ev_meta = (larcv::EventMeta*)mgr.get_data("meta",_output_meta_producer);      
    }
    ev_meta->clear();    
    
    // ----------------------------------------------------------------

    // Run, subrun, event
    int run    = ev_in_adc->run();
    int subrun = ev_in_adc->subrun();
    int event  = ev_in_adc->event();
    
    // crop corresponding flow and visibility images from cropped images
    const int src_plane = 2;
    const larcv::ImageMeta& src_meta = img_v[2].meta();
    int ncrops = cropped_v.size()/3;
    int nsaved = 0;
    //std::cout << "UBCropLArFlow processing " << ncrops << " input crops" << std::endl;

    std::vector<larcv::Image2D> overlap_img; // store an image used to keep track of previously cropped pixels
    if ( _limit_overlap ) {
      larcv::Image2D overlap( img_v[src_plane].meta() );
      overlap.paint(0.0);
      overlap_img.emplace_back( std::move(overlap) );
    }
    
    for (int icrop=0; icrop<ncrops; icrop++) {

      // this is a pointer to crop_v
      std::vector<const larcv::Image2D*> crop_v; // I guess everyone hates bare pointers now. could use reference wrapper?
      std::vector<larcv::Image2D*> mod_crop_v;
      for (int i=0; i<3; i++) {
	crop_v.push_back( &(cropped_v.at( 3*icrop+i )) );
	mod_crop_v.push_back( &(cropped_v.at( 3*icrop+i )) );	
      }

      // if limiting overlap, check that we are not overlapping too many pixels
      if ( _limit_overlap ) {
	const larcv::ImageMeta& cropped_src_meta = crop_v[src_plane]->meta();
	float npix_overlap = 0.0;
	int crop_r_start = src_meta.row( cropped_src_meta.min_y() );
	int crop_c_start = src_meta.col( cropped_src_meta.min_x() );
	for (int r=crop_r_start; r<crop_r_start+(int)cropped_src_meta.rows(); r++) {
	  for (int c=crop_c_start; c<crop_c_start+(int)cropped_src_meta.cols(); c++) {
	    if ( overlap_img[0].pixel(r,c)>0 )
	      npix_overlap+=1.0;
	  }
	}
	float frac_overlap = npix_overlap/float(cropped_src_meta.rows()*cropped_src_meta.cols());
	if ( frac_overlap>_max_overlap_fraction ) {
	  LARCV_INFO() << "Skipping overlapping image. Frac overlap=" << frac_overlap << "." << std::endl;
	  continue;
	}
	//std::cout << "Overlap fraction: " << frac_overlap << std::endl;
      }

      // Crop of Truth Labels (only if using MC file)
      bool passes_check_filter = true;
      std::vector<larcv::Image2D> cropped_flow; // stores cropped flow image set for this crop
      std::vector<larcv::Image2D> cropped_visi; // stores cropped visi image set for this crop
      std::vector<float> check_results(5,0);    // stores metrics for checking crop quality
      if ( _is_mc ) {
      
	LARCV_DEBUG() << "Start crop of Flow and Visibility images of image #" << icrop << std::endl;

	make_cropped_flow_images( src_plane, src_meta,
				  crop_v, *flo_v, *vis_v,
				  _thresholds_v,
				  _do_maxpool, _row_downsample_factor, _col_downsample_factor,
				  cropped_flow, cropped_visi,
				  &logger() );
      
	// check the quality of the crop
	if ( _check_flow ) {
	  check_results = check_cropped_images( src_plane, crop_v, _thresholds_v, cropped_flow, cropped_visi, _make_check_image, &(logger()), 0 );
	  UBCropLArFlow::_check_img_counter++;
	  
	  // check filter: has minimum visible pixels
	  if ( check_results[1]>=100 && check_results[2]>=100 ) {
	    passes_check_filter = true;
	  }
	  else {
	    passes_check_filter = false;
	  }
	}
      }//end of is mc
      
      if ( passes_check_filter || !_require_min_goodpixels ) {


	// if we are limiting overlaps, we need to mark overlap image
	if ( _limit_overlap ) {
	  const larcv::ImageMeta& cropped_src_meta = crop_v[src_plane]->meta();
	  int crop_r_start = src_meta.row( cropped_src_meta.min_y() );
	  int crop_c_start = src_meta.col( cropped_src_meta.min_x() );
	  for (int r=crop_r_start; r<crop_r_start+(int)cropped_src_meta.rows(); r++) {
	    for (int c=crop_c_start; c<crop_c_start+(int)cropped_src_meta.cols(); c++) {
	      overlap_img[0].set_pixel(r,c,1.0);
	    }
	  }
	}
	
	//LARCV_DEBUG() << "Store LArFlow Crop" << std::endl;
	if ( _save_output ) {
	  for ( auto& pimg : mod_crop_v )
	    ev_out_adc->emplace( std::move(*pimg) );
	  if ( _is_mc ) {
	    ev_vis_adc->emplace( std::move(cropped_visi) );
	    ev_flo_adc->emplace( std::move(cropped_flow) );
	  }
	}
	else {
	  for ( auto& pimg : mod_crop_v )
	    ev_out_adc->emplace( std::move(*pimg) );

	  if ( _is_mc ) {
	    for ( auto& img : cropped_visi )
	      ev_vis_adc->emplace( std::move(img) );
	    for ( auto& img : cropped_flow )
	      ev_flo_adc->emplace( std::move(img) );
	  }
	}
	//std::cout << "Store LArFlow Crops (nstored=" << ev_out_adc->image2d_array().size() << ")" << std::endl;
	
	// save meta
	ev_meta->store("nabove",int(check_results[0]));
	std::vector<int> nvis_v(2);
	nvis_v[0] = int(check_results[1]);
	nvis_v[1] = int(check_results[2]);
	ev_meta->store("nvis",nvis_v);
	std::vector<double> ncorrect_v(2);
	ncorrect_v[0] = check_results[3];
	ncorrect_v[1] = check_results[4];
	ev_meta->store("ncorrect",ncorrect_v);
	
	if ( _save_output ) {
	  foutIO->set_id( run, subrun, 100*event+icrop );
	  foutIO->save_entry();
	}
	nsaved++;
      }//end of if passes_check_fiter

      if ( _max_images>0 && nsaved>=_max_images )
	break;
    }
    
    return true;
  }
  
  void UBCropLArFlow::make_cropped_flow_images( const int src_plane,
						const larcv::ImageMeta& srcmeta,
						const std::vector<const larcv::Image2D*>& croppedadc_v,
						const std::vector<larcv::Image2D>& srcflow,
						const std::vector<larcv::Image2D>& srcvisi,
						const std::vector<float>& thresholds,
						const bool do_maxpool, const int row_ds_factor, const int col_ds_factor,
						std::vector<larcv::Image2D>& cropped_flow,
						std::vector<larcv::Image2D>& cropped_visi,
						const larcv::logger* log ) {

    // make the cropped flow images
    // we do not use any member variables of UBCropFlow instance in order to
    //   let this function be static, so other code can use it
    //
    // inputs
    // ------
    // src_plane: source image plane ID
    // srcmeta: meta of wholeview source image
    // croppedadc_v: vector of cropped adc image, provides meta for target images
    // srcflow: uncropped flow images
    // srcvisi: uncropped visibility images
    // thresholds: threshold below which pixel not considered matchable
    // maxpool: if true, we do max-pooling of cropped image
    // row_ds_factor: if maxpool=true, the downsampling factor for the row dim
    // col_ds_factor: if maxpool=true, the downsampling factor for the col dim    
    //
    // output
    // ------
    // croppedadc_v: if maxpool, then will modify croppedadc_v to maxpool
    // cropped_flow: cropped lar flow images. one for each target plane (tot 2)
    // cropped_visi: cropped visibility images. one for each target plane (tot 2)
    //
    // misc
    // ----
    // log: pointer to logger of calling processor

    // target planes given source plane
    const int targetplanes[3][2] = { {1,2},
				     {0,2},
				     {0,1} };
    // index of source flow and visibility image in srcflow and srcvisi, respectively
    const int targetindex[3][2] = { {0,1},
				    {2,3},
				    {4,5} };

    // get source adc and target image/meta
    const larcv::Image2D& adcimg = *(croppedadc_v[src_plane]);
    const larcv::ImageMeta& meta = adcimg.meta();

    const larcv::ImageMeta* target_meta[2];
    target_meta[0] = &(croppedadc_v[ targetplanes[src_plane][0] ]->meta());
    target_meta[1] = &(croppedadc_v[ targetplanes[src_plane][1] ]->meta());

    // make output flow/vis images
    // make two per source image
    cropped_flow.clear();
    cropped_visi.clear();

    // the meta for these are the same as the source image
    for (int i=0; i<2; i++) {
      larcv::Image2D flo( meta );
      flo.paint(UBCropLArFlow::_NO_FLOW_VALUE_ );
      cropped_flow.emplace_back( std::move(flo) );
      
      larcv::Image2D visi( meta );
      visi.paint(0.0);
      cropped_visi.emplace_back( std::move(visi) );
    }
        
    // loop over the two targets
    for (int i=0; i<2; i++) {

      int trgt_idx = targetindex[src_plane][i];
      int trgt_pl  = targetplanes[src_plane][i];
      const larcv::Image2D& target_adc = *croppedadc_v[trgt_pl];
      const larcv::Image2D& source_vis = srcvisi[trgt_idx];
      const larcv::Image2D& source_flo = srcflow[trgt_idx];
      larcv::Image2D& out_vis = cropped_visi[i];
      larcv::Image2D& out_flo = cropped_flow[i];

      // first, copy over values from the whole-view flo and vis images
      out_vis.copy_region( source_vis );
      out_flo.copy_region( source_flo );

      // VECTORIZED VERSION
      if ( _fusevector ) {
	std::vector<float> colindex(out_flo.meta().cols());
	for ( int i=0; i<(int)colindex.size(); i++ )
	  colindex[i] = i;
	
	// now, we must update, because
	//  (1) some pixels might not be visible
	//  (2) we need to remove the offset
	
	// determine flow offset
	// the difference in pixels of the x-origin of the source and target image
	float _flowoffset = (float)srcmeta.col(out_flo.meta().min_x()) - (float)srcmeta.col(target_meta[i]->min_x());
	float target_xmin = target_meta[i]->min_x();
	float target_xmax = target_meta[i]->max_x();
	//std::cout << "Flow offset=" << _flowoffset << std::endl;
	
	// mask below threshold
	const std::vector<float>& adcimg_asvec = adcimg.as_vector();
	std::transform( adcimg_asvec.begin(), adcimg_asvec.end(),
			out_flo.as_vector().begin(),
			out_flo.as_mod_vector().begin(),
			MaskBelowThreshold( thresholds[src_plane], UBCropLArFlow::_NO_FLOW_VALUE_ ) );
	
	// rows are contiguous
	std::vector<float> col_axis  = out_vis.meta().xaxis();      
	for ( size_t c=0; c<out_vis.meta().cols(); c++) {
	  std::transform( out_flo.row_start(c), out_flo.row_end(c),
			  out_vis.row_start(c), 
			  out_vis.row_start(c),
			  ModVisibility(target_xmin,target_xmax, col_axis[c] ) );
	}
	
	// we can update the flow in one shot
	std::vector<float>& flo_img_v = out_flo.as_mod_vector();
	std::transform( flo_img_v.begin(), flo_img_v.end(),
			flo_img_v.begin(),
			FlowOffset(_flowoffset) );
	
	// finally mask the flow value where visibility is now-zero
	std::transform( out_vis.as_vector().begin(), out_vis.as_vector().end(),
			out_flo.as_vector().begin(),
			out_flo.as_mod_vector().begin(),
			MaskBelowThreshold( 1.0, UBCropLArFlow::_NO_FLOW_VALUE_ ) );
	
	// need to mask when target pixel adc is below threshold. how to vectorize?
	// first get adc values from target image by following flow
	std::vector<float> flowadc( out_flo.as_vector().size(), UBCropLArFlow::_NO_FLOW_VALUE_ );
	for ( int c=0; c<(int)out_flo.meta().cols(); c++ ) {
	  std::vector<float> target_adc_v( target_adc.row_start_const(c), target_adc.row_end_const(c) );
	  std::transform( out_flo.row_start(c), out_flo.row_end(c),
			  flowadc.begin()+c*out_flo.meta().rows(),
			  FollowFlow( target_adc_v, c ) );
	}
	// next, mask flow values when targetadc is below threshold
	std::transform( flowadc.begin(), flowadc.end(),
			out_vis.as_vector().begin(),
			out_vis.as_mod_vector().begin(),
			MaskBelowThreshold( thresholds[trgt_pl], 0.0 ) );

      }
      else {
	// Use scalar code
      
	// we try to avoid this loop
	// loop over rows and cols
	for (int r=0; r<(int)meta.rows(); r++) {
	  
	  float tick  = meta.pos_y(r);
	  int src_row = srcmeta.row(tick);
	  
	  for (int c=0; c<(int)meta.cols(); c++) {
	    
	    float wire  = meta.pos_x(c);
	    int src_col = srcmeta.col(wire);
	    
	    // flow and visi values (in the source coordinate system)
	    float visi = source_vis.pixel( src_row, src_col );
	    float flow = source_flo.pixel( src_row, src_col );
	    
	    // is the target pixel in the cropped image?
	    float target_wire = (wire + flow);
	    if ( target_wire<target_meta[i]->min_x() || target_wire>=target_meta[i]->max_x() ) {
	      // outside the target cropped image
	      out_flo.set_pixel( r, c, UBCropLArFlow::_NO_FLOW_VALUE_ );
	      out_vis.set_pixel( r, c, 0.0 );
	      // this shouldn't happen unless a mistake?
	    }
	    else {
	      
	      // inside. need to adjust the flow.
	      int target_col = target_meta[i]->col( target_wire );
	      float target_adc = croppedadc_v[ targetplanes[src_plane][i] ]->pixel( r, target_col );
	      
	      float target_flow = target_col - c;
	      out_flo.set_pixel( r, c, target_flow );
	      out_vis.set_pixel( r, c, visi );
	      
	      if ( target_adc<thresholds[ targetplanes[src_plane][i] ] )
		out_vis.set_pixel( r, c, 0.0 );
	      
	    }
	    
	  }//end of loop over cols
	}//end of loop over rows
      }//end of if vectorize/scalar code
      
      // // max pool, if desired
      // if ( do_maxpool && (row_ds_factor>1 || col_ds_factor>1) ) {
      // 	// allocate empty images
      // 	larcv::Image2D ds_src_adc;
      // 	larcv::Image2D ds_target_adc;
      // 	larcv::Image2D ds_flow;
      // 	larcv::Image2D ds_visi;
      // 	maxPool( row_ds_factor, col_ds_factor,
      // 		 adcimg, target_adc, out_flo, out_vis,
      // 		 thresholds,
      // 		 ds_src_adc, ds_target_adc, ds_flow, ds_visi );
      // 	// store source maxpool
      // 	maxpooled_source.emplace_back(std::move(ds_src_adc));
      // 	// store target_maxpool
      // 	maxpooled_target.emplace_back(std::move(ds_target_adc));
	
      // 	// swap the flow and visi
      // 	std::swap( out_flo, ds_flow );
      // 	std::swap( out_vis, ds_visi );
      // }
      
    }// end of loop over target images
      

    // if (do_maxpool) {
    //   // swap out downsampled adc images
    //   std::swap( croppedadc_v[src_plane], maxpooled_source[0] );
    //   for (int i=0; i<2; i++) {
    // 	int trgt_pl  = targetplanes[src_plane][i];
    // 	std::swap( croppedadc_v[trgt_pl], maxpooled_target[i] );
    //   }
    // }
    
    return;
  }

  void UBCropLArFlow::downsample_crops( const std::vector<larcv::Image2D*>& cropped_adc_v,
					const std::vector<larcv::Image2D>& cropped_flow_v,
					const std::vector<larcv::Image2D>& cropped_visi_v,
					std::vector<larcv::Image2D>& downsampled_adc_v,
					std::vector<larcv::Image2D>& downsampled_flow_v,
					std::vector<larcv::Image2D>& downsampled_visi_v ) {
    return; // to do
  }

  std::vector<float> UBCropLArFlow::check_cropped_images( const int src_plane,
							  const std::vector<const larcv::Image2D*>& cropped_adc_v,
							  const std::vector<float>& thresholds,
							  const std::vector<larcv::Image2D>& cropped_flow,
							  const std::vector<larcv::Image2D>& cropped_visi,
							  const bool visualize_flow,
							  const larcv::logger* log, const int verbosity ) {

    // we follow the flow to the target image.
    // correct if
    //  source adc above threshold has flow value
    //  and ( source visi=1+target image pixel above threshold OR target image pixel below and source visi is 0 )
    
    // inputs
    // ------
    // src_plane: index of source plane
    // cropped_adc_v: ADC images for all 3 planes
    // thresholds: ADC thresholds for each plane
    // cropped_flow: cropped flow image (2 per source plane)
    // cropped_visi: cropped visi image (2 per source plane)
    // visualize_flow: will dump out pngs using ROOT
    // log: pointer to a larbys::logger object. (default=NULL)
    // verbosity: unused right now
    //
    // outputs
    // -------
    // vector<float>: result of checks. entries:
    //   [0]: number of source values above threshold
    //   [1,2]: number of visible pixels (vis=1) in source
    //   [3,4]: num correct

    
    const int targetplanes[3][2] = { {1,2},
				     {0,2},
				     {0,1} };
    
    int nabove[2] = {0};              // pixels above threshold
    int nvis[2] = {0};                // pixels with vis=1
    int nwrong_flow2nothing[2] = {0}; // flow goes to target pixel with adc below thresh
    int nwrong_flowob[2] = {0};       // flow goes out of target image bounds
    int nwrong_visbelow[2] = {0};     // vis=1 for below thresh pixel
    int nwrong_badvisi[2] = {0};      // vis=0 but flow goes to above thresh target pixel
    int nwrong_nolabel[2] = {0};      // vis=1 and above thresh, but no flow value
    int ncorrect[2] = {0};            // total correct pixels
    
    const larcv::Image2D& src_adc = *(cropped_adc_v[src_plane]);
    const larcv::ImageMeta& meta = src_adc.meta();


    TH2D* hflow[2] = {NULL};
    TH2D* hvisi[2] = {NULL};
    TH2D* hcheck_flow[2] = {NULL};     // visualize images
    TH2D* hcheck_vismatch[2] = {NULL}; // visualize images
    if ( visualize_flow ) {
      
      const larcv::ImageMeta* tar_meta[2];
      for (int i=0; i<2; i++) {
	tar_meta[i] = &(cropped_adc_v[ targetplanes[src_plane][i] ]->meta());
      }
	
      // plot the flow (in source coordinates)
      std::stringstream ss_y2uflow;
      ss_y2uflow << "hflow_" << src_plane << "to" << targetplanes[src_plane][0] << "_" << _check_img_counter;      
      hflow[0] = new TH2D( ss_y2uflow.str().c_str(), "", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );

      std::stringstream ss_y2vflow;
      ss_y2vflow << "hflow_" << src_plane << "to" << targetplanes[src_plane][1] << "_" << _check_img_counter;      
      hflow[1] = new TH2D( ss_y2vflow.str().c_str(), "", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );

      // plot the visibility (in source coordinates)
      std::stringstream ss_y2uvis;;
      ss_y2uvis << "hvisi_" << src_plane << "to" << targetplanes[src_plane][0] << "_" << _check_img_counter;      
      hvisi[0] = new TH2D( ss_y2uvis.str().c_str(), "", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );
      
      std::stringstream ss_y2vvis;
      ss_y2vvis << "hvisi_" << src_plane << "to" << targetplanes[src_plane][1] << "_" << _check_img_counter;      
      hvisi[1] = new TH2D( ss_y2vvis.str().c_str(), "", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );
      
      // we follow the flow and mark values (in the target coordinate system)
      std::stringstream ss1;
      ss1 << "hcheck_" << src_plane << "to" << targetplanes[src_plane][0] << "_" << _check_img_counter;
      std::stringstream ss1_t;
      ss1_t << "Flow/Vis Check: plane" << src_plane << " to " << targetplanes[src_plane][0] << " #" << _check_img_counter;
      hcheck_flow[0] = new TH2D( ss1.str().c_str(), ss1_t.str().c_str(),
				 tar_meta[0]->cols(), tar_meta[0]->min_x(), tar_meta[0]->max_x(),
				 tar_meta[0]->rows(), tar_meta[0]->min_y(), tar_meta[0]->max_y() );

      std::stringstream ss2;
      ss2 << "hcheck_" << src_plane << "to" << targetplanes[src_plane][1] << "_" << _check_img_counter;
      std::stringstream ss2_t;
      ss2_t << "Flow/Vis Check: plane" << src_plane << " to " << targetplanes[src_plane][1] << " #" << _check_img_counter;
      hcheck_flow[1] = new TH2D( ss2.str().c_str(), ss2_t.str().c_str(),
				 tar_meta[1]->cols(), tar_meta[1]->min_x(), tar_meta[1]->max_x(),
				 tar_meta[1]->rows(), tar_meta[1]->min_y(), tar_meta[1]->max_y() );

      // we indicate if match found in target image (in source coordinates)
      std::stringstream ss3;
      ss3 << "hvismatch_" << src_plane << "to" << targetplanes[src_plane][0] << "_" << _check_img_counter;
      hcheck_vismatch[0] = new TH2D( ss3.str().c_str(), ss3.str().c_str(), meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );

      std::stringstream ss4;
      ss4 << "hvismatch_" << src_plane << "to" << targetplanes[src_plane][1] << "_" << _check_img_counter;
      hcheck_vismatch[1] = new TH2D( ss4.str().c_str(), ss4.str().c_str(), meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );
    }
    
	
    for (int i=0; i<2; i++) {
      int trgt_img = targetplanes[src_plane][i];
      
      
      for (int r=0; r<(int)meta.rows(); r++) {
	for (int c=0; c<(int)meta.cols(); c++) {
	  // loop over target plane

	  float visi = cropped_visi[i].pixel(r,c);
	  float flow = cropped_flow[i].pixel(r,c);
	  if ( visualize_flow ) {
	    hvisi[i]->SetBinContent(c+1,r+1,visi);
	    hflow[i]->SetBinContent(c+1,r+1,flow);
	  } 
	  
	  if ( visi>0.5 )
	    nvis[i]++;
	  
	  if ( src_adc.pixel(r,c)<thresholds[src_plane] ) {
	    if ( visi>0.5 ) {
	      nwrong_visbelow[i]++;
	      if ( visualize_flow ) {
		hcheck_vismatch[i]->SetBinContent(c+1,r+1,-5.0);
	      }
	    }
	    continue; // continue, as source below threshold now, not interesting
	  }

	  // above threshold pixels
	  nabove[i]++;
	  	  
	  int targetc = c+flow;
	  if ( flow<=UBCropLArFlow::_NO_FLOW_VALUE_ ) {
	    if ( visi<0.5) {
	      ncorrect[i]++;
	      if ( visualize_flow )
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, 1.0 );
	    }
	    else {
	      nwrong_nolabel[i]++;
	      if ( visualize_flow )	      
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, -1.0 );
	    }
	    continue; // no valid value of flow, continue
	  }	  
	
	  float targetadc = 0;
	  bool  goodtarget = false;
	  try {
	    targetadc = cropped_adc_v[ trgt_img ]->pixel( r, targetc );
	    goodtarget = true;
	  }
	  catch (std::exception& e) {
	    goodtarget = false;
	    std::cout << __PRETTY_FUNCTION__ << ":" << __FILE__ << "." << __LINE__ << ": good vis pixel flow out of bounds. "
		      << "flow: " << c << "->" << targetc 
		      << "\n\t"
		      << e.what()
		      << std::endl;
	  }

	  
	  if (!goodtarget) {
	    nwrong_flowob[i]++;
	    if ( visualize_flow ) 
	      hcheck_vismatch[i]->SetBinContent( c+1, r+1, -4.0 );
	    continue; // flow out of bounds
	  }
	  	  
	  if ( visi>0.5 ) {
	    // should be vis
	    if ( targetadc>=thresholds[ trgt_img ] ) {
	      ncorrect[i]++;
	      if ( visualize_flow ) {
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, 3.0 );
		hcheck_flow[i]->SetBinContent( targetc+1, r+1, 3.0 );
	      }
	    }
	    else {
	      nwrong_flow2nothing[i]++;
	      if ( visualize_flow ) {
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, -3.0 );
		hcheck_flow[i]->SetBinContent( targetc+1, r+1, -3.0 );		
	      }
	    }
	    
	  }
	  else {
	    // should be invisible
	    if ( targetadc<thresholds[ trgt_img ] ) {
	      ncorrect[i]++;
	      if ( visualize_flow ) {
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, 2.0 );
		hcheck_flow[i]->SetBinContent( targetc+1, r+1, 2.0 );
	      }
	    }
	    else {
	      nwrong_badvisi[i]++;
	      if ( visualize_flow ) {
		hcheck_vismatch[i]->SetBinContent( c+1, r+1, -2.0 );
		hcheck_flow[i]->SetBinContent( targetc+1, r+1, -2.0 );
	      }
	    }
	  }
	}
      
      }
    }

    if ( log!=NULL ) {
      for (int i=0; i<2; i++) {
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "      	
	  << "[source plane " << src_plane << "-> target plane " << targetplanes[src_plane][i]  << "]" << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)	  
	  << "  nabovethresh=" << nabove[i] << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  << "  nvis=" << nvis[i] << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  << "  ncorrect=" << ncorrect[i] << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  << "  ncorrect/nabove=" << float(ncorrect[i])/float(nabove[i]) << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "
	  << "  badvisi/nabove: "      << float(nwrong_badvisi[i])/float(nabove[i]) << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "
	  << "  visbelow/nabove: "      << float(nwrong_visbelow[i])/float(nvis[i]) << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "      
	  << "  flow2nothing/nabove: " << float(nwrong_flow2nothing[i])/float(nabove[i]) << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)      
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "      
	  << "  flow2OB/nabove: " << float(nwrong_flowob[i])/float(nabove[i]) << std::endl;
	(*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	  //std::cout << __FILE__ << "." << __LINE__ << "::" << __FUNCTION__ << ": "      
	  << "  nolabel/nabove: " << float(nwrong_nolabel[i])/float(nabove[i]) << std::endl;
      }
    }
      
    // dump canvas and clean up
    if ( visualize_flow ) {
      for (int i=0; i<2; i++) {
	TCanvas c("c", "", 1600, 1800 );
	c.Divide(2,3);
	c.cd(1);

	// SOURCE ADC
	std::stringstream ss1;
	ss1 << "hsource" << i << "_p" << src_plane << "_" << _check_img_counter;
	TH2D hsrc = as_th2d( src_adc, ss1.str() );
	hsrc.SetTitle("Source ADC");
	hsrc.SetMaximum(50.0);
	hsrc.SetMinimum( 0.0);
	setRainbowPalette();
	hsrc.Draw("COLZ");
      

	// TARGET ADC
	c.cd(2);
	std::stringstream ss2;
	ss2 << "htarget" << i << "_" << src_plane << "to" << targetplanes[src_plane][i] << "_" << _check_img_counter;
	TH2D htar1 = as_th2d( *cropped_adc_v[targetplanes[src_plane][i]], ss2.str() );
	htar1.SetTitle("TargetADC");
	htar1.SetMaximum(50.0);
	htar1.SetMinimum( 0.0);
	setRainbowPalette();      
	htar1.Draw("COLZ");

	// FLOW (drawn in source coordinates)
	c.cd(3);
	hflow[i]->SetTitle("Flow");
	hflow[i]->SetMaximum( 900);
	hflow[i]->SetMinimum(-900);
	setBWRPalette();
	hflow[i]->Draw("COLZ");

	// VISIBILITY (drawn in source coordinates)
	c.cd(6);
	hvisi[i]->SetTitle("Y2U Visibility");
	hvisi[i]->Draw("COLZ");
      
	c.cd(5);
	// check vis (drawn in source coordinates)
	hcheck_vismatch[i]->SetTitle("Check Visibility Match");      
	hcheck_vismatch[i]->SetMaximum( 5.0);
	hcheck_vismatch[i]->SetMinimum(-5.0);
	setBWRPalette();      
	hcheck_vismatch[i]->Draw("COLZ");

	c.cd(4);
	// check flow (drawn in target coordinates)
	hcheck_flow[i]->SetMaximum(5.0);
	hcheck_flow[i]->SetMinimum(-5.0);
	setBWRPalette();
	hcheck_flow[i]->Draw("COLZ");

	// save
	std::stringstream css1;
	css1 << "ccheck_" << src_plane << "to" << targetplanes[src_plane][i] << "_" << _check_img_counter << ".png";
	c.SaveAs( css1.str().c_str() );
      }
      
      for (int i=0; i<2; i++) {
	delete hflow[i];
	delete hvisi[i];
	delete hcheck_flow[i];
	delete hcheck_vismatch[i];
      }
    }

    std::vector<float> results(5);
    results[0] = nabove[0];
    results[1] = nvis[0];
    results[2] = nvis[1];    
    results[3] = ncorrect[0];
    results[4] = ncorrect[1];

    return results;
  }
  
  void UBCropLArFlow::finalize()
  {
    if ( _save_output )
      foutIO->finalize();
  }

  void UBCropLArFlow::maxPool( const int row_downsample_factor, const int col_downsample_factor,
			       const larcv::Image2D& src_adc, const larcv::Image2D& target_adc,
			       const larcv::Image2D& flow, const larcv::Image2D& visi,
			       const std::vector<float>& thresholds,
			       larcv::Image2D& ds_src_adc, larcv::Image2D& ds_target_adc,
			       larcv::Image2D& ds_flow, larcv::Image2D& ds_visi ) {
    
    // to make things easy, we require that the downsampling factors divide evenly into the rows and cols
    if ( src_adc.meta().rows()%row_downsample_factor!=0 ) {
      std::stringstream ss;
      ss << "Rows of image (" << src_adc.meta().rows() << ") are not a multiple of the row downsample factor " << row_downsample_factor << std::endl;
      throw std::runtime_error( ss.str() );
    }
    if ( src_adc.meta().cols()%col_downsample_factor!=0 ) {
      std::stringstream ss;
      ss << "Cols of image (" << src_adc.meta().cols() << ") are not a multiple of the col downsample factor " << col_downsample_factor << std::endl;
      throw std::runtime_error( ss.str() );
    }


    // make the output images

    // source ADC
    const larcv::ImageMeta& meta_src = src_adc.meta();
    larcv::ImageMeta out_meta_src( meta_src.min_x(), meta_src.min_y(), meta_src.max_x(), meta_src.max_y(),
				   meta_src.rows()/row_downsample_factor, meta_src.cols()/col_downsample_factor,
				   meta_src.id(), meta_src.unit() );
    larcv::Image2D out_src( out_meta_src );

    // target ADC
    const larcv::ImageMeta& meta_target = target_adc.meta();
    larcv::ImageMeta out_meta_target( meta_target.min_x(), meta_target.min_y(), meta_target.max_x(), meta_target.max_y(),
				   meta_target.rows()/row_downsample_factor, meta_target.cols()/col_downsample_factor,
				   meta_target.id(), meta_target.unit() );
    larcv::Image2D out_target( out_meta_target );
    out_target.paint(-1.0);
    
    // flow ADC
    const larcv::ImageMeta& meta_flow = flow.meta();
    larcv::ImageMeta out_meta_flow( meta_flow.min_x(), meta_flow.min_y(), meta_flow.max_x(), meta_flow.max_y(),
				   meta_flow.rows()/row_downsample_factor, meta_flow.cols()/col_downsample_factor,
				   meta_flow.id(), meta_flow.unit() );
    larcv::Image2D out_flow( out_meta_flow );

    // visi ADC
    const larcv::ImageMeta& meta_visi = visi.meta();
    larcv::ImageMeta out_meta_visi( meta_visi.min_x(), meta_visi.min_y(), meta_visi.max_x(), meta_visi.max_y(),
				   meta_visi.rows()/row_downsample_factor, meta_visi.cols()/col_downsample_factor,
				   meta_visi.id(), meta_visi.unit() );
    larcv::Image2D out_visi( out_meta_visi );
    out_visi.paint(0.0);

    // we loop through output rows and cols
    // we ask which pixel in the neighborhood of the src and target image is largest
    // for the source image, that pixel position is used to determine what value of the flow
    //  and visibility are chosen.
    // the flow is corrected to be on the scale of the downsampled images
    const larcv::ImageMeta& out_src_meta = out_src.meta();
    const larcv::ImageMeta& out_tar_meta = out_target.meta();        
    for ( int ro=0; ro<(int)out_src_meta.rows(); ro++) {
      for (int co=0; co<(int)out_src_meta.cols(); co++) {

	int ri_start = ro*row_downsample_factor;
	int ri_end   = (ro+1)*row_downsample_factor;
	int ci_start = co*col_downsample_factor;
	int ci_end   = (co+1)*col_downsample_factor;

	// for the block of input pixels, we determine maximum ADC deposited (above threshold)
	bool oneabove = false;
	float maxadc = thresholds[ (int)(meta_src.id()) ];
	int max_row = -1;
	int max_col = -1;
	float max_flow = UBCropLArFlow::_NO_FLOW_VALUE_;
	float max_visi = UBCropLArFlow::_NO_FLOW_VALUE_;
	for (int ri=ri_start; ri<ri_end; ri++) {
	  for (int ci=ci_start; ci<ci_end; ci++) {

	    try {
	    
	      if ( src_adc.pixel( ri, ci )>=maxadc
		   && flow.pixel(ri,ci)>UBCropLArFlow::_NO_FLOW_VALUE_ ) {
		maxadc = src_adc.pixel(ri,ci);
		oneabove = true;
		max_flow = flow.pixel(ri,ci);
		max_visi = visi.pixel(ri,ci);
		max_row = ri;
		max_col = ci;
	      }
	    }
	    catch ( std::exception& e ) {
	      std::stringstream ss;
	      ss << __PRETTY_FUNCTION__ << "@" << __FILE__ << "." << __LINE__ << ": error retrieving source info\n\t"
		 << e.what() << std::endl;
	      throw std::runtime_error( ss.str() );
	    }
	  }
	}
	
	if ( !oneabove ) {
	  try {
	    out_src.set_pixel( ro, co, 0.0 );
	    out_flow.set_pixel( ro, co, UBCropLArFlow::_NO_FLOW_VALUE_ );
	    out_visi.set_pixel( ro, co, 0.0 );
	  }
	  catch ( std::exception& e ) {
	    std::stringstream ss;
	    ss << __PRETTY_FUNCTION__ << "@" << __FILE__ << "." << __LINE__ << ": error setting empty output pixel\n\t"
	       << out_src.meta().dump() << "\n\t"
	       << out_flow.meta().dump() << "\n\t"
	       << out_visi.meta().dump() << "\n\t"    
	       << e.what() << std::endl;
	    throw std::runtime_error( ss.str() );
	  }
	  continue;
	}

	// now we have to adjust the flow to get to the correct target pixel
	int orig_out_col = max_col + max_flow;
	int outcol       = orig_out_col/col_downsample_factor;
	int oflow        = outcol-co;
	if ( outcol>=(int)out_tar_meta.cols() || outcol<0 ) {
	  std::stringstream ss;
	  ss << __PRETTY_FUNCTION__ << "." << __LINE__ << ": adjusted flow " << co << "->" << outcol << " does not fall within max-pooled image. "
	     << " original flow=" << max_flow << " dsflow=" << oflow << "\n"
	     << out_tar_meta.dump()
	     << std::endl;
	  throw std::runtime_error( ss.str() );
	}

	// valid flow
	// whats the max adc value in target region
	bool targetabove = false;
	float target_maxadc = 0;
	for ( int tci=outcol*col_downsample_factor; tci<(outcol+1)*col_downsample_factor; tci++) {
	  for (int tri=ro*row_downsample_factor; tri<(ro+1)*row_downsample_factor; tri++) {
	    float target_flow_adc = target_adc.pixel( tri, tci );
	    if ( target_flow_adc>thresholds[out_tar_meta.id()] ) {
	      targetabove = true;
	      target_maxadc = target_flow_adc;
	    }
	  }
	}
	
	// checks out. fill the outputs
	try {
	  out_src.set_pixel( ro, co, src_adc.pixel( max_row, max_col ) );
	  out_flow.set_pixel( ro, co, oflow );	  
	  if ( targetabove ) {
	    out_visi.set_pixel( ro, co, 1.0 );
	    out_target.set_pixel( ro, outcol, target_maxadc );
	  }
	  else {
	    out_visi.set_pixel( ro, co, 0.0 );
	    out_target.set_pixel( ro, outcol, 0.0 );
	  }
	}
	catch ( std::exception& e ) {
	  std::stringstream ss;
	  ss << __PRETTY_FUNCTION__ << "@" << __FILE__ << "." << __LINE__ << ": error setting pooled output pixel.  "
	     << "pooled flow=" << co << "->" << oflow << " max_flow=" << max_col << "->" << max_col+max_flow << "\n\t"
	     << e.what() << std::endl;
	  throw std::runtime_error( ss.str() );
	}
      }//end of loop over output columns
    }//end of loop over output rows

    // the target image will be incomplete for pixels where source did not project into it
    // so we fill the remainder here
    for ( int ro=0; ro<(int)out_tar_meta.rows(); ro++) {
      for (int co=0; co<(int)out_tar_meta.cols(); co++) {
	if ( out_target.pixel(ro,co)>=0.0 )
	  continue; // because already filled
	
	int ri_start = ro*row_downsample_factor;
	int ri_end   = (ro+1)*row_downsample_factor;
	int ci_start = co*col_downsample_factor;
	int ci_end   = (co+1)*col_downsample_factor;

	// for the block of input pixels, we determine maximum ADC deposited (above threshold)
	bool oneabove = false;
	float maxadc = thresholds[ (int)(meta_target.id()) ];
	for (int ri=ri_start; ri<ri_end; ri++) {
	  for (int ci=ci_start; ci<ci_end; ci++) {
	    if ( target_adc.pixel( ri, ci )>maxadc ) {
	      oneabove = true;
	      maxadc = target_adc.pixel(ri,ci);
	    }
	  }
	}
	if ( oneabove ) {
	  try {
	    out_target.set_pixel( ro, co, maxadc );
	  }
	  catch (std::exception& e) {
	    std::stringstream ss;
	    ss << __PRETTY_FUNCTION__ << "@" << __FILE__ << "." << __LINE__ << ": error max-pooling output. \n\t"
	       << e.what() << std::endl;
	    throw std::runtime_error(ss.str());
	  }
	}
      }//end of output col loop
    }//end of output roow loop
    
    // output by a swap to avoid a copy
    // object passed to us gets destroyed when function ends
    std::swap( ds_src_adc,    out_src );
    std::swap( ds_target_adc, out_target );
    std::swap( ds_flow,       out_flow );
    std::swap( ds_visi,       out_visi );
    
  }
  
  void UBCropLArFlow::setBWRPalette() {
    
    // A colour palette that goes blue->white->red, useful for
    // correlation matrices
    const int NRGBs = 3;
    const int n_color_contours = 999;
    
    if ( UBCropLArFlow::_colors==NULL ) {
      _colors=new int[n_color_contours];
      
      Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
      Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
      Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
      Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
      int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
      for(uint i=0; i<n_color_contours; ++i) UBCropLArFlow::_colors[i]=colmin+i;
    }
    
    gStyle->SetNumberContours(n_color_contours);
    gStyle->SetPalette(n_color_contours, UBCropLArFlow::_colors);
  }
  
  void UBCropLArFlow::setRainbowPalette() {
    gStyle->SetPalette( kRainBow );
  }


}
#endif
