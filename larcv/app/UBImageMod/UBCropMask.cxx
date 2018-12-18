#ifndef __UBCROPMASK_CXX__
#define __UBCROPMASK_CXX__

#include <sstream>

#include "UBCropMask.h"
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

  static UBCropMaskProcessFactory __global_UBCropMaskProcessFactory__;
  int   UBCropMask::_check_img_counter = 0;
  const float UBCropMask::_NO_FLOW_VALUE_ = -4000;

  UBCropMask::UBCropMask(const std::string name)
    : ProcessBase(name)
  {}

  void UBCropMask::configure(const PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");
    _input_bbox_producer    = cfg.get<std::string>("InputBBoxProducer");
    _input_adc_producer     = cfg.get<std::string>("InputADCProducer");
    _input_cropped_producer = cfg.get<std::string>("InputCroppedADCProducer");
    _input_masks_producer  = cfg.get<std::string>("InputMasksProducer");

    //
    //
    _output_adc_producer    = cfg.get<std::string>("OutputCroppedADCProducer");
    //I added these -J
    _output_masks_producer = cfg.get<std::string>("OutputMasksProducer");
    _output_weights_producer= cfg.get<std::string>("OutputWeightsProducer");
    //

    _output_meta_producer   = cfg.get<std::string>("OutputCroppedMetaProducer");
    _output_filename        = cfg.get<std::string>("OutputFilename");

    _max_images             = cfg.get<int>("MaxImages",-1);
    _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) );








    // sparsity requirement: prevent cropping over the same regions
    _limit_overlap        = cfg.get<bool>("LimitOverlap",false);
    _max_overlap_fraction = cfg.get<float>("MaxOverlapFraction", 0.5 );

    // debug options
    _check_flow             = cfg.get<bool>("CheckFlow",false);      // output to screen, checks of cropped images
    _make_check_image       = cfg.get<bool>("MakeCheckImage",false); // dump png of image checks

    //Here we get the offsets for the weight image classes


    if ( _make_check_image )
      gStyle->SetOptStat(0);

    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();

  }

  void UBCropMask::initialize()
  {}

  bool UBCropMask::process(IOManager& mgr)
  {
    LARCV_NORMAL() << "\n\nBeginning Process!\n\n\n"   << std::endl;
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // reset counters
    nCosmic = 0;
    nElectron = 0;
    nProton = 0;
    nOther = 0;
    nNeutrino = 0;
    nNeutron = 0;



    // input ADC
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_v = ev_in_adc->image2d_array();
    const std::vector<larcv::ImageMeta> meta_orig_v = {img_v.at(0).meta(), img_v.at(1).meta(), img_v.at(2).meta() };

    // // input Labels
    // auto ev_in_labels  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_labels_producer));
    // if (!ev_in_labels) {
    //   LARCV_CRITICAL() << "No Input Labels Image2D found with a name: " << _input_labels_producer << std::endl;
    //   throw larbys();
    // }
    // const std::vector< larcv::Image2D >& labels_v = ev_in_labels->image2d_array();

    // input masks
    auto ev_in_masks  = (larcv::EventClusterMask*)(mgr.get_data("clustermask", _input_masks_producer));
    if (!ev_in_masks) {
      LARCV_CRITICAL() << "No Input ClusterMasks found with a name: " << _input_masks_producer << std::endl;
      throw larbys();
    }
    const std::vector< std::vector< larcv::ClusterMask > >& masks_vv = ev_in_masks->clustermask_array();



    // input BBox of crop
    auto  ev_in_bbox  = (larcv::EventBBox2D*)(mgr.get_data("bbox2d", _input_bbox_producer));
    if (!ev_in_bbox) {
      LARCV_CRITICAL() << "No Input BBox2D found with a name: " << _input_bbox_producer << std::endl;
      throw larbys();
    }


    // cropped input ADC
    auto ev_in_cropped  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_cropped_producer));
    if (!ev_in_cropped) {
      LARCV_CRITICAL() << "No Input Cropped ADC Image2D found with a name: " << _input_cropped_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& cropped_v = ev_in_cropped->image2d_array();


    // ----------------------------------------------------------------

    // Output  containers
    larcv::EventImage2D* ev_out_adc  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_adc_producer);
    larcv::EventClusterMask* ev_out_masks  = (larcv::EventClusterMask*)foutIO->get_data("clustermask",_output_masks_producer);
    // larcv::EventImage2D* ev_out_weights  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_weights_producer);
    ev_out_adc->clear();
    ev_out_masks->clear();
    // ev_out_weights->clear();

    // Output Meta containers
    larcv::EventMeta*    ev_meta     = (larcv::EventMeta*)foutIO->get_data("meta",_output_meta_producer);
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

    std::vector<larcv::Image2D> overlap_img; // store an image used to keep track of previously cropped pixels
    if ( _limit_overlap ) {
      larcv::Image2D overlap( img_v[src_plane].meta() );
      overlap.paint(0.0);
      overlap_img.emplace_back( std::move(overlap) );
    }

    for (int icrop=0; icrop<ncrops; icrop++) {

      // this is a copy. not great. could swap if needed ...
      std::vector<larcv::Image2D> crop_v;
      for (int i=0; i<3; i++) {
	        crop_v.push_back( cropped_v.at( 3*icrop+i ) );
      }

      // if limiting overlap, check that we are not overlapping too many pixels
      if ( _limit_overlap ) {
        	const larcv::ImageMeta& cropped_src_meta = crop_v[src_plane].meta();
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
        	  LARCV_NORMAL() << "Skipping overlapping image. Frac overlap=" << frac_overlap << "." << std::endl;
        	  continue;
        	}
        	std::cout << "Overlap fraction: " << frac_overlap << std::endl;
      }

      LARCV_DEBUG() << "Start crop of Flow and Visibility images of image #" << icrop << std::endl;
      std::vector<larcv::Image2D> cropped_labels;
      // std::vector<larcv::Image2D> cropped_weights ={};

      std::cout << "Number Saved Here: "<< nsaved<< std::endl;
      // std::vector<larcv::ImageMeta> meta_orig_v = {img_v.at(0).meta(), img_v.at(1).meta(), img_v.at(2).meta() }
      std::vector<larcv::ImageMeta> meta_crop_v = {crop_v.at(0).meta(), crop_v.at(1).meta(), crop_v.at(2).meta() };

      //Replace this chunk with making the crops
      AdjustMasksUsingBBox2D( *ev_in_bbox,
        meta_orig_v,
        meta_crop_v,
				masks_vv,
        icrop,
				*ev_out_masks
      ); //,&logger()
      // const std::vector<std::vector< larcv::ClusterMask >>& cropped_masks_vv = ev_out_masks->clustermask_array();





      // check the quality of the crop
      bool passes_check_filter = true;
       std::vector<float> check_results(5,0);
  //     if ( _check_flow ) {
	// check_results = check_cropped_images( src_plane, crop_v, _thresholds_v, cropped_flow, cropped_visi, _make_check_image, &(logger()), 0 );
	// UBCropMask::_check_img_counter++;
  //
	// // check filter: has minimum visible pixels
	// if ( check_results[1]>=100 && check_results[2]>=100 ) {
	//   passes_check_filter = true;
	// }
	// else {
	//   passes_check_filter = false;
	// }
  //     }

      if ( passes_check_filter ) {


      	// if we are limiting overlaps, we need to mark overlap image
      	if ( _limit_overlap ) {
      	  const larcv::ImageMeta& cropped_src_meta = crop_v[src_plane].meta();
      	  int crop_r_start = src_meta.row( cropped_src_meta.min_y() );
      	  int crop_c_start = src_meta.col( cropped_src_meta.min_x() );
      	  for (int r=crop_r_start; r<crop_r_start+(int)cropped_src_meta.rows(); r++) {
      	    for (int c=crop_c_start; c<crop_c_start+(int)cropped_src_meta.cols(); c++) {
      	      overlap_img[0].set_pixel(r,c,1.0);
      	    }
      	  }
      	}

        // //INSERT WEIGHTED IMAGE making
        // larcv::Image2D weight_img(cropped_labels_v[2].meta());
        //
        // for (int plane = 0 ; plane<3; plane++){
        //
        //
        //
        //
        //   //std::cout << "plane!" << std::endl;
        //   int nbackground =0;
        //   int ntrack =0;
        //   int nshower =0;
        //   int nend =0;
        //
        //   int rows = cropped_labels_v[2].meta().rows();
        //   int cols = cropped_labels_v[2].meta().cols();
        //
        //
        //   weight_img.paint(0.0);
        //   for (int row = 0 ; row < rows; row++) {
        //     for (int col = 0 ; col < cols; col++) {
        //       int pix =cropped_labels_v[plane].pixel(row,col);
        //       if (pix == 0) {nbackground++;}
        //       if (pix == 1) {ntrack++;}
        //       if (pix == 2) {nshower++;}
        //       if (pix == 3) {nend++;}
        //
        //
        //     }//End of cols loop (Counting Totals)
        //
        //   }//End of rows loop (Count totals)
        //   //std::cout << "Finish RCloops1!" << std::endl;
        //   for (int row = 0 ; row < rows; row++) {
        //     for (int col = 0 ; col < cols; col++) {
        //       int pix =cropped_labels_v[plane].pixel(row,col);
        //       if (pix == 0) {weight_img.set_pixel(row,col,1.0/(nbackground+_background_offset));}
        //       else if (pix == 1) {weight_img.set_pixel(row,col,1.0/(ntrack+_track_offset));}
        //       else if (pix == 2) {weight_img.set_pixel(row,col,1.0/(nshower+_shower_offset));}
        //       else if (pix == 3) {weight_img.set_pixel(row,col,1.0/(nend+_end_offset));}
        //       else {std::cout << "You've got a problem, pixel labeled something strange!" << std::endl;}
        //
        //     }//End of cols loop (divide by totals)
        //
        //   }//End of rows loop (divide by totals)
        //
        //
        // cropped_weights.emplace_back(weight_img);
        // //std::cout << "emplace1!" << std::endl;
        //
        // }//End of plane loop
        //END WEIGHTED IMAGE making



      	ev_out_adc->emplace( std::move(crop_v) );

      	//ev_out_labels->emplace( std::move(cropped_labels) );
      	// ev_out_weights->emplace( std::move(cropped_weights) );

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

      	foutIO->set_id( run, subrun, 100*event+icrop );
      	foutIO->save_entry();
      	nsaved++;
      }

      if ( _max_images>0 && nsaved>=_max_images )
	     break;
    }
    return true;
  }

  void UBCropMask::cropUsingBBox2D( const larcv::EventBBox2D& bbox_vec,
           const std::vector<larcv::Image2D>& img_v,
           int num_calls,
           larcv::EventImage2D& output_imgs ) {
           // inputs
           // ------
           // bbox_vec, vector of bounding boxes for (u,v,y)
           // img_v, source adc images
           // y1, y2: range of y-wires fully covered by U,V
           // fill_y_image: if true, we fill the entire cropped y-image.
           //               else we only fill the region of y-wires that
           //               are entirely covered by the u,v cropped images
           // minpixfrac: if value is >0, we enforce a minimum value on the
           //               number of pixels occupied in the Y-image
           //
           // outputs
           // --------
           // output_imgs, cropped output image2d instances filled into eventimage2d container

           // get bounding boxes
           const larcv::BBox2D& bbox_u = bbox_vec[0+num_calls*3];
           const larcv::BBox2D& bbox_v = bbox_vec[1+num_calls*3];
           const larcv::BBox2D& bbox_y = bbox_vec[2+num_calls*3];

           //Crop y plane
           larcv::Image2D crop_yimg = img_v[2].crop( bbox_y );
           //Crop u plane
           larcv::Image2D crop_up = img_v[0].crop( bbox_u );
           //Crop z plane
           larcv::Image2D crop_vp = img_v[1].crop( bbox_v );
           output_imgs.emplace( std::move(crop_up) );
           output_imgs.emplace( std::move(crop_vp) );
           output_imgs.emplace( std::move(crop_yimg) );

    return;
  }

  void UBCropMask::AdjustMasksUsingBBox2D( const larcv::EventBBox2D& bbox_vec,
           const std::vector<larcv::ImageMeta>& meta_orig_v,
           const std::vector<larcv::ImageMeta>& meta_crop_v,
           const std::vector<std::vector<larcv::ClusterMask>>& masks_vv,
           int num_calls,
           larcv::EventClusterMask& output_masks
           ) {
           // inputs
           // ------
           // bbox_vec, vector of bounding boxes for (u,v,y)
           // masks_vv, source vector of masks on a plane for each of 3 planes (vector of vector of masks)
           //
           // outputs
           // --------
           // output_masks, adjusted vector of vector of masks based on what the crop contains

           // get bounding boxes
           // const larcv::BBox2D& bbox_u = bbox_vec[0+num_calls*3];
           // const larcv::BBox2D& bbox_v = bbox_vec[1+num_calls*3];
           // const larcv::BBox2D& bbox_y = bbox_vec[2+num_calls*3];

           //Go through the masks in each plane and adjust if they are in the bounding box
           for (int plane=0; plane < masks_vv.size(); plane++){
             //Get bounding box
             const larcv::BBox2D& bbox_crop = bbox_vec[plane+num_calls*3];
             //set up empty vec of masks
             std::vector<ClusterMask> masks_v(0, ClusterMask());
             std::cout << "Total Masks Available in plane "<< plane<< ": "<< masks_vv[plane].size() << "\n";
             for (int m=0; m<masks_vv[plane].size(); m++){
               // std::cout << "Mask\n";
               // std::cout << masks_vv[plane][m].box.min_x() <<"  Min X     " << masks_vv[plane][m].box.max_x()  << "  Max X     \n";
               // std::cout << masks_vv[plane][m].box.min_y() <<"  Min Y     " << masks_vv[plane][m].box.max_y()  << "  Max Y     \n";
               // std::cout << "Crop\n";
               // std::cout << bbox_crop.min_x() <<"  Min X     " << bbox_crop.max_x()  << "  Max X     \n";
               // std::cout << bbox_crop.min_y() <<"  Min Y     " << bbox_crop.max_y()  << "  Max Y     \n";

               // std::cout << masks_vv[plane][m].box.min_x() <<"," << masks_vv[plane][m].box.max_x()<< ","  <<masks_vv[plane][m].box.min_y() <<"," << masks_vv[plane][m].box.max_y()<<","<<bbox_crop.min_x() <<"," << bbox_crop.max_x()<< ","  <<bbox_crop.min_y() <<"," << bbox_crop.max_y()  << "\n";


               if (masks_vv[plane][m].box.min_x() >= bbox_crop.max_x()) {continue;}// Mask all to right of crop, don't use
               if (masks_vv[plane][m].box.max_x() <= bbox_crop.min_x()) {continue;}// Mask all to left of crop, don't use
               if (masks_vv[plane][m].box.min_y() >= bbox_crop.max_y()) {continue;}// Mask all above of crop, don't use
               if (masks_vv[plane][m].box.max_y() <= bbox_crop.min_y()) {continue;}// Mask all below of crop, don't use
               //Okay mask is in crop, let's find it:
               std::cout << "We've Got a good Crop:     ";
               int final_type= -9999;
               //Add Type Counter:
               if (masks_vv[plane][m].type == 13) {
                 nCosmic++;
                 final_type = 0;
               }
               else if (masks_vv[plane][m].type == 2112) {
                 nNeutron++;
                 final_type =1;
               }
               else if (masks_vv[plane][m].type == 2212) {
                 nProton++;
                 final_type =2;
               }
               else if (masks_vv[plane][m].type == 11) {
                 nElectron++;
                 final_type =3;
               }
               else if (masks_vv[plane][m].type == -1) {
                 nNeutrino++;
                 final_type =4;
               }
               else {
                 nOther++;
                 final_type =5;
               }




               //Check if mask is fully contained in crop:
               if (masks_vv[plane][m].box.min_x() >= bbox_crop.min_x() && masks_vv[plane][m].box.min_y() >= bbox_crop.min_y() && masks_vv[plane][m].box.max_x() <= bbox_crop.max_x() && masks_vv[plane][m].box.max_y() <= bbox_crop.max_y()) {
                 std::vector<Point2D> pts_v(0,Point2D());
                 for (int pt=0; pt < masks_vv[plane][m].points_v.size(); pt++ ){
                   double w = meta_orig_v[plane].pos_x(masks_vv[plane][m].points_v[pt].x);
                   double t = meta_orig_v[plane].pos_y(masks_vv[plane][m].points_v[pt].y);
                   double col_new = meta_crop_v[plane].col(w);
                   double row_new = meta_crop_v[plane].row(t);
                   Point2D this_point(col_new, row_new);
                   pts_v.push_back(this_point);
                 }
                 masks_v.push_back(ClusterMask(masks_vv[plane][m].box, meta_crop_v[plane], pts_v, final_type));
                 std::cout << "Fully Contained\n";

               }
               //Mask must be partially obscured by crop, make new mask
               else{
                 std::cout << "Partial Containment:   ";

                 //Create new bbox for adjusted mask
                 double min_x = 99999;
                 double min_y = 99999;
                 double max_x = -1;
                 double max_y = -1;
                 bool box_filled =0;


                 //Determine bounding box coordiantes
                 // if (masks_vv[plane][m].box.min_x() <= bbox_crop.min_x()) {
                 //   min_x = bbox_crop.min_x() ;//- bbox_crop.min_x()/img_v[plane].meta().pixel_width();
                 //   // min_x = 0;
                 // }
                 // else{
                 //   min_x = masks_vv[plane][m].box.min_x() ;//- bbox_crop.min_x()/img_v[plane].meta().pixel_width();
                 // }
                 // if (masks_vv[plane][m].box.max_x() >= bbox_crop.max_x()) {
                 //   max_x = bbox_crop.max_x() ;//- bbox_crop.min_x()/img_v[plane].meta().pixel_width();
                 // }
                 // else{
                 //   max_x = masks_vv[plane][m].box.max_x() ;//- bbox_crop.min_x()/img_v[plane].meta().pixel_width();
                 // }
                 // if (masks_vv[plane][m].box.min_y() <= bbox_crop.min_y()) {
                 //   min_y = bbox_crop.min_y() ;//- bbox_crop.min_y()/img_v[plane].meta().pixel_height();
                 //   // min_y = 0;
                 // }
                 // else{
                 //   min_y = masks_vv[plane][m].box.min_y() ;//- bbox_crop.min_y()/img_v[plane].meta().pixel_height();
                 // }
                 // if (masks_vv[plane][m].box.max_y() >= bbox_crop.max_y()) {
                 //   max_y = bbox_crop.max_y() ;//- bbox_crop.min_y()/img_v[plane].meta().pixel_height();
                 // }
                 // else{
                 //   max_y = masks_vv[plane][m].box.max_y() ;//- bbox_crop.min_y()/img_v[plane].meta().pixel_height();
                 // }

                 // // std::cout << min_x<<  " "<<max_x << " "<< min_y <<" "<< max_y<< "\n\n\n";
                 // larcv::BBox2D adjusted_box(min_x, min_y, max_x, max_y);
                 //Create New vector of points for adjusted mask
                 std::vector<Point2D> pts_v(0,Point2D());
                 for (int pt=0; pt<masks_vv[plane][m].points_v.size() ; pt++){
                   double w = meta_orig_v[plane].pos_x(masks_vv[plane][m].points_v[pt].x);
                   double t = meta_orig_v[plane].pos_y(masks_vv[plane][m].points_v[pt].y);
                   //Check if point in crop
                   //if (bbox_crop.contains(Point2D(w,t))){
                   if (w > bbox_crop.min_x() && t > bbox_crop.min_y() && w < bbox_crop.max_x() && t < bbox_crop.max_y()){
                     if (box_filled == 0){
                       box_filled=1;
                       std::cout << "     FILLED AT LEAST ONCE\n";
                     }
                     if (w < min_x) {
                       min_x = w;
                     }
                     if (w > max_x) {
                       max_x = w;
                     }
                     if (t < min_y) {
                       min_y = t;
                     }
                     if (t > max_y) {
                       max_y = t;
                     }
                     double col_new = meta_crop_v[plane].col(w);
                     double row_new = meta_crop_v[plane].row(t);
                     Point2D this_point(col_new, row_new);
                     pts_v.push_back(this_point);
                   }

                 }//End loop through ancestor points in mask
                 if (box_filled ==0) {
                   std::cout << "\n\n\n\n\n\n\n\nCluster Bounds Crossed, but no points in Crop, Num Masks in plane:    "<< masks_v.size()<< "\n\n\n\n\n\n\n\n";
                   continue;
                 }
                 larcv::BBox2D adjusted_box(min_x, min_y, max_x, max_y);
                 masks_v.push_back(ClusterMask(adjusted_box, meta_crop_v[plane], pts_v, final_type));
                 std::cout << "Num masks in plane "<< plane << ": "<< masks_v.size()<<"\n";
               }


             }//end masks in this plane loop
             //Push back clusters contained in this plane.
             output_masks.emplace(std::move(masks_v));
           }//end planes loop

    return;
  }

  void UBCropMask::finalize()
  {
    std::cout << "\n\nFinalizing!!\n\n\n";

    foutIO->finalize();
  }



}
#endif
