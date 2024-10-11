
//_____________________________________________________________________________
//
// Tracks collection
//
//_____________________________________________________________________________

#include "TagTracks.h"

//_____________________________________________________________________________
TagTracks::Coll::Coll() {

  //data members of Track to be written in the output from Coll of Track objects
  DATA_ADD_UNIT_ATTR( x )
  DATA_ADD_UNIT_ATTR( y )
  DATA_ADD_UNIT_ATTR( theta_x )
  DATA_ADD_UNIT_ATTR( theta_y )
  DATA_ADD_UNIT_ATTR( chi2_x )
  DATA_ADD_UNIT_ATTR( chi2_y )
  DATA_ADD_UNIT_ATTR( chi2_xy )
  DATA_ADD_UNIT_ATTR( is_prim )
  DATA_ADD_UNIT_ATTR( itrk )
  DATA_ADD_UNIT_ATTR( pdg )
  DATA_ADD_UNIT_ATTR( prim_id )
  DATA_ADD_UNIT_ATTR( num_diff_itrk )

  DATA_ADD_UNIT_ATTR( is_rec )
  DATA_ADD_UNIT_ATTR( rec_en )
  DATA_ADD_UNIT_ATTR( rec_theta )
  DATA_ADD_UNIT_ATTR( rec_phi )
  DATA_ADD_UNIT_ATTR( rec_en_err )
  DATA_ADD_UNIT_ATTR( rec_theta_err )
  DATA_ADD_UNIT_ATTR( rec_phi_err )
  DATA_ADD_UNIT_ATTR( rec_Q2 )
  DATA_ADD_UNIT_ATTR( ninp )
  DATA_ADD_UNIT_ATTR( ilay )

  DATA_ADD_UNIT_ATTR( has_cal )
  DATA_ADD_UNIT_ATTR( cal_x )
  DATA_ADD_UNIT_ATTR( cal_y )
  DATA_ADD_UNIT_ATTR( cal_en )
  DATA_ADD_UNIT_ATTR( cal_ext_x )
  DATA_ADD_UNIT_ATTR( cal_ext_y )

  DATA_ADD_UNIT_ATTR( has_mcp )
  DATA_ADD_UNIT_ATTR( mcp_en )
  DATA_ADD_UNIT_ATTR( mcp_theta )
  DATA_ADD_UNIT_ATTR( mcp_phi )
  DATA_ADD_UNIT_ATTR( mcp_Q2 )

}//Coll

