
void make_axis(double, double, double, double, double, double);
void make_point(double, double, double, TGeoVolume*, TGeoVolume*);

//_____________________________________________________________________________
void run_geo() {

  //initialize
  gSystem->Load("libGeom");
  TGeoManager *geo = new TGeoManager("geo", "geo");
  TGeoMedium *med = new TGeoMedium("mat",1, new TGeoMaterial("mat"));
  TGeoVolume *top = geo->MakeBox("top", med, 100, 100, 100);
  geo->SetTopVolume(top);

  TGeoVolume *p = geo->MakeSphere("p", med, 0, 0.18);
  p->SetLineColor(kGreen);
  p->SetTransparency(0);



  FiberYZ fib(10, 8, 1.1);

  //slice loop
  for(size_t i=0; i<fib.GetNSlice(); i++) {

    FiberYZ::slice& slc = fib.GetSlice(i);

    //point loop
    for(size_t ip=0; ip<slc.points.size(); ip++) {

      //cout <<slc.points.at(ip).at(2) << " " << slc.points.at(ip).at(1) << " " << slc.points.at(ip).at(0) << endl;

      make_point(slc.points.at(ip).at(2), slc.points.at(ip).at(1), slc.points.at(ip).at(0)-1, p, top);
    }//point loop
  }//slice loop


  //end of geometry definition
  geo->CloseGeometry();

  //add graphics and text objects
  TEveManager::Create();

  TEveGeoTopNode *node = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  node->SetVisLevel(4);
  node->GetNode()->GetVolume()->SetVisibility(kFALSE);
  gEve->AddGlobalElement(node);

  //z-axis
  make_axis(12, 0, 0, -1, 0, 0);

  //y-axis
  make_axis(0, 12, 0, 0, -1, 0);

  //x-axis
  make_axis(0, 0, -4, 0, 0, 2);

  //finish and draw
  gEve->FullRedraw3D(kTRUE);

}//run_geo

//_____________________________________________________________________________
void make_point(double x, double y, double z, TGeoVolume *p, TGeoVolume *top) {

  top->AddNode(p, 1, new TGeoTranslation(x, y, z) );

}//make_point

//_____________________________________________________________________________
void make_axis(double x, double y, double z, double xs, double ys, double zs) {

  double rL = 0.12;

  double L = TMath::Sqrt( (x-xs)*(x-xs) + (y-ys)*(y-ys) + (z-zs)*(z-zs) );

  TEveArrow *axis = new TEveArrow(x, y, z, xs, ys, zs);
  axis->SetMainColor(kOrange-3);
  axis->SetTubeR(rL/L);
  axis->SetConeL(5*rL/L);
  axis->SetConeR(3*rL/L);
  gEve->AddGlobalElement(axis);

/*
  TEveText *label = new TEveText("x");
  label->SetFontSize( 30 );
  label->SetMainColor( kGreen );
  label->RefMainTrans().SetPos(x, y, z-4);
  gEve->AddGlobalElement(label);
*/

}//make_axis































