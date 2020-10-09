/*------------------------------------------------------------------------------*
 * File Name:				 													*
 * Creation: 																	*
 * Purpose: OriginC Source C file												*
 * Copyright (c) ABCD Corp.	2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010		*
 * All Rights Reserved															*
 * 																				*
 * Modification Log:															*
 *------------------------------------------------------------------------------*/
 
////////////////////////////////////////////////////////////////////////////////////
// Including the system header file Origin.h should be sufficient for most Origin
// applications and is recommended. Origin.h includes many of the most common system
// header files and is automatically pre-compiled when Origin runs the first time.
// Programs including Origin.h subsequently compile much more quickly as long as
// the size and number of other included header files is minimized. All NAG header
// files are now included in Origin.h and no longer need be separately included.
//
// Right-click on the line below and select 'Open "Origin.h"' to open the Origin.h
// system header file.
#include <Origin.h>
////////////////////////////////////////////////////////////////////////////////////

//#pragma labtalk(0) // to disable OC functions for LT calling.

////////////////////////////////////////////////////////////////////////////////////
// Include your own header files here.

////////////////////////////////////////////////////////////////////////////////////
// Start your functions here.


//Use this function for CW - clockwise rotation of active matrix layer
void skew_cw(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	skew_1(angle, mbase1, mbase2);
	skew_4(angle, mbase2);	
}

//Use this function for CC - counterclockwise rotation of active matrix layer
void skew_cc(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	skew_2(angle, mbase1, mbase2);
	skew_3(angle, mbase2);
}

void skew_1(int angle, matrixbase &mbase1, matrixbase &mbase2){
	angle = 45 / angle;
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;

	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(col_size);
		mbase1.GetColumn(v1,i);
		v2.SetSubVector(v1, abs((128-(i))/angle));
		mbase2.SetColumn(v2,i);
		
		mbase1.GetColumn(v1,255-i);
		v1.GetSubVector(v2, abs((128-(i))/angle));
		mbase2.SetColumn(v2,255-i);
	}	
	return;	
}

void skew_2(int angle, matrixbase &mbase1, matrixbase &mbase2){
	angle = 45 / angle;
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(col_size);
		mbase1.GetColumn(v1,255-i);
		v2.SetSubVector(v1, abs((128-(i))/angle));
		mbase2.SetColumn(v2,255-i);
		
		mbase1.GetColumn(v1,i);
		v1.GetSubVector(v2, abs((128-(i))/angle));
		mbase2.SetColumn(v2,i);
	}
	return;
}

void skew_3(int angle, matrixbase &mbase2){
	angle = 45 / angle;
	int col_size = mbase2.GetNumCols();
	int row_size = mbase2.GetNumRows();
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(row_size);
		mbase2.GetRow(v1,i);
		v2.SetSubVector(v1, abs((128-(i))/angle));
		mbase2.SetRow(v2,i);
		
		mbase2.GetRow(v1,255-i);
		v1.GetSubVector(v2, abs((128-(i))/angle));
		mbase2.SetRow(v2,255-i);
	}
	return;
}

void skew_4(int angle, matrixbase &mbase2){
	angle = 45 / angle;
	int col_size = mbase2.GetNumCols();
	int row_size = mbase2.GetNumRows();
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(row_size);
		mbase2.GetRow(v1,255-i);
		v2.SetSubVector(v1, abs((128-(i))/angle));
		mbase2.SetRow(v2,255-i);	
		mbase2.GetRow(v1,i);
		v1.GetSubVector(v2, abs((128-(i))/angle));
		mbase2.SetRow(v2,i);
	}
	return;	
}

//Use theese 4 functions for skew matrix lines of active matrix layer
void skew_1(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;

	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(col_size);
		mbase1.GetColumn(v1,i);
		v2.SetSubVector(v1, abs((128-(i))/3));
		mbase2.SetColumn(v2,i);
		
		mbase1.GetColumn(v1,255-i);
		v1.GetSubVector(v2, abs((128-(i))/3));
		mbase2.SetColumn(v2,255-i);
	}	
	return;	
}

void skew_2(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(col_size);
		mbase1.GetColumn(v1,255-i);
		v2.SetSubVector(v1, abs((128-(i))/3));
		mbase2.SetColumn(v2,255-i);
		
		mbase1.GetColumn(v1,i);
		v1.GetSubVector(v2, abs((128-(i))/3));
		mbase2.SetColumn(v2,i);
	}
	return;
}

void skew_3(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(row_size);
		mbase1.GetRow(v1,i);
		v2.SetSubVector(v1, abs((128-(i))/3));
		mbase2.SetRow(v2,i);
		
		mbase1.GetRow(v1,255-i);
		v1.GetSubVector(v2, abs((128-(i))/3));
		mbase2.SetRow(v2,255-i);
	}
	return;	
}

void skew_4(int angle){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(0);  // The 1st matrix object
	matrixbase &mbase1 = mobject1.GetDataObject();  // Get data from matrix object
	// new a matrix window
	MatrixPage matPage;
	matPage.Create("Origin");
	MatrixLayer mlayer2 = matPage.Layers(); // get active matrix sheet
	MatrixObject mobject2 = mlayer2.MatrixObjects(0);  // Get matrix object
	matrixbase &mbase2 = mobject2.GetDataObject();  // Get data object
	int col_size = mbase1.GetNumCols();
	int row_size = mbase1.GetNumRows();
	mbase2.SetSize(row_size, col_size);  // Set size 2 rows x 3 columns
	vector v1,v2;
	int i;
	
	v2.SetSize(col_size);
	for(i=0; i<col_size/2; i++){
		v2.RemoveAll();
		v2.SetSize(row_size);
		mbase1.GetRow(v1,255-i);
		v2.SetSubVector(v1, abs((128-(i))/3));
		mbase2.SetRow(v2,255-i);	
		mbase1.GetRow(v1,i);
		v1.GetSubVector(v2, abs((128-(i))/3));
		mbase2.SetRow(v2,i);
	}
	return;
	
	
}