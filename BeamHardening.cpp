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

//Extract BeamHardening curve from defined x and y pixels
//Matrices must be in an active MatrixObject in different layers. Layers LongName must represent a thickness of irradiated sample
void Extract_beam_Hardening_curve(int x, int y){
	int i=0;
	vector<int> vector_x, vector_y,vm;
	MatrixLayer mlayer = Project.ActiveLayer();  // Active matrix sheet
    while(1){
		if(!mlayer.SetActive(i)) break;			//After we select the last active object break
		Matrix m(mlayer);
    	MatrixObject mobject = mlayer.MatrixObjects(i);
		vector_x.Add(atoi(mobject.GetLongName()));
		vector_y.Add(m.GetCellValue(x,y));
    	i++;
    }
	
    Worksheet wksResult;
    wksResult.Create("Origin");
    DataRange dr;
    dr.Add(wksResult,0,"X");
    dr.Add(wksResult,1,"Y");
    dr.SetData(&vector_y,&vector_x);
}

//
void Beam_Hardening_correction_linear(float koeff){
	MatrixLayer activelayer = Project.ActiveLayer();  							//Get active layer
	MatrixObject RAW_object = activelayer.MatrixObjects(-1);
	matrix & RAW_matrix=RAW_object.GetDataObject();
	//int activelayer_index = activelayer.GetActive();
	MatrixPage matPage = activelayer.GetPage();									//Get active page from layer
	MatrixLayer first_sheet = matPage.Layers(0);								//Acess first sheet in the page
	MatrixLayer second_sheet = matPage.Layers(1);								//Acess second sheet in the page
	int BH_corrected_matrix = first_sheet.Insert(1, -1, -1);					//Insert new matrix
	MatrixObject mobject = first_sheet.MatrixObjects(BH_corrected_matrix);		//Select this new matrix object
	matrix & mat=mobject.GetDataObject();										//connect it to a matrix two dimensional array, This is the corrected matrix
	
	vector thickness,cc;
	int l=0;
	int m[10];
	
	while(1){
		if(!second_sheet.SetActive(l)) break;			//After we select the last active object break
    	MatrixObject BH_object = second_sheet.MatrixObjects(l);
		thickness.Add(atoi(BH_object.GetLongName()));
    	l++;
    }
	
	//Linear approximation
	float c,y0;
	for(int i=0;i<256;i++){
		for(int j=0;j<256;j++){
			for(int k=0;l>k;k++){
				MatrixObject BH_object = second_sheet.MatrixObjects(k);
				matrix & BH_matrix=BH_object.GetDataObject();
				float d1=koeff*RAW_matrix[i][j];	//for debug
				float d2=BH_matrix[i][j];			//for debug
				if(k==l){
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					c=(BH_matrix[i][j]-BH_matrixx[i][j])/(thickness[k]-thickness[k-1]);
					y0=BH_matrix[i][j]-c*thickness[k];
					float d3=BH_matrixx[i][j];		//for debug
 					break;
				}				
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k == 0){
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k+1);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					float f = BH_matrixx[i][j];
					c=(BH_matrixx[i][j]-BH_matrix[i][j])/(thickness[k+1]-thickness[k]);
					y0=BH_matrix[i][j]-c*thickness[k];
					float d3=BH_matrixx[i][j];		//for debug
					break;
				}
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k != l){
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					float f = BH_matrixx[i][j];
					c=(BH_matrix[i][j]-BH_matrixx[i][j])/(thickness[k]-thickness[k-1]);
					y0=BH_matrix[i][j]-c*thickness[k];
					float d3=BH_matrixx[i][j];		//for debug
					break;
				}

			}
			mat[i][j]=(koeff*RAW_matrix[i][j]-y0)/c;
		}
	}
}

//
void Beam_Hardening_correction_exponential(float koeff){

	MatrixLayer activelayer = Project.ActiveLayer();  							//Get active layer
	MatrixObject RAW_object = activelayer.MatrixObjects(-1);
	matrix & RAW_matrix=RAW_object.GetDataObject();
	//int activelayer_index = activelayer.GetActive();
	MatrixPage matPage = activelayer.GetPage();									//Get active page from layer
	MatrixLayer first_sheet = matPage.Layers(0);								//Acess first sheet in the page
	MatrixLayer second_sheet = matPage.Layers(1);								//Acess second sheet in the page
	int BH_corrected_matrix = first_sheet.Insert(1, -1, -1);					//Insert new matrix
	MatrixObject mobject = first_sheet.MatrixObjects(BH_corrected_matrix);		//Select this new matrix object
	matrix & mat=mobject.GetDataObject();										//connect it to a matrix two dimensional array, This is the corrected matrix
	
	vector thickness,cc;
	int l=0;
	int m[10];
	
	float a,A,O;
	float x1,x2,y1,y2;
	while(1){
		if(!second_sheet.SetActive(l)) break;			//After we select the last active object break
    	MatrixObject BH_object = second_sheet.MatrixObjects(l);
		thickness.Add(atoi(BH_object.GetLongName()));
    	l++;
    }
	//exponential approximation
	float c,y0;
	for(int i=0;i<256;i++){
		for(int j=0;j<256;j++){
			for(int k=0;l>k;k++){
				MatrixObject BH_object = second_sheet.MatrixObjects(k);
				matrix & BH_matrix=BH_object.GetDataObject();
				MatrixObject BH_object_last = second_sheet.MatrixObjects(l-1);
				matrix & BH_matrix_last=BH_object_last.GetDataObject();
				
				float d1=koeff*RAW_matrix[i][j];	//for debug
				float d2=BH_matrix[i][j];			//for debug
				if(k==(l-1)){ //its is not number 1 its the letter 1!!!! 				If RAW matrix values are smaller than any values from BH matrixes
					MatrixObject BH_objectx_y1 = second_sheet.MatrixObjects(k-2);
					matrix & BH_matrix_y1=BH_objectx_y1.GetDataObject();
					MatrixObject BH_objectx_y2 = second_sheet.MatrixObjects(k-1);
					matrix & BH_matrix_y2=BH_objectx_y2.GetDataObject();
					
					O=BH_matrix[i][j];
					x1=thickness[k-2], x2=thickness[k-1];
					y1=BH_matrix_y1[i][j], y2=BH_matrix_y2[i][j];
					float ln1=ln(y1-O);
					float ln2=ln(y2-O);
					A=exp((((x1/x2)*ln2)-ln1)/((x1/x2)-1));
					float ln3=ln((y2-O)/A);
					a=(1/x2)*ln3;
					break;
				}				
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k == 0){				//If RAW matrix values are bigger than any values from BH matrixes
					MatrixObject BH_objectx = second_sheet.MatrixObjects(3);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					MatrixObject BH_object_next = second_sheet.MatrixObjects(k+1);
					matrix & BH_matrix_next=BH_object_next.GetDataObject();
					O=BH_matrixx[i][j];
					x1=thickness[0], x2=thickness[1];
					y1=BH_matrix[i][j], y2=BH_matrix_next[i][j];
					float ln1=ln(y1-O);
					float ln2=ln(y2-O);
					A=exp((((x1/x2)*ln2)-ln1)/((x1/x2)-1));
					float ln3=ln((y2-O)/A);
					a=(1/x2)*ln3;
					break;
				}
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k != l){				//This is case when RAW matrix values are between BH matrices
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
					matrix & BH_matrix_previous=BH_objectx.GetDataObject();
					MatrixObject BH_object_next = second_sheet.MatrixObjects(k+1);
					matrix & BH_matrix_next=BH_object_next.GetDataObject();
					if(k<l-2){		//If the O point (the constant background) is the nearest point int the BH curve (the last point is where is the smallest counts is)
						MatrixObject BH_object_last = second_sheet.MatrixObjects(k+2);
						matrix & BH_matrix_last=BH_object_last.GetDataObject();
						O=BH_matrix_last[i][j];
					}
					else O=BH_matrix_last[i][j];
					//O=BH_matrix_last[i][j];		//If the O point (the constant background) is the last point in the BH curve (last point is where the smallest counts is)
					x1=thickness[k-1], x2=thickness[k];
					y1=BH_matrix_previous[i][j], y2=BH_matrix[i][j];
					float thickness_ratio=x1/x2;
					float ln1=ln(y1-O);
					float ln2=ln(y2-O);
					
					float B=ln1-ln2;
					float e=((x1/x2)*ln2-ln1)/((x1/x2)-1);
					A=exp((((x1/x2)*ln2)-ln1)/((x1/x2)-1));
					float ln3=ln((y2-O)/A);
					a=(1/x2)*ln3;
					break;
				}

			}
			float y=RAW_matrix[i][j];
			mat[i][j]=(1/a)*ln((y-O)/A);
			float result = mat[i][j];
		}
	}
}

//
void Beam_Hardening_correction_exponential2(float koeff){		//from paper Data processing and image reconstruction

	MatrixLayer activelayer = Project.ActiveLayer();  							//Get active layer
	MatrixObject RAW_object = activelayer.MatrixObjects(-1);
	matrix & RAW_matrix=RAW_object.GetDataObject();
	//int activelayer_index = activelayer.GetActive();
	MatrixPage matPage = activelayer.GetPage();									//Get active page from layer
	MatrixLayer first_sheet = matPage.Layers(0);								//Acess first sheet in the page
	MatrixLayer second_sheet = matPage.Layers(1);								//Acess second sheet in the page
	int BH_corrected_matrix = first_sheet.Insert(1, -1, -1);					//Insert new matrix
	MatrixObject mobject = first_sheet.MatrixObjects(BH_corrected_matrix);		//Select this new matrix object
	matrix & mat=mobject.GetDataObject();										//connect it to a matrix two dimensional array, This is the corrected matrix
	
	vector thickness,cc;
	int l=0;
	int m[10];
	
	float a1,A1,a2,A2,O;
	float x1,x2,x3,y,y1,y2,y3;
	while(1){
		if(!second_sheet.SetActive(l)) break;			//After we select the last active object break
    	MatrixObject BH_object = second_sheet.MatrixObjects(l);
		thickness.Add(atoi(BH_object.GetLongName()));
    	l++;
    }
	//Linear approximation
	float c,y0;
	for(int i=0;i<256;i++){
		for(int j=0;j<256;j++){
			for(int k=0;l>k;k++){
				MatrixObject BH_object = second_sheet.MatrixObjects(k);
				matrix & BH_matrix=BH_object.GetDataObject();

				float d1=koeff*RAW_matrix[i][j];	//for debug
				float d2=BH_matrix[i][j];			//for debug
				if(k==l){ //its is not number 1 its the letter 1!!!! 				If RAW matrix values are smaller than any values from BH matrixes
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					c=(BH_matrix[i][j]-BH_matrixx[i][j])/(thickness[k]-thickness[k-1]);
					y0=BH_matrix[i][j]-c*thickness[k];
					float d3=BH_matrixx[i][j];		//for debug
 					break;
				}				
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k == 0){				//If RAW matrix values are bigger than any values from BH matrixes
					MatrixObject BH_objectx = second_sheet.MatrixObjects(k+1);
					matrix & BH_matrixx=BH_objectx.GetDataObject();
					float f = BH_matrixx[i][j];
					MatrixObject BH_object_next = second_sheet.MatrixObjects(k+1);
					matrix & BH_matrix_next=BH_object_next.GetDataObject();
					O=BH_matrixx[i][j];
					x1=thickness[k], x2=thickness[k+1];
					y1=BH_matrix[i][j], y2=BH_matrix_next[i][j];
					float thickness_ratio=thickness[k]/thickness[k+1];;
					A1=exp((((x1/x2)*ln(y2-O))-ln(y1-O))/((x1/x2)-1));
					a1=(1/x2)*ln((y2-O)/A1);
					float d3=BH_matrixx[i][j];		//for debug
					break;
				}
				if(koeff*RAW_matrix[i][j]>BH_matrix[i][j] && k != 1 && k<l-2){				//This is case when RAW matrix values are between BH matrices
					if (k<l-2){
						MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
						matrix & BH_matrix_previous=BH_objectx.GetDataObject();
						MatrixObject BH_object_next = second_sheet.MatrixObjects(k+1);
						matrix & BH_matrix_next=BH_object_next.GetDataObject();
						MatrixObject BH_object_last = second_sheet.MatrixObjects(k+2);
						matrix & BH_matrix_last=BH_object_last.GetDataObject();
						O=BH_matrix_last[i][j];
						//x1=thickness[k-1], x2=thickness[k];
						x1=thickness[k-1], x2=thickness[k]; x3=thickness[k+1];
						
						y1=BH_matrix_previous[i][j], y2=BH_matrix[i][j], y3=BH_matrix_next[i][j];
						float thickness_ratio=x1/x2;
						float ln1=ln(y1-y3);
						float ln2=ln(y2-y3);
						
						float B=ln1-ln2;
						float e=((x1/x2)*ln2-ln1)/((x1/x2)-1);
						A1=exp((((x1/x2)*ln2)-ln1)/((x1/x2)-1));
						float ln3=ln((y2-O)/A1);
						a1=(1/x2)*ln3;
						ln1=ln(y2-O);
						ln2=ln(y3-O);
						A2=exp((((x2/x3)*ln2)-ln1)/((x2/x3)-1));
						ln3=ln((y3-O)/A2);
						a2=(1/x3)*ln3;
						break;
					}
					else {
						MatrixObject BH_objectx = second_sheet.MatrixObjects(k-1);
						matrix & BH_matrix_previous=BH_objectx.GetDataObject();
						MatrixObject BH_object_next = second_sheet.MatrixObjects(k+1);
						matrix & BH_matrix_next=BH_object_next.GetDataObject();
						x1=thickness[k-1], x2=thickness[k];
						y1=BH_matrix_previous[i][j], y2=BH_matrix[i][j];
						float thickness_ratio=x1/x2;
						float ln1=ln(y1-O);
						float ln2=ln(y2-O);
						
						float B=ln1-ln2;
						float e=((x1/x2)*ln2-ln1)/((x1/x2)-1);
						A1=exp((((x1/x2)*ln2)-ln1)/((x1/x2)-1));
						float ln3=ln((y2-O)/A1);
						a1=(1/x2)*ln3;
						break;
					}
				}

			}
			//mat[i][j]=(1/a)*ln((y-O)/A);
			y=RAW_matrix[i][j];
			if(k<l-2)	mat[i][j]=((y-y2)/(a2*(y1-y2)))*ln((y-O)/A2) + ((y1-y)/(a1*(y1-y2)))*ln((y-y3)/A1);
			else		mat[i][j]=(1/a1)*ln((RAW_matrix[i][j]-O)/A1);
			//float x1=((y-y2)/(a2*(y1-y2)))*ln((y-O)/A2);
			//float x2=((y1-y)/(a1*(y1-y2)))*ln((y-y3)/A1);
			float result = mat[i][j];
		}
	}
}
