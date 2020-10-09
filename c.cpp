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
#include <wks2mat.h>

////////////////////////////////////////////////////////////////////////////////////

//#pragma labtalk(0) // to disable OC functions for LT calling.

////////////////////////////////////////////////////////////////////////////////////
// Include your own header files here.


////////////////////////////////////////////////////////////////////////////////////
// Start your functions here.

//from matrix extract COLUMNS "from" and "to" and put them in worsheet. Last column in worksheet is the avarage value of each rows.
void extract_columns_from_matrix(int from, int to){
	MatrixLayer mlayer = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject = mlayer.MatrixObjects(1);  // The 1st matrix object, YOU CHANGE THIS NUMBER ACCORDING TO MATRIX YOU WANT TO EXTRAKT FROM
	matrixbase &mbase = mobject.GetDataObject();  // Get data from matrix object
	vector v1;
	
	Worksheet wks;
	wks.Create();
	wks.SetSize(0,0);
	wks.AddCol("y");
	int r=wks.AddCol("sum");
	wks.Columns("y").SetType(OKDATAOBJ_DESIGNATION_X);
	string s = "int((";
	for(int i=from; i<=to; i++){
		wks.AddCol((string)i);
		mbase.GetColumn(v1,i-1);
		int x = (i-from);
		Dataset dsX(wks, i-from + 2);
		dsX.Append(v1);
		s += "Col(" + i + ")+";
	}
	s += "0)/" + (to-from+1) + ")";
	Column col_sum, col_i;
	col_sum.Attach(wks, r);
	col_sum.SetFormula(s, AU_AUTO);
	col_sum.ExecuteFormula();
	col_i.Attach(wks, 0);
	col_i.SetFormula("i", AU_AUTO);
	col_i.ExecuteFormula();
	
}

void get_peak(double a=10, double b=10, int smooth1=1, int smooth2=1){
	MatrixLayer mlayer = Project.ActiveLayer();
	MatrixObject mobject = mlayer.MatrixObjects(0);
	vector<int> z;
	Worksheet wks1, wks2, wks3, result("result");
	GraphPage gp;
	gp.Create("Line");
	vector res;
	z = matrix_to_xyz(mobject, &wks1);
	double min, max;
	int nBinSize_maxvalue = 10000;
	z.GetMinMax(min, max);
	int nBinSize_startvalue = max;
	if (nBinSize_startvalue > nBinSize_maxvalue) nBinSize_startvalue = nBinSize_maxvalue;
	int nBinSize = nBinSize_startvalue;
	int flag = 1;
	while(1){
		vector vBinCenters(nBinSize), vAbsoluteCounts(nBinSize), vCumulativeCounts(nBinSize), normalized_counts(nBinSize);
		frequency_count(z, &wks1, nBinSize, &vBinCenters, &vAbsoluteCounts, &vCumulativeCounts, 0);
		
		double min, max;
		vAbsoluteCounts.GetMinMax(min, max);
		for(int i=0;i<nBinSize;i++) normalized_counts[i] = vAbsoluteCounts[i]/max;
		double pUSS, pCSS;
		ocmath_basic_summary_stats(nBinSize, normalized_counts, NULL, NULL , NULL, NULL , NULL, NULL , NULL, NULL , &pUSS, &pCSS);	//statistics from which the smootnes is determined
		//if(res.GetSize() <= 4) {
		printf("\nBinSize: %d, Peaks: %d", nBinSize, res.GetSize());
		printf("\nBinSize: %d, Uncorrected Sum of Squares: %f", nBinSize, pUSS);
		nBinSize = nBinSize - (nBinSize * 1.0/6.0);
		if(pUSS < a) {
			write_to_worksheet_and_plot(&gp, &wks1, vBinCenters, vAbsoluteCounts, vCumulativeCounts, "countsmatrixone");
			PEAK_vectors(&wks1, vBinCenters, vAbsoluteCounts, 0);
			if (smooth1)	curve_smooth_sg(vAbsoluteCounts, 5,5);
			write_to_worksheet_and_plot(&gp, &wks1, vBinCenters, vAbsoluteCounts, vCumulativeCounts, "matrixonesmoothed");
			res = PEAK_vectors(&wks1, vBinCenters, vAbsoluteCounts, 0);
			break;
		}

	}
	
	printf("\n-----------------------------------------------------------------------------------\n");
	
	z = create_wks3(&wks1, &wks2, &wks3, mlayer, res[0]);					//Create new worksheet from data in wks1 and in matrix2 (variable mobject1)
	nBinSize = nBinSize_startvalue;
	flag=3;
	while(1){
		vector vBinCenters(nBinSize), vAbsoluteCounts(nBinSize), vCumulativeCounts(nBinSize), normalized_counts(nBinSize);			//This is the output data from frequency_count function
		frequency_count(z, &wks2, nBinSize, &vBinCenters, &vAbsoluteCounts, &vCumulativeCounts, 1);	//res is a vector. The size is the number of peak, then this number is 1, the res[0] is the position of the searched peak.
		
		double min, max;
		vAbsoluteCounts.GetMinMax(min, max);
		for(int i=0;i<nBinSize;i++) normalized_counts[i] = vAbsoluteCounts[i]/max;		//first noramilize graph
		double pUSS, pCSS;
		ocmath_basic_summary_stats(nBinSize, normalized_counts, NULL, NULL , NULL, NULL , NULL, NULL , NULL, NULL , &pUSS, &pCSS);	//seconnd calculate statistics on normalized graph
		
		printf("\nBinSize: %d, Peaks: %d", nBinSize, res.GetSize());
		printf("\nBinSize: %d, Uncorrected Sum of Squares: %f, Corrected Sum of Squares: %f", nBinSize, pUSS, pCSS);
		nBinSize = nBinSize - (nBinSize * 1.0/5.0);
		if(pUSS < b) {																									//third if pUSS parameter of statistics met the given value, use bin size
		//if(res.GetSize() <= 1){
			write_to_worksheet_and_plot(&gp, &wks2, vBinCenters, vAbsoluteCounts, vCumulativeCounts, "countsmatrixtwo");
			PEAK_vectors(&wks2, vBinCenters, vAbsoluteCounts, 1);
			if (smooth2)	curve_smooth_sg(vAbsoluteCounts, 5, 5);
			write_to_worksheet_and_plot(&gp, &wks2, vBinCenters, vAbsoluteCounts, vCumulativeCounts, "matrixtwosmoothed");
			PEAK_vectors(&wks2, vBinCenters, vAbsoluteCounts, 1);
			break;
		}
	}
	format_line_plot(&gp);
	FWHM(&wks2, 0, 2, 1);
	float FWHM_value = FWHM(&wks2, 0, 2, 2);
	Dataset dsX(result, 0), dsY(result, 2);
    dsX.Append((vector)res[0]); dsY.Append((vector)FWHM_value);
    //wks1.GetPage().Rename("mat1");
    //wks2.GetPage().Rename("mat1 + mat2");
    //wks3.GetPage().Rename("mat2");
}

vector create_wks3(Worksheet *wks1, Worksheet *wks2, Worksheet *wks3, MatrixLayer mlayer, int peak_pos){
	wks2->Create("Origin");		//Create Worksheet2
	int    iNumCols = wks2->GetNumCols();
    for(int ii = iNumCols ; ii >= 0 ; ii--){
    	wks2->DeleteCol(ii);
    }
	vector<int> vnLabelTypes;
	DWORD dwLabel = RCLT_LONG_NAME;
	vnLabelTypes.Add(dwLabel);
	dwLabel = RCLT_UNIT;
	vnLabelTypes.Add(dwLabel);
	MatrixObject mobject1 = mlayer.MatrixObjects(1);
	matrix_to_xyz(mobject1, wks3);
	wks1->CopyTo(*wks2, 0, -1, 0, -1, 0, -1, CPYT_COPY_COLUMN_FORMAT | CPYT_COPY_COLUMN_DESIGNATIONS | CPYT_COPY_COLUMN_LABELS, 0, &vnLabelTypes);
	wks3->CopyTo(*wks2, 2, -1, 0, -1, 3, -1, CPYT_COPY_COLUMN_FORMAT | CPYT_COPY_COLUMN_DESIGNATIONS | CPYT_COPY_COLUMN_LABELS, 0, &vnLabelTypes);
	wks2->AddCol("E");
	wks2->Columns("C").SetType(OKDATAOBJ_DESIGNATION_Z);
	wks2->Columns("D").SetType(OKDATAOBJ_DESIGNATION_Z);
	wks2->Columns("E").SetType(OKDATAOBJ_DESIGNATION_Z);
	Column colE;
	colE.Attach(*wks2, 4);
	string s = "int(Col(D)/Col(C)*" + peak_pos + ")";	
	colE.SetFormula(s, AU_AUTO);
	colE.ExecuteFormula();
	return colE.GetDataObject(); //return the result from column E
}

// Source matrix is the active layer (Origin matrix window).
//The function creates xyz worksheet with converted matrix. Return is a vector z.
vector matrix_to_xyz(MatrixObject mo, Worksheet *wks)
{
	// Create a MatrixLayer, and set its dimension and number of rows & cols.
    int nRows = mo.GetNumRows();
    int nCols = mo.GetNumCols();

    double dXmin, dXmax, dYmin, dYmax;
    mo.GetXY(dXmin, dYmin, dXmax, dYmax);
    
    Matrix& mat = mo.GetDataObject();
   
    vector<double> a;
    vector<double> b;
    vector<double> c;
    
    a.SetSize(nRows * nCols);
    b.SetSize(nRows * nCols);
    c.SetSize(nRows * nCols);

    // Note: data in mat is in column order. If set bColOrder = false, mat will be transposed and
    // copied to vector c.
    int nRet = ocmath_mat_to_regular_xyz(mat, nRows, nCols, dXmin, dXmax, dYmin, dYmax, a, b, c, false);
	
    if (nRet > 0)
    {
        // Create a worksheet to show XYZ scatters.
        //Worksheet wks;
        wks->Create("Origin");
        wks->SetSize(-1, 3);
        Dataset dsX(*wks, 0), dsY(*wks, 1), dsZ(*wks, 2);
        dsX.Append(a); dsY.Append(b); dsZ.Append(c);
    }    
    return c;
}

vector frequency_count(vector vData, Worksheet *wks,int nBinSize, vector *vBinCenters, vector *vAbsoluteCounts, vector *vCumulativeCounts, int i){
	double min, max;
	vData.Replace(0, NANUM, MATREPL_TEST_EQUAL);
	vData.GetMinMax(min, max);
	//int x = (int)max%nBinSize;			
	FreqCountOptions fcoOptions;
	int x = ((int)max - (int)min)%nBinSize;
	fcoOptions.FromMin = min;//min + (nBinSize - (int)max%nBinSize);
	fcoOptions.ToMax = max+(nBinSize - x);//max + (nBinSize - (int)max%nBinSize);	//this is need for integer bin size.
	fcoOptions.StepSize = nBinSize;
	fcoOptions.IncludeLTMin = 0;
	fcoOptions.IncludeGEMax = 0;
	
	int nOption = FC_NUMINTERVALS; // to extend last bin if not a full bin

	int nRet = ocmath_frequency_count(
    vData, vData.GetSize(), &fcoOptions,
    *vBinCenters, nBinSize, *vAbsoluteCounts, nBinSize,
    *vCumulativeCounts, nBinSize, nOption);

	if( STATS_NO_ERROR == nRet )
		//out_str("Done");
	return 0;
	return PEAK_vectors(wks, *vBinCenters, *vAbsoluteCounts, i);
	
}

//Return number of peaks and positions from vector x,y
vector PEAK_vectors(Worksheet *wks, vector vxData, vector vyData, int i)
{

        
    uint nDataSize = vxData.GetSize();
    int iSize = vxData.GetSize();
    vector vxPeaks, vyPeaks;
    vector<int> vnIndices;
 
    vxPeaks.SetSize(nDataSize);
    vyPeaks.SetSize(nDataSize);
    vnIndices.SetSize(nDataSize);
	int nRet;
    if (i == 0)	nRet = ocmath_find_peaks_1st_derivative( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,false);
	if (i == 1)	nRet = ocmath_find_peaks_by_local_maximum( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,11);
    if( nRet < OE_NOERROR )
    {
        printf("error code1: %d\n", nRet);
        vxPeaks.SetSize(999);	//This is the error value
        return vxPeaks;
    }
    
    double dMinHeight = 10;
    int nPeakNum = nDataSize;
    nRet=ocmath_test_peaks_by_height( &nPeakNum, vxPeaks, vyPeaks, vnIndices, dMinHeight);
     if( nRet < OE_NOERROR )
    {
        printf("error code2: %d\n", nRet);
        vxPeaks.SetSize(999);	//This is the error value
        return vxPeaks;
    }
    vxPeaks.SetSize(nPeakNum);		//delete NANUM values
    vyPeaks.SetSize(nPeakNum);		//delete NANUM values
    vnIndices.SetSize(nPeakNum);	//delete NANUM values
    printf("\nThe peak position for %s is at %.1f number of Peak found is: %d", wks->GetPage().GetName(), vxPeaks[0], vxPeaks.GetSize());
    return vxPeaks;
}


void write_to_worksheet_and_plot(GraphPage *gp, Worksheet *wks, vector vBinCenters, vector vAbsoluteCounts, vector vCumulativeCounts, string s){
	int index = (wks->GetPage()).AddLayer("freq_counts1");
	Worksheet freq_counts1 = (wks->GetPage()).Layers(index);
	int    iNumCols = freq_counts1.GetNumCols();
    for(int ii = iNumCols ; ii >= 0 ; ii--){
    	freq_counts1.DeleteCol(ii);
    }
	freq_counts1.AddCol("BinCenter");
	freq_counts1.Columns("BinCenter").SetType(OKDATAOBJ_DESIGNATION_X);
	freq_counts1.AddCol("BinEnd");
	freq_counts1.AddCol(s);
	freq_counts1.AddCol("CumulCounts");
	Dataset BinCenter(freq_counts1, 0), BinEnd(freq_counts1, 1), Counts(freq_counts1, 2), CumulCounts(freq_counts1, 3);
	BinCenter.Append(vBinCenters), Counts.Append(vAbsoluteCounts); CumulCounts.Append(vCumulativeCounts);
	
	plot(gp, &freq_counts1);
}

void plot(GraphPage *gp, Worksheet *freq_counts1){
	GraphLayer gl = gp->Layers(0);
	Curve crv(*freq_counts1,0,2);
	int nPlot = gl.AddPlot(crv, IDM_PLOT_LINE);
	gl.Rescale();
}

void format_line_plot(GraphPage *gp)
{
    GraphLayer gl = gp->Layers(0);
    if( gl )
    {
		Tree tr;
		tr.Root.Curves.Curve2.Line.Color.nVal = 1; // red
		tr.Root.Curves.Curve3.Line.Color.nVal = 2; // red
		tr.Root.Curves.Curve4.Line.Color.nVal = 3; // red
		if( 0 == gl.UpdateThemeIDs(tr.Root) ) // if no err
		{
			gl.ApplyFormat(tr, true, true);
		}
    }
}

float FWHM(Worksheet *wks, int x, int y, int layer)
{
	// Assumes the active worksheet contains column A and column B with data
	if(wks)
	{
		Worksheet freq_counts = (wks->GetPage()).Layers(layer);
		Curve crv; // Use default constructor to create unattached Origin C Curve object
		crv.Attach( freq_counts, x, y ); // Attach Origin C Curve object to internal Origin data sets
		string ActWksName  = wks->GetPage().GetName();
		vector<int> v;
		//Column cc(wks, 0);
		//vector<ushort> vRowMap(cc);
		//wks.ExtractOneGroup(v,1, vRowMap, 0, 0, vRowMap.GetSize()-1, 1);

		
		float dYoffset = 0.5;
		float dWx = fwhm(crv,dYoffset);  // Demonstration of fwhm
		printf("\nThe peak width of %s at the half max (FWHM)= %g (with offset=%g)\n", ActWksName, dWx, dYoffset);
		return dWx;
	}
}

void create_matrix(){
	MatrixPage matPage;
	matPage.Create("double_matrix");

}



//Convert two matrices in one worksheet to xyzzz format, for the last format set forumula "int(Col(D)/Col(C)*1)
void two_matrices_to_xyzzz(){
	MatrixLayer mlayer = Project.ActiveLayer();
	MatrixObject mobject1 = mlayer.MatrixObjects(0);
	MatrixObject mobject2 = mlayer.MatrixObjects(1);
	Worksheet wks;
	wks.Create();
	wks.SetSize(0,0);
	wks.AddCol("A");
	wks.AddCol("B");
	wks.Columns("A").SetType(OKDATAOBJ_DESIGNATION_X);
	wks.Columns("B").SetType(OKDATAOBJ_DESIGNATION_Y);
	matrix2_to_xyz(mobject1, mobject2, &wks);
	//matrix_to_xyz(mobject, &wks)
	//wks2.AddCol("E");
	//wks2->Columns("C").SetType(OKDATAOBJ_DESIGNATION_Z);
	// Create a MatrixLayer, and set its dimension and number of rows & cols.
    Column colE;
	colE.Attach(wks, 4);
	colE.SetFormula("int(Col(D)/Col(C)*1)", AU_AUTO);
	wks.Columns("C").SetType(OKDATAOBJ_DESIGNATION_Z);
	wks.Columns("D").SetType(OKDATAOBJ_DESIGNATION_Z);
	wks.Columns("E").SetType(OKDATAOBJ_DESIGNATION_Z);
}

//This function is used by function matrix_to_worksheet()
void matrix2_to_xyz(MatrixObject mo1, MatrixObject mo2, Worksheet *wks)
{
	// Create a MatrixLayer, and set its dimension and number of rows & cols.
    int nRows1 = mo1.GetNumRows();	int nRows2 = mo2.GetNumRows();
    int nCols1 = mo1.GetNumCols();	int nCols2 = mo2.GetNumCols();


    double dXmin1, dXmax1, dYmin1, dYmax1, dXmin2, dXmax2, dYmin2, dYmax2;
    mo1.GetXY(dXmin1, dYmin1, dXmax1, dYmax1);    
    mo2.GetXY(dXmin2, dYmin2, dXmax2, dYmax2);
    
    Matrix& mat1 = mo1.GetDataObject();
    Matrix& mat2 = mo2.GetDataObject();
   
    vector<double> a1,a2;
    vector<double> b1,b2;
    vector<double> c1,c2;
    
    a1.SetSize(nRows1 * nCols1);	a2.SetSize(nRows2 * nCols2);
    b1.SetSize(nRows1 * nCols1);	b2.SetSize(nRows2 * nCols2);
    c1.SetSize(nRows1 * nCols1);	c2.SetSize(nRows2 * nCols2);

    // Note: data in mat is in column order. If set bColOrder = false, mat will be transposed and
    // copied to vector c.
    int nRet1 = ocmath_mat_to_regular_xyz(mat1, nRows1, nCols1, dXmin1, dXmax1, dYmin1, dYmax1, a1, b1, c1, false);
	int nRet2 = ocmath_mat_to_regular_xyz(mat2, nRows2, nCols2, dXmin2, dXmax2, dYmin2, dYmax2, a2, b2, c2, false);
	
    if (nRet1 > 0 && nRet2 > 0)
    {
        // Create a worksheet to show XYZ scatters.
        //Worksheet wks;
        wks->SetSize(-1, 5);
        Dataset dsX(*wks, 0), dsY(*wks, 1), dsZ1(*wks, 2), dsZ2(*wks, 3);
        dsX.Append(a1); dsY.Append(b1); dsZ1.Append(c1);	dsZ2.Append(c2);
    }
    else printf("Error in matrix2_to_xyz");
}



void eltol(int from, int to){
	MatrixLayer mlayer1 = Project.ActiveLayer();  // Active matrix sheet
	MatrixObject mobject1 = mlayer1.MatrixObjects(1);  // The 1st matrix object
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
	v2.SetSize(col_size);
	for(int i=0; i<col_size; i++){
		mbase1.GetColumn(v1,i);
		v2.SetSubVector(v1,abs(128-i));
		v1.SetSize(col_size-i);
		mbase2.SetColumn(v1,i);
		}
	return;
	
	
}








//These functions are not used!!!

void FWHM(int x, int y)
{
	// Assumes the active worksheet contains column A and column B with data
	Worksheet wks = Project.ActiveLayer();
	if(wks)
	{
		Curve crv; // Use default constructor to create unattached Origin C Curve object
		crv.Attach( wks, x, y ); // Attach Origin C Curve object to internal Origin data sets
		string ActWksName  = wks.GetPage().GetName();
		vector<int> v;
		//Column cc(wks, 0);
		//vector<ushort> vRowMap(cc);
		//wks.ExtractOneGroup(v,1, vRowMap, 0, 0, vRowMap.GetSize()-1, 1);

		
		float dYoffset = 0.5;
		float dWx = fwhm(crv,dYoffset);  // Demonstration of fwhm
		printf("The peak width of %s at the half max (FWHM)= %g (with offset=%g)\n", ActWksName, dWx, dYoffset);
	}
}

//Return a number of found peaks and its positions from curve from active worksheet column x and y (PEAK_worksheet(x,y))
void PEAK_worksheet(int column_number_x, int column_number_y)
{
	Worksheet wks = Project.ActiveLayer();
    Column col1 = wks.Columns(column_number_x);  // 1st column
	Column col2 = wks.Columns(column_number_y);  // 2nd column
    vector vxData(col1, -1, -1, WRITEBACK_DELETE_ON_SHRINK | WRITEBACK_INSERT_ON_EXPAND); //transform column to vector
    vector vyData(col2, -1, -1, WRITEBACK_DELETE_ON_SHRINK | WRITEBACK_INSERT_ON_EXPAND);
        
    uint nDataSize = vxData.GetSize();
    int iSize = vxData.GetSize();
    vector vxPeaks, vyPeaks;
    vector<int> vnIndices;
 
    vxPeaks.SetSize(nDataSize);
    vyPeaks.SetSize(nDataSize);
    vnIndices.SetSize(nDataSize);
 
    //int nRet = ocmath_find_peaks_2nd_derivative( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,false);
	int nRet = ocmath_find_peaks_by_local_maximum( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,11);
    if( nRet < OE_NOERROR )
    {
        printf("error code: %d\n", nRet);
        return;
    }
    
    double dMinHeight = 10;
    int nPeakNum = nDataSize;
    nRet=ocmath_test_peaks_by_height( &nPeakNum, vxPeaks, vyPeaks, vnIndices, dMinHeight);
     if( nRet < OE_NOERROR )
    {
        printf("error code: %d\n", nRet);
        return;
    }
    vxPeaks.SetSize(nPeakNum);		//delete NANUM values
    vyPeaks.SetSize(nPeakNum);		//delete NANUM values
    vnIndices.SetSize(nPeakNum);	//delete NANUM values
    printf("The peak position for %s is at %.1f number of Peak found is: %d", wks.GetPage().GetName(), vxPeaks[0], vxPeaks.GetSize());
}

//This function writes found peaks into new worksheet and plots it into the analyzed graph. Graph must be selected!!!
void PEAK_graph()
{
    GraphLayer gl = Project.ActiveLayer();
    if (!gl)
    {
        return;
    }
 
    //get data from the first data plot
    DataPlot dp = gl.DataPlots(0);        
    DataRange dr;
    vector vxData, vyData;
    if(dp.GetDataRange(dr))
    {
        DWORD dwPlotID;
        if(dr.GetData(DRR_GET_DEPENDENT | DRR_NO_FACTORS, 0, &dwPlotID, NULL, &vyData, &vxData) < 0)
        {
            printf("get_plot_data failed GetData");
            return;
        }
    }
 
    uint nDataSize = vxData.GetSize();
    int iSize = vxData.GetSize();
 
    vector vxPeaks, vyPeaks;
    vector<int> vnIndices;
 
    vxPeaks.SetSize(nDataSize);
    vyPeaks.SetSize(nDataSize);
    vnIndices.SetSize(nDataSize);
 
   //int nRet = ocmath_find_peaks_2nd_derivative( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,false);
	int nRet = ocmath_find_peaks_by_local_maximum( &nDataSize, vxData, vyData, vxPeaks, vyPeaks, vnIndices, POSITIVE_DIRECTION | NEGATIVE_DIRECTION,11);
    if( nRet < OE_NOERROR )
    {
        printf("error code: %d\n", nRet);
        return;
    }
    
    double dMinHeight = 10;
    int nPeakNum = nDataSize;
    nRet=ocmath_test_peaks_by_height( &nPeakNum, vxPeaks, vyPeaks, vnIndices, dMinHeight);
     if( nRet < OE_NOERROR )
    {
        printf("error code: %d\n", nRet);
        return;
    }
    vxPeaks.SetSize(nPeakNum);		//delete NANUM values
    vyPeaks.SetSize(nPeakNum);		//delete NANUM values
    vnIndices.SetSize(nPeakNum);	//delete NANUM values
    
    //new a worksheet to output the result
    WorksheetPage wksPage;
    wksPage.Create("Origin");
    Worksheet wksResult = wksPage.Layers(0);
    int nIndCol, nXCol, nYCol;
    nIndCol = wksResult.AddCol("Indices");
    nXCol = wksResult.AddCol("X Coordinate");
    nYCol = wksResult.AddCol("Y Coordinate");
    wksResult.Columns(nIndCol).SetType(OKDATAOBJ_DESIGNATION_X);
    wksResult.Columns(nXCol).SetType(OKDATAOBJ_DESIGNATION_X);
    wksResult.Columns(nYCol).SetType(OKDATAOBJ_DESIGNATION_Y);
    
    DataRange drOut;
    drOut.Add("X", wksResult, 0, nIndCol, -1, nIndCol);
    drOut.Add("Y", wksResult, 0, nXCol, -1, nXCol);
    drOut.Add("Z", wksResult, 0, nYCol, -1, nYCol);
    drOut.SetData(&vyPeaks, &vxPeaks, &vnIndices);
    
    //show the result in the data plot
    XYRange plotRange;
    plotRange.Add("X", wksResult, 0, nXCol, -1, nXCol);
    plotRange.Add("Y", wksResult, 0, nYCol, -1, nYCol);
    gl.AddPlot(plotRange, IDM_PLOT_SCATTER);
}



//From there only examples

void Worksheet_ExtractOneGroup_Ex1(int nGroupVal = 1, int nGroupCol = 0, int nDataCol = 1)
{
    Worksheet wks = Project.ActiveLayer();
    if(!wks)
        return;
    
    Column cc(wks, nGroupCol);
    vector<ushort> vRowMap(cc);// assume group col already contain group numeric values like 1,2,3 etc
    vector vResult;
    if(wks.ExtractOneGroup(vResult, nDataCol, vRowMap, 0, 0, -1, nGroupVal))
    {
        // save result into new wks to view
        Worksheet wResult("Test");
        if(wResult == NULL)
        {
            wResult.Create("Origin");
            //wResult.GetPage().Rename("Test");
            wResult.SetSize(-1, 0);
        }
        int nCol = wResult.AddCol();
        DataRange dr;
        dr.Add("X", wResult, 0, nCol, -1, nCol);
        dr.SetData(vResult);
    }
    else
        printf("Fail to extract grouped vectors from column in worksheet!");
}

int Worksheet_GetSelectedColumns_Ex1()
{
	Worksheet wks1 = Project.ActiveLayer();
	WorksheetPage wksPage = wks1.GetPage();
    //WorksheetPage wksPage = Project.WorksheetPages(0);
    if(!wksPage)
        return -1;
    
    Worksheet wks(wksPage.GetName());
    vector<int> v;
    BOOL bOK = wks.GetSelectedColumns(v);
    if(bOK)
    {
        for( int ii = 0 ; ii < v.GetSize() ; ii++ )
            printf("Column %u has a selection\n", v[ii] + 1);
        return (int) v.GetSize();
    }
    else
        return -1;
}