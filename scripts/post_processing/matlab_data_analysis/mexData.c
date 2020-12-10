//mexData.C
//   can efficiently extract data from a file which contains ASCII data and text.
//   written by Mike.Zhang, 5-5/28/2002, Copyright reserved.
//   run under MatLab enviroment and compiled by MatLab embeded Lcc C compiler.
//   Copyright (C) 2002-2020 ABINIT group (Mike Zhang)
//   This file is distributed under the terms of the
//   GNU General Public License, see ~abinit/COPYING
//   or http://www.gnu.org/copyleft/gpl.txt .
//   For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

#include "mex.h"
#include <stdio.h>
#include <string.h>

mxArray *data;
double *Data;
int dataDims;
int M, N, L; // Support max dimensions  is up to 3
char strSign[128] = "Eigenvalues (   eV  ) for nkpt=",
     strSign2[128] = "",
     strInnerInterval[128] = "coord)",
     strOuterInterval[128] = "";

int fstrstr(FILE *fp, const char* Str, const int ScanLen)
{ 
 int ScanChar;
 long curFP, bytes=0;
 
 if (ScanLen == 0)
   return 1;  // null string always return TRUE.
   
 curFP = ftell(fp);
 do {
   ScanChar = fgetc(fp);
   if ((char)ScanChar == Str[0])
   {
     int i;
     for (i=1; i<ScanLen; i++)
       if ((ScanChar = fgetc(fp)) != Str[i])
	 break;
     if (i == ScanLen) // Success!
       return 1;
     else 
       fseek(fp, -i, SEEK_CUR);
   }
 } while (ScanChar != EOF);
 // Search again from the beginning of the file
 fseek(fp, 0, SEEK_SET);
 do {
   ScanChar = fgetc(fp);
   if ((char)ScanChar == Str[0])
   {
     int i;
     for (i=1; i<ScanLen; i++)
       if ((ScanChar = fgetc(fp)) != Str[i])
	 break;
     if (i == ScanLen) // Success!
       return 1;
     else 
       fseek(fp, -i, SEEK_CUR);
   }
 } while (++bytes < curFP);
 // If search failed, the fp will return old position.
 return 0;
}

int readData(char *filename)
{
 FILE *fp;
 int row, col, z;
 int ScanLength = strlen(strSign);
 int ScanLength2 = strlen(strSign2);
 int InnerIntervalLength = strlen(strInnerInterval);
 int OuterIntervalLength = strlen(strOuterInterval);
  
 if (!(fp = fopen(filename, "rb")))
   return -1;

 // Now scan the buffer string:
 fseek(fp, 0L, SEEK_SET);
 if (!fstrstr(fp, strSign, ScanLength))
   return -3;
 if (!fstrstr(fp, strSign2, ScanLength2))
   return -7;
   
 for (z=0; z<L; z++)
 { 
   if (!fstrstr(fp, strOuterInterval, OuterIntervalLength)) 
     return -9;
   for (row=0; row<M; row++)
   {     
     if (!fstrstr(fp, strInnerInterval, InnerIntervalLength)) // Bypass "kpt#...coord)" parts
       return -5;
     for (col=0; col<N; col++)
       fscanf(fp, "%lf", &Data[row+M*col+z*M*N]);
   }   
 }

 fclose(fp);
 return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 char name[128];
 int dims[3]={1, 1, 1}, i, status;
 double *reDIMS; 
 
 if (nrhs <= 1)
 {
   mexPrintf("\nData Extractor\n"
             "Usage: mexData([M N L], ReadFileNameString, SignatureString, InnerIntervalString, OuterIntervalString, SubSignatureString);\n"
             "       mexData extracts data from ASCII code file which contains useful data as well as text infomation.\n\n"
             "       As for two-dimension data, mexData searches Signature_String which indicates the useful data will appear,\n"
             "       then searches the Inner_Interval_String to locate the accurate position of the useful data.\n"
             "       After the N columns data are read, mexData again searches Interval_String to prepare reading another row data.\n"
             "       When all M-row by N-column data are extracted, mexData returns a double float M-by-N matrix and finishes its work.\n\n"
             "       As for three-dimension M-by-N-by-L data, Outer_Interval_String is used to differentiate every L dimension.\n"
             "       mexData returns M-by-N-by-L matrix.\n\n"
             "       Modified at 7/30/2002, add OuterIntervalString, SubSignatureString.\n"
             "       Modified at 8/11/2002, correct L dimension disposition(use Data[row+M*col+z*M*N] instead of Data[L*row+M*L*col+z]).\n"
             "       Used by your own warranty, Mike.Zhang(mz24cn@hotmail.com), June 2002.\n\n");
   return;
 }
 reDIMS = mxGetPr(prhs[0]);
 dataDims = mxGetN(prhs[0]);
 for (i=0; i<dataDims; i++)
   dims[i] = (int)reDIMS[i];
 M = dims[0]; 
 N = dims[1];
 L = dims[2];
 data = mxCreateNumericArray(dataDims, dims, mxDOUBLE_CLASS, mxREAL);
 Data = mxGetPr(data);
 
 if (nrhs >= 2)
 {
   mxGetString(prhs[1], name, 127);
   if (nrhs >= 3)
     mxGetString(prhs[2], strSign, 127);
   if (nrhs >= 4)
     mxGetString(prhs[3], strInnerInterval, 127);
   if (nrhs >= 5)
     mxGetString(prhs[4], strOuterInterval, 127);
   if (nrhs >= 6)
     mxGetString(prhs[5], strSign2, 127);
 }
 mexPrintf("Data Extractor\n"
           "Read file              :      '%s'\n"
           "Signature string       :      '%s'\n"
           "InnerInterval string   :      '%s'\n"
           "OuterInterval string   :      '%s'\n"
           "2nd Signature string   :      '%s'\n"
           "Output array dimensions:      %d x %d x %d\n", 
           name, strSign, strInnerInterval, strOuterInterval, strSign2, M, N, L);
 status = readData(name);
 if (status == 1)
   mexPrintf("Data extracted success.\n");
 else if (status == -1)
   mexPrintf("Error: Can't open file '%s'.\n", name);
 else if (status == -3)
   mexPrintf("Error: Can't find string '%s'.\n", strSign);
 else if (status == -5)
   mexPrintf("Error: Can't find string '%s'.\n", strInnerInterval);
 else if (status == -7)
   mexPrintf("Error: Can't find string '%s'.\n", strSign2);
 else if (status == -9)
   mexPrintf("Error: Can't find string '%s'.\n", strOuterInterval);
   
 plhs[0] = data;
} 
