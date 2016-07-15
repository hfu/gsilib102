/*------------------------------------------------------------------------------
* gsiqc.c : quality check header
*
*    Copyright (C) 2014 by Geospatial Information Authority of Japan,
*    All rights reserved.
*    
*    Released under the BSD, and GPL Licenses.
*
*
*  Original software: RTKLIB ver.2.4.2 p4
*
*    Copyright (C) 2007-2013 by T.Takasu, All rights reserved.
*
*
* references :
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#ifndef GSIQC_H
#define GSIQC_H

#define OUTFILE_S        0    /* output file number for S-file */
#define OUTFILE_AZ       1    /* output file number for azimuth */
#define OUTFILE_EL       2    /* output file number for elevation */
#define OUTFILE_ION      3    /* output file number for ion */
#define OUTFILE_IOD      4    /* output file number for iod */
#define OUTFILE_ION5     5    /* output file number for ion5 */
#define OUTFILE_IOD5     6    /* output file number for iod5 */
#define OUTFILE_MP12     7    /* output file number for mp12 */
#define OUTFILE_MP21     8    /* output file number for mp21 */
#define OUTFILE_MP15     9    /* output file number for mp15 */
#define OUTFILE_MP51    10    /* output file number for mp51 */
#define OUTFILE_MP25    11    /* output file number for mp25 */
#define OUTFILE_MP52    12    /* output file number for mp52 */
#define OUTFILE_S1      13    /* output file number for s1 */
#define OUTFILE_S2      14    /* output file number for s2 */

#define OUTFILE_MAX    13    /* output file total number*/

#define ASCII_SYMBOL1      "C"
#define ASCII_SYMBOL2      "m"
#define ASCII_SYMBOL3      "I"
#define ASCII_SYMBOL4      "M"
#define ASCII_SYMBOL5      "1"
#define ASCII_SYMBOL6      "2"
#define ASCII_SYMBOL7      "Z"
#define ASCII_SYMBOL8      "5"
#define ASCII_SYMBOL9      "-"
#define ASCII_SYMBOL10     "L"
#define ASCII_SYMBOL11     "?"
#define ASCII_SYMBOL12     "+"
#define ASCII_SYMBOL13     "^"
#define ASCII_SYMBOL14     "."
#define ASCII_SYMBOL15     "c"
#define ASCII_SYMBOL16     ":"
#define ASCII_SYMBOL17     "="
#define ASCII_SYMBOL18     "z"
#define ASCII_SYMBOL19     "~"
#define ASCII_SYMBOL20     "*"
#define ASCII_SYMBOL21     ","
#define ASCII_SYMBOL22     "a"
#define ASCII_SYMBOL23     ";"
#define ASCII_SYMBOL24     "e"
#define ASCII_SYMBOL25     "s"
#define ASCII_SYMBOL26     "o"
#define ASCII_SYMBOL27     "y"
#define ASCII_SYMBOL28     "N"
#define ASCII_SYMBOL29     "_"

#define ASCII_SYMBOL1_NO      -13
#define ASCII_SYMBOL2_NO      -12
#define ASCII_SYMBOL3_NO      -4
#define ASCII_SYMBOL4_NO      -10
#define ASCII_SYMBOL5_NO      -8
#define ASCII_SYMBOL6_NO      -7
#define ASCII_SYMBOL7_NO      -6
#define ASCII_SYMBOL8_NO      -5
#define ASCII_SYMBOL9_NO      -3
#define ASCII_SYMBOL10_NO     -11
#define ASCII_SYMBOL11_NO     16
#define ASCII_SYMBOL12_NO     -9
#define ASCII_SYMBOL13_NO     -1
#define ASCII_SYMBOL14_NO     14
#define ASCII_SYMBOL15_NO     10
#define ASCII_SYMBOL16_NO     12
#define ASCII_SYMBOL17_NO      1
#define ASCII_SYMBOL18_NO      5
#define ASCII_SYMBOL19_NO      9
#define ASCII_SYMBOL20_NO     11
#define ASCII_SYMBOL21_NO     13
#define ASCII_SYMBOL22_NO      6
#define ASCII_SYMBOL23_NO      7
#define ASCII_SYMBOL24_NO      2
#define ASCII_SYMBOL25_NO      8
#define ASCII_SYMBOL26_NO      4
#define ASCII_SYMBOL27_NO      3
#define ASCII_SYMBOL28_NO     15
#define ASCII_SYMBOL29_NO     -2

#endif /* GSIQC_H */
