@echo off
SET error=0
cd D:\Data\Jenkins\mRNADynamics-tests

REM ExportDataForFISH
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testExport_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testExport_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testExport_2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry.log" -r "runTestCase(testExport_2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testExport_2017_10_13_20170921_4.log" -r "runTestCase(testExport_2017_10_13_20170921_4);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testExport_2018_06_05_A140P_MSE_30uW_550V.log" -r "runTestCase(testExport_2018_06_05_A140P_MSE_30uW_550V);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testExportAndTifs_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testExportAndTifs_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1

REM SegmentSpots
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testSegmentSpots_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testSegmentSpots_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1

REM SegmentSpotsML
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testTifsCreation_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testTifsCreation_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testGenerateDogsWeka_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testGenerateDogsWeka_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\testSegmentSpotsML_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testSegmentSpotsML_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1

REM trackmRNADynamics
C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\testTrackmRNADynamics_2015_07_25_P2P_75uW_bi_short.log" -r "runTestCase(testTrackmRNADynamics_2015_07_25_P2P_75uW_bi_short);"
if %errorlevel% == 1 set error=1

REM C:\"Program Files"\MATLAB\R2017b\bin\matlab.exe -wait -nodesktop -nosplash -sd "D:\Data\Jenkins\mRNADynamics-tests" -logfile "D:\Data\Jenkins\logs\tests.log" -r "runTestCase();"
exit /b %error%