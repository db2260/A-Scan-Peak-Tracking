#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <memory>
#include <algorithm>
#include <map>
#include <iterator>
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <cassert>

//extern "C" {
	#include "../gnuplot-cpp/gnuplot_i.c"
//}

#include "Peaks.h"
#include "savgol.hpp"
#include "polyfit.cpp"
#include "gaussian_fit.cpp"

using namespace std;
//using namespace Peaks;

void fitting(double A2, double A1, double A0, gnuplot_ctrl *h, double* x) {
        double sigma=sqrt(-1/(2*A2));
        double mean=A1*sigma*sigma;
        double A=exp(A0+mean*mean/(2*sigma*sigma));

        double y[201];
        for(int i=0; i<201; i++){
                y[i] = (A*exp(-1*pow((x[i]-mean), 2))/(2*pow(sigma, 2)));
        }

        gnuplot_setstyle(h, "lines");
        gnuplot_plot_xy(h, x, y, 201, "fit");
}

namespace Peaks
{
        // Computes the area of the given peak.
        void Peaks::computeArea(float* data, size_t dataLen, GraphPeak& currentPeak)
        {
                currentPeak.area = (float)0.0;

                if (currentPeak.leftTrough.x < currentPeak.rightTrough.x)
                {
                        for (size_t index = (size_t)currentPeak.leftTrough.x + 1; index <= currentPeak.rightTrough.x; ++index)
                        {
                                float b = data[index] + data[index - 1];
                                currentPeak.area += ((float)0.5 * b);
                        }
                }
        }

        // Returns a list of peaks in the given array of numeric values. Only peaks that go above the given threshold will be counted.
        GraphPeakList Peaks::findPeaksOverThreshold(float* data, size_t dataLen, size_t* numPeaks, float threshold)
        {
                std::vector<GraphPeak> peaks;
                GraphPeak currentPeak;

                for (size_t x = 0; x < dataLen; ++x)
                {
                        float y = data[x];

                        if (y < threshold)
                        {
                                // Have we found a peak? If so, add it and start looking for the next one.
                                if (currentPeak.rightTrough.x > 0)
                                {
                                        // Still descending
                                        if (y <= currentPeak.rightTrough.y)
                                        {
                                                currentPeak.rightTrough.x = x;
                                                currentPeak.rightTrough.y = y;
                                        }

                                        // Rising
                                        else
                                        {
                                                Peaks::computeArea(data, dataLen, currentPeak);
                                                peaks.push_back(currentPeak);
                                                currentPeak.clear();
                                        }
                                }

                                // Are we looking for a left trough?
                                else if (currentPeak.leftTrough.x == 0)
                                {
                                        currentPeak.leftTrough.x = x;
                                        currentPeak.leftTrough.y = y;
                                }

                                // If we have a left trough and an existing peak, assume this is the right trough - for now.
                                else if ((currentPeak.peak.x > currentPeak.leftTrough.x) && (currentPeak.leftTrough.x > 0))
                                {
                                        currentPeak.rightTrough.x = x;
                                        currentPeak.rightTrough.y = y;
                                }
                                else
                                {
                                        currentPeak.leftTrough.x = x;
                                        currentPeak.leftTrough.y = y;
                                }
                        }
                        else if (currentPeak.leftTrough.x > 0) // Left trough is set.
                        {
                                // Are we looking for a peak or is this bigger than the current peak?
                                if (currentPeak.peak.x == 0 || y >= currentPeak.peak.y)
                                {
                                        currentPeak.peak.x = x;
                                        currentPeak.peak.y = y;
                                }
                        }
                        else if (currentPeak.rightTrough.x > 0) // Right trough is set.
                        {
                                Peaks::computeArea(data, dataLen, currentPeak);
                                peaks.push_back(currentPeak);
                                currentPeak.clear();
                        }
                        else // Nothing is set, but the value is above the threshold.
                        {
                                currentPeak.leftTrough.x = x;
                                currentPeak.leftTrough.y = y;
                        }
                }

                return peaks;
        }

        // Returns a list of peaks in the given array of numeric values. Only peaks that go above the given sigma line will be counted.
        GraphPeakList Peaks::findPeaksOverStd(float* data, size_t dataLen, size_t* numPeaks, float sigmas)
        {
                float mean = Peaks::average(data, dataLen);
                float stddev = sigmas * Peaks::standardDeviation(data, dataLen, mean);
                float threshold = mean + stddev;
                return Peaks::findPeaksOverThreshold(data, dataLen, numPeaks, threshold);
        }

        // Computes the area of the given peak.
        void Peaks::computeArea(const std::vector<float>& data, GraphPeak& currentPeak)
        {
                currentPeak.area = (float)0.0;

                if (currentPeak.leftTrough.x < currentPeak.rightTrough.x)
                {
                        for (size_t index = (size_t)currentPeak.leftTrough.x + 1; index <= currentPeak.rightTrough.x; ++index)
                        {
                                float b = data.at(index) + data.at(index - 1);
                                currentPeak.area += ((float)0.5 * b);
                        }
                }
        }

        // Returns a list of peaks in the given array of numeric values. Only peaks that go above the given threshold will be counted.
        GraphPeakList Peaks::findPeaksOverThreshold(const std::vector<float>& data, float threshold)
        {
                std::vector<GraphPeak> peaks;
                GraphPeak currentPeak;

                float x = 0;

                for (auto iter = data.begin(); iter < data.end(); ++iter, ++x)
                {
                        float y = *iter;

                        if (y < threshold)
                        {
                                // Have we found a peak? If so, add it and start looking for the next one.
                                if (currentPeak.rightTrough.x > 0)
                                {
                                        // Still descending
                                        if (y <= currentPeak.rightTrough.y)
                                        {
                                                currentPeak.rightTrough.x = x;
                                                currentPeak.rightTrough.y = y;
                                        }

                                        // Rising
                                        else
                                        {
                                                Peaks::computeArea(data, currentPeak);
                                                peaks.push_back(currentPeak);
                                                currentPeak.clear();
                                        }
                                }

                                // Are we looking for a left trough?
                                else if (currentPeak.leftTrough.x == 0)
                                {
                                        currentPeak.leftTrough.x = x;
                                        currentPeak.leftTrough.y = y;
                                }

                                // If we have a left trough and an existing peak, assume this is the right trough - for now.
                                else if ((currentPeak.peak.x > currentPeak.leftTrough.x) && (currentPeak.leftTrough.x > 0))
                                {
                                        currentPeak.rightTrough.x = x;
                                        currentPeak.rightTrough.y = y;
                                }
                                else
                                {
                                        currentPeak.leftTrough.x = x;
                                        currentPeak.leftTrough.y = y;
                                }
                        }
                        else if (currentPeak.leftTrough.x > 0) // Left trough is set.
                        {
                                // Are we looking for a peak or is this bigger than the current peak?
                                if (currentPeak.peak.x == 0 || y >= currentPeak.peak.y)
                                {
                                        currentPeak.peak.x = x;
                                        currentPeak.peak.y = y;
                                }
                        }
                        else if (currentPeak.rightTrough.x > 0) // Right trough is set.
                        {
                                Peaks::computeArea(data, currentPeak);
                                peaks.push_back(currentPeak);
                                currentPeak.clear();
                        }
                        else // Nothing is set, but the value is above the threshold.
                        {
                                currentPeak.leftTrough.x = x;
                                currentPeak.leftTrough.y = y;
                        }
                }

                return peaks;
        }

        // Returns a list of peaks in the given array of numeric values. Only peaks that go above the given sigma line will be counted.
        GraphPeakList Peaks::findPeaksOverStd(const std::vector<float>& data, float sigmas)
        {
                float mean = Peaks::average(data);
                float stddev = sigmas * Peaks::standardDeviation(data, mean);
                float threshold = mean + stddev;
                return Peaks::findPeaksOverThreshold(data, threshold);
        }

        float Peaks::average(const float* data, size_t numPoints)
        {
                float sum = 0;

                for (auto index = 0; index < numPoints; index++)
                        sum = sum + data[index];
                return sum / (float)numPoints;
        }

        float Peaks::average(const std::vector<float>& data)
        {
                float sum = 0;

                for (auto iter = data.begin(); iter != data.end(); ++iter)
                        sum = sum + (*iter);
                return sum / (float)data.size();
        }

        float Peaks::average(const GraphLine& data)
        {
                float sum = (float)0.0;

                for (auto iter = data.begin(); iter != data.end(); ++iter)
                        sum = sum + (*iter).y;
                return sum / (float)data.size();
        }

        float Peaks::variance(const float* data, size_t numPoints, float mean)
        {
                float numerator = (float)0.0;

                for (auto index = 0; index < numPoints; index++)
                        numerator = numerator + ((data[index] - mean) * (data[index] - mean));
                return numerator / (float)(numPoints - 1);
        }

        float Peaks::variance(const std::vector<float>& data, float mean)
        {
                float numerator = (float)0.0;

                for (auto iter = data.begin(); iter != data.end(); ++iter)
                        numerator = numerator + ((*iter - mean) * (*iter - mean));
                return numerator / (float)(data.size() - 1);
        }

        float Peaks::variance(const GraphLine& data, float mean)
        {
                float numerator = (float)0.0;

                for (auto iter = data.begin(); iter != data.end(); ++iter)
                        numerator = numerator + (((*iter).y - mean) * ((*iter).y - mean));
                return numerator / (float)(data.size() - 1);
        }

        float Peaks::standardDeviation(const float* data, size_t numPoints, float mean)
        {
                float var = variance(data, numPoints, mean);
                return sqrt(var);
        }

        float Peaks::standardDeviation(const std::vector<float>& data, float mean)
        {
                float var = variance(data, mean);
                return sqrt(var);
        }

        float Peaks::standardDeviation(const GraphLine& data, float mean)
        {
                float var = variance(data, mean);
                return sqrt(var);
        }

        // Computes the area of the given peak.
        void Peaks::computeArea(const GraphLine& data, GraphPeak& currentPeak)
        {
                currentPeak.area = (float)0.0;

                if (currentPeak.leftTrough.x < currentPeak.rightTrough.x)
                {
                        auto pointIter = std::find(data.begin(), data.end(), currentPeak.leftTrough);

                        if (pointIter != data.end())
                        {
                                GraphPoint prevPoint = (*pointIter);
                                ++pointIter;

                                for (; pointIter != data.end(); ++pointIter)
                                {
                                        float b = (*pointIter).y + prevPoint.y;
                                        currentPeak.area += ((float)0.5 * b);
                                        prevPoint = (*pointIter);
                                }
                        }
                }
        }

        // Returns a list of peaks in the given array of graph points. Only peaks that go above the given threshold will be counted.
        GraphPeakList Peaks::findPeaksOverThreshold(const GraphLine& data, float threshold)
        {
                std::vector<GraphPeak> peaks;
                GraphPeak currentPeak;

                for (auto iter = data.begin(); iter < data.end(); ++iter)
                {
                        const GraphPoint& pt = *iter;

                        if (pt.y < threshold)
                        {
                                // Have we found a peak? If so, add it and start looking for the next one.
                                if (currentPeak.rightTrough.x > 0) // Right trough is set
                                {
                                        // Still descending
                                        if (pt.y <= currentPeak.rightTrough.y)
                                        {
                                                currentPeak.rightTrough.x = pt.x;
                                                currentPeak.rightTrough.y = pt.y;
                                        }

                                        // Rising
                                        else
                                        {
                                                Peaks::computeArea(data, currentPeak);
                                                peaks.push_back(currentPeak);
                                                currentPeak.clear();
                                        }
                                }

                                // Are we looking for a left trough?
                                else if (currentPeak.leftTrough.x == 0) // Left trough is not set.
                                {
                                        currentPeak.leftTrough = pt;
                                }

                                // If we have a left trough and an existing peak, assume this is the right trough - for now.
                                else if ((currentPeak.peak.x > currentPeak.leftTrough.x) && (currentPeak.leftTrough.x > 0))
                                {
                                        currentPeak.rightTrough = pt;
                                }
                                else
                                {
                                        currentPeak.leftTrough = pt;
                                }
                        }
                        else if (currentPeak.leftTrough.x > 0) // Left trough is set.
                        {
                                // Are we looking for a peak or is this bigger than the current peak, making it the real peak?
                                if (currentPeak.peak.x == 0 || pt.y >= currentPeak.peak.y)
                                {
                                        currentPeak.peak = pt;
                                }
                        }
                        else if (currentPeak.rightTrough.x > 0) // Right trough is set.
                        {
                                Peaks::computeArea(data, currentPeak);
                                peaks.push_back(currentPeak);
                                currentPeak.clear();
                        }
                        else // Nothing is set, but the value is above the threshold.
                        {
                                currentPeak.leftTrough = pt;
                        }
                }

                return peaks;
        }

        // Returns a list of peaks in the given array of graph points. Only peaks that go above the given sigma line will be counted.
        GraphPeakList Peaks::findPeaksOverStd(const GraphLine& data, float sigmas)
        {
                float mean = average(data);
                float stddev = sigmas * standardDeviation(data, mean);
                float threshold = mean + stddev;
                return Peaks::findPeaksOverThreshold(data, threshold);
        }


         typedef std::vector<float> NumVec;
        // Loads the peak data test file. Expected format is timestamp, x, y, z.
        std::vector<NumVec> Peaks::readThreeAxisDataFromCsv(const std::string& fileName)
        {
                std::vector<NumVec> columns;
                std::vector<float> xList, yList, zList;

                std::string line;
                std::ifstream infile(fileName);

        std::getline(infile, line, '\n');
        std::getline(infile, line, '\n');
        std::getline(infile, line, '\n');

                while (std::getline(infile, line, '\n'))
                {
                        std::istringstream iss(line);
                        std::string value;
                        //uint64_t ts = 0;
                        float x = 0.0, y = 0.0, z = 0.0;

                        //iss >> ts;
                        //iss.ignore(256, ',');
                        iss >> x;
                        iss.ignore(256, ',');
                        iss >> y;
                        iss.ignore(256, ',');
                        iss >> z;

                        //tsList.push_back(ts);
                        xList.push_back(x);
                        yList.push_back(y);
                        zList.push_back(z);
                }

                //columns.push_back(tsList);
                columns.push_back(xList);
                columns.push_back(yList);
                columns.push_back(zList);
                return columns;
        }


        std::vector<GraphPoint> Peaks::findPeaks(const std::vector<NumVec>& csvData, float threshold)
        {

                auto csvIter = csvData.begin();

                GraphPeakList peaks;

                for (; csvIter != csvData.end(); ++csvIter)
                {
                        peaks = Peaks::findPeaksOverThreshold((*csvIter), threshold);
                        //peaks2 = Peaks::Sigmas((*csvIter), sigmas);
                }


                std::vector<Peaks::GraphPoint> vec;
                GraphPeakList::iterator i = peaks.begin();
                for(; i != peaks.end(); ++i) {
                        if((*i).peak.y >= 1.75 && (*i).peak.x > 124504 && (*i).peak.x < 175215) {
                                vec.push_back((*i).peak);
                        }
                }

                for(int i=0; i<vec.size(); i++) {
                        vec.at(i).x = vec.at(i).x/500;
                }

                return vec;
        }


                std::vector<GraphLine> Peaks::print(const std::vector<NumVec>& csvData)
        {

        std::vector<GraphPoint> curve1;
        std::vector<GraphPoint> curve2;
        std::vector<GraphPoint> curve3;

        int a = 0;
        cout << "size " << csvData.at(0).size() << endl;
        for(double i=21103; i<=21303; i++) {
            GraphPoint v;
            v.x = csvData.at(0).at(i);
            v.y = csvData.at(1).at(i);
//            cout << "before1\n" << endl;
    //        cout << v.x << " " << v.y << endl;
            curve1.push_back(v);
  //          cout << "after1\n" << endl;
          //  cout << static_cast<double>(curve1.at(a).x) << " " << static_cast<double>(curve1.at(a).y) << endl;
            a++;
        }
        printf("%ld\n\n\n", a);
        a=0;
        for(double j=124436; j<=124636; j++) {
            GraphPoint v;
            v.x = csvData.at(0).at(j);
            v.y = csvData.at(1).at(j);
            curve2.push_back(v);
            //cout << curve2.at(b+124436).x << " " << curve1.at(b+124436).y << endl;
           // b++;
           //cout << static_cast<double>(curve2.at(a).x) << " " << static_cast<double>(curve2.at(a).y) << endl;
            a++;
        }
        printf("%ld\n\n\n", a);
        a=0;
        for(double k=163963; k<=164163; k++) {
            GraphPoint v;
            v.x = csvData.at(0).at(k);
            v.y = csvData.at(1).at(k);
            curve3.push_back(v);
            //cout << curve3.at(c+163963).x << " " << curve1.at(c+163963).y << endl;
           // c++;
           //cout << static_cast<double>(curve3.at(a).x) << " " << static_cast<double>(curve3.at(a).y) << endl;
            a++;
        }
        printf("%ld\n\n\n", a);
        a=0;

        std::vector<std::vector<GraphPoint>> curves;

        curves.push_back(curve1);
        curves.push_back(curve2);
        curves.push_back(curve3);

        return curves;
        }
}

int main(int argc, char** argv) {

        gnuplot_ctrl * h;
        h = gnuplot_init();

    char file_path[] = "../deepta/ascans/1 (2).csv";
    #include "Peaks.h"
        #include "savgol.hpp"
        Peaks::Peaks inst;

    std::ifstream read(file_path);
    auto csvData = inst.readThreeAxisDataFromCsv(file_path);

    int s = csvData.at(1).size();
    double a[s], b[s];
    for(int i=0; i<s; i++) {
        a[i] = csvData.at(0).at(i);
        b[i] = csvData.at(1).at(i);
    }
        //float t = 1.0;
    std::vector<Peaks::GraphPoint> axisPeaks =inst.findPeaks(csvData, 1.0);


        int s2 = axisPeaks.size();
        double vecx[s2], vecy[s2];
        for(int j=0; j<s2; j++) {
                vecx[j] = axisPeaks.at(j).x;
                vecy[j] = axisPeaks.at(j).y;
                std::cout << "x: " << vecx[j] << " y: " << vecy[j] << std::endl;
        }

        const int windowsize = 25;
        std::vector<float> ydash1, y2;
        filter::savgol(csvData.at(2).begin(), csvData.at(2).end(), std::back_inserter(ydash1), filter::SmoothQuarticQuintic(windowsize));
        filter::savgol(ydash1.begin(), ydash1.end(), std::back_inserter(y2), filter::SmoothQuarticQuintic(windowsize));


    std::vector<Peaks::GraphLine> curves = inst.print(csvData);
        //cout << "after print\n" << endl;

        //GUASS CURVE 1
        double curvex[201], curvey[201];
        for(int i=0; i<=200; i++) {
                curvex[i] = (static_cast<double>(curves.at(0).at(i).x));
                curvey[i] = (static_cast<double>(curves.at(0).at(i).y));
        }

        std::vector<double> curvexv, curveyv;
        for(int i=0; i<=200; i++) {
                curvexv.push_back(curvex[i]);
                curveyv.push_back(curvey[i]);
        }
     
       // curve_fit(curvexv, curveyv, h,1);
        double A2=-0.0159539;
        double A1=1.3531;
        double A0=-28.6898;
      //  fitting(A2, A1, A0, h, curvex);
        

        //GUASS CURVE 2
        double curvex2[201], curvey2[201];
        for(int i=0; i<=200; i++) {
                curvex2[i] = (static_cast<double>(curves.at(1).at(i).x));
                curvey2[i] = (static_cast<double>(curves.at(1).at(i).y));
                //cout << curvex.at(i) << " " << curvey.at(i) << endl;
        }
        std::vector<double> curvex2v, curvey2v;
        for(int i=0; i<=200; i++) {
                curvex2v.push_back(curvex2[i]);
                curvey2v.push_back(curvey2[i]);
        }
        
        curve_fit(curvex2v, curvey2v, h,2);
        A2=-1.73322;
        A1=860.321;
        A0=-107122;
        //fitting(A2, A1, A0, h, curvex2);
        


        //GUASS CURVE 3
        double curvex3[201], curvey3[201];
        for(int i=0; i<=200; i++) {
                curvex3[i] = (static_cast<double>(curves.at(2).at(i).x));
                curvey3[i] = (static_cast<double>(curves.at(2).at(i).y));
                //cout << curvex.at(i) << " " << curvey.at(i) << endl;
        }
        std::vector<double> curvex3v, curvey3v;
        for(int i=0; i<=200; i++) {
                curvex3v.push_back(curvex3[i]);
                curvey3v.push_back(curvey3[i]);
        }
     
        curve_fit(curvex3v, curvey3v, h,3);
        A2=-0.0152269;
        A1=9.88596;
        A0=-1671.58;
        //fitting(A2, A1, A0, h, curvex3);
         


   gnuplot_setstyle(h, "lines");
   gnuplot_plot_xy(h, a, b, s, "data");
 //   gnuplot_resetplot(h);
   gnuplot_setstyle(h, "points");
   gnuplot_plot_xy(h, vecx, vecy, s2, "peak");

       sleep(200);
      gnuplot_close(h);


        return 0;
}