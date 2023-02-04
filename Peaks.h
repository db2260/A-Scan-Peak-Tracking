#pragma once

#ifndef _PEAKS_
#define _PEAKS_

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

namespace Peaks
{
        /**
         * Defines a point. X values are integers. Y values are floating point.
         */
        class GraphPoint
        {
        public:
                float x;
                float y;

                GraphPoint() { x = 0; y = (float)0.0; }
                GraphPoint(float newX, float newY) { x = newX; y = newY; }
                GraphPoint(const GraphPoint& rhs) { x = rhs.x; y = rhs.y; }

                GraphPoint& operator=(const GraphPoint& rhs)
                {
                        x = rhs.x;
                        y = rhs.y;
                        return *this;
                }

                bool roughlyEqual(float a, float b, float epsilon) const
                {
                        float absA = fabs(a);
                        float absB = fabs(b);
                        return fabs(a - b) <= ( (absA < absB ? absB : absA * epsilon) );
                }

                bool operator==(const GraphPoint& rhs) const
                {
                        bool xEqual = roughlyEqual(x, rhs.x, (float)0.0001);
                        bool yEqual = roughlyEqual(y, rhs.y, (float)0.0001);
                        return xEqual && yEqual;
                }

                void clear()
                {
                        x = (float)0.0;
                        y = (float)0.0;
                }



        };

        typedef std::vector<GraphPoint> GraphLine;

        /**
         * List of points.
         */

        /**
         * Defines a peak. A peak is described by three points: a left trough, a peak, and a right trough.
         */
        class GraphPeak
        {
        public:
                GraphPoint leftTrough;
                GraphPoint peak;
                GraphPoint rightTrough;
                float area;

                GraphPeak() { clear(); }

                GraphPeak(const GraphPeak& rhs)
                {
                        leftTrough = rhs.leftTrough;
                        peak = rhs.peak;
                        rightTrough = rhs.rightTrough;
                        area = rhs.area;
                }

                GraphPeak& operator=(const GraphPeak& rhs)
                {
                        leftTrough = rhs.leftTrough;
                        peak = rhs.peak;
                        rightTrough = rhs.rightTrough;
                        area = rhs.area;
                        return *this;
                }

                bool operator==(const GraphPeak& rhs) const
                {
                        return (leftTrough == rhs.leftTrough) && (peak == rhs.peak) && (rightTrough == rhs.rightTrough);
                }

                bool operator < (const GraphPeak& rhs) const { return (area < rhs.area); }
                bool operator > (const GraphPeak& rhs) const { return (area > rhs.area); }

                void clear()
                {
                        leftTrough.clear();
                        peak.clear();
                        rightTrough.clear();
                        area = (float)0.0;
                }
        };

        /**
         * List of peaks.
         */
        typedef std::vector<GraphPeak> GraphPeakList;

        /**
         * Collection of peak finding algorithms.
         */
        class Peaks
        {

        public:
         typedef std::vector<float> NumVec;
                 typedef GraphPoint GraphPoint;
                /**
                 * Returns a list of all statistically significant peaks in the given waveform.
                 * These are defined as peaks that rise more than one standard deviation above the mean for at least three points on the x axis.
                 */
                static GraphPeakList findPeaksOverThreshold(float* data, size_t dataLen, size_t* numPeaks, float threshold = 0.0);
                static GraphPeakList findPeaksOverStd(float* data, size_t dataLen, size_t* numPeaks, float sigmas = 1.0);
                static GraphPeakList findPeaksOverThreshold(const std::vector<float>& data, float threshold = 0.0);
                static GraphPeakList findPeaksOverStd(const std::vector<float>& data, float sigmas = 1.0);
                static GraphPeakList findPeaksOverThreshold(const GraphLine& data, float threshold = 0.0);
                static GraphPeakList findPeaksOverStd(const GraphLine& data, float sigmas = 1.0);
        std::vector<NumVec> readThreeAxisDataFromCsv(const std::string& fileName);
        std::vector<GraphPoint> Peaks::findPeaks(const std::vector<NumVec>& csvData, float threshold);
        std::vector<std::vector<GraphPoint>> Peaks::print(const std::vector<NumVec>& csvData);

        private:
                static float average(const float* data, size_t numPoints);
                static float average(const std::vector<float>& data);
                static float average(const GraphLine& data);

                static float variance(const float* data, size_t numPoints, float mean);
                static float variance(const std::vector<float>& data, float mean);
                static float variance(const GraphLine& data, float mean);

                static float standardDeviation(const float* data, size_t numPoints, float mean);
                static float standardDeviation(const std::vector<float>& data, float mean);
                static float standardDeviation(const GraphLine& data, float mean);

                static void computeArea(float* data, size_t dataLen, GraphPeak& currentPeak);
                static void computeArea(const std::vector<float>& data, GraphPeak& currentPeak);
                static void computeArea(const GraphLine& data, GraphPeak& currentPeak);
        };
};

#endif