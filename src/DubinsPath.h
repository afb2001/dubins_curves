#ifndef SRC_DUBINSPATH_H
#define SRC_DUBINSPATH_H

#include <vector>

namespace dubins {

    class DubinsPath {
    public:
        DubinsPath() = default;

        /**
         * Find the shortest Dubins path between q1 = {x1, y1, theta1} and q2 = {x2, y2, theta2} with the given
         * turning radius.
         * @param x1
         * @param y1
         * @param theta1
         * @param x2
         * @param y2
         * @param theta2
         * @param radius
         */
        DubinsPath(double x1, double y1, double theta1, double x2, double y2, double theta2, double radius);

        bool isValid() const { return m_Radius > 0; } // basically for documentation

        void sample(double& x, double& y, double& theta, double& distance) const;

    private:
        enum PathType {
            LRL,
            LSL,
            LSR,
            RLR,
            RSL,
            RSR,
            INVALID
        };

        enum SegmentType {
            L,
            R,
            S,
        };

        struct DubinsParameters {
            double length1, length2, length3;
            PathType pathType;

            double totalLength() const { return length1 + length2 + length3; }
        };

        struct LineSegment {
            double x1{}, y1{}, x2{}, y2{};
            bool valid = false;
            double length() const { return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)); }
        };

        double m_X1, m_Y1, m_Theta1, m_X2, m_Y2, m_Theta2, m_Radius = -1;
        DubinsParameters m_DubinsParameters{0, 0, 0, INVALID};

        // Centers of circles to which q1 is tangent. C1 starts with a left turn, C2 starts with a right.
        double m_Center1X, m_Center1Y, m_Center2X, m_Center2Y;
        // Centers of circles to which q2 is tangent. C3 ends with a left turn, C4 ends with a right.
        double m_Center3X, m_Center3Y, m_Center4X, m_Center4Y;

        DubinsParameters computeLRL() const;

        DubinsParameters computeLSL() const;

        DubinsParameters computeLSR() const;

        DubinsParameters computeRLR() const;

        DubinsParameters computeRSL() const;

        DubinsParameters computeRSR() const;

        inline void potentiallyReplaceParameters(const DubinsParameters& parameters);

        static inline SegmentType getSegmentType(int index, PathType pathType);

        static double arcLength(double x1, double y1, double x2, double y2, double radius, bool directionIsLeft);

        static std::vector<LineSegment> computePossibleTangents(double x1, double y1, double x2, double y2,
                                                         double radius);

        /**
         * Floating point modulus suitable for rings
         *
         * fmod doesn't behave correctly for angular quantities, this function does
         */
        static double fmodr(double x, double y)
        {
            return x - y*floor(x/y);
        }

        static double mod2pi(double theta)
        {
            return fmodr(theta, 2 * M_PI);
        }

    };

}

#endif //SRC_DUBINSPATH_H
