#include <cmath>
#include <cfloat>
#include <iostream>
#include "DubinsPath.h"

using namespace dubins;

DubinsPath::DubinsPath(double x1, double y1, double theta1, double x2, double y2, double theta2, double radius) :
    m_X1(x1), m_Y1(y1), m_Theta1(theta1), m_X2(x2), m_Y2(y2), m_Theta2(theta2), m_Radius(radius)
{
    // compute tangent circles
    auto perpTheta1 = mod2pi(m_Theta1 + M_PI_2);
    auto perpDx1 = radius * cos(perpTheta1);
    auto perpDy1 = radius * sin(perpTheta1);
    m_Center1X = m_X1 + perpDx1;
    m_Center1Y = m_Y1 + perpDy1;
    m_Center2X = m_X1 - perpDx1;
    m_Center2Y = m_Y1 - perpDy1;

    auto perpTheta2 = mod2pi(m_Theta2 + M_PI_2);
    auto perpDx2 = radius * cos(perpTheta2);
    auto perpDy2 = radius * sin(perpTheta2);
    m_Center3X = m_X2 + perpDx2;
    m_Center3Y = m_Y2 + perpDy2;
    m_Center4X = m_X2 - perpDx2;
    m_Center4Y = m_Y2 - perpDy2;

    // compute and store parameters
    // break ties the same way every time (in this order)
    potentiallyReplaceParameters(computeLRL());
    potentiallyReplaceParameters(computeLSL());
    potentiallyReplaceParameters(computeLSR());
    potentiallyReplaceParameters(computeRLR());
    potentiallyReplaceParameters(computeRSL());
    potentiallyReplaceParameters(computeRSR());
}

std::vector<DubinsPath::LineSegment> DubinsPath::computePossibleTangents(double x1, double y1, double x2, double y2,
                                                             double radius) {
    // Heavily influenced by https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Tangents_between_two_circles
    std::vector<LineSegment> result(4);
    double dSq = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    if (dSq <= 0) return result;

    double d = sqrt(dSq);
    double vx = (x2 - x1) / d;
    double vy = (y2 - y1) / d;

    int i = 0;

    // Let A, B be the centers, and C, D be points at which the tangent
    // touches first and second circle, and n be the normal vector to it.
    //
    // We have the system:
    //   n * n = 1          (n is a unit vector)
    //   C = A + r1 * n
    //   D = B +/- r2 * n
    //   n * CD = 0         (common orthogonality)
    //
    // n * CD = n * (AB +/- r2*n - r1*n) = AB*n - (r1 -/+ r2) = 0,  <=>
    // AB * n = (r1 -/+ r2), <=>
    // v * n = (r1 -/+ r2) / d,  where v = AB/|AB| = AB/d
    // This is a linear equation in unknown vector n.

    for (int sign1 = +1; sign1 >= -1; sign1 -= 2) {
        double c = (radius - sign1 * radius) / d;

        // Now we're just intersecting a line with a circle: v*n=c, n*n=1

        if (c*c > 1.0) continue;
        double h = sqrt(fmax(0.0, 1.0 - c*c));

        for (int sign2 = +1; sign2 >= -1; sign2 -= 2) {
            double nx = vx * c - sign2 * h * vy;
            double ny = vy * c + sign2 * h * vx;

            auto& a = result[i++];
            a.x1 = x1 + radius * nx;
            a.y1 = y1 + radius * ny;
            a.x2 = x2 + sign1 * radius * nx;
            a.y2 = y2 + sign1 * radius * ny;
            a.valid = true;
        }
    }

//    std::cerr << "Found " << i << " tangents" << std::endl;

    return result;
}

double DubinsPath::arcLength(double x1, double y1, double x2, double y2, double radius, bool directionIsLeft) {
    auto theta = atan2(y2, x2) - atan2(y1, x1);
    if (theta > 0 && directionIsLeft) theta += 2 * M_PI;
    else if (theta > 0 && !directionIsLeft) theta -= 2 * M_PI;
    return fabs(theta * radius);
}

DubinsPath::DubinsParameters DubinsPath::computeLRL() const {
    DubinsParameters result {0, 0, 0, INVALID};
    auto d = sqrt((m_Center3X - m_Center1X) * (m_Center3X - m_Center1X) + (m_Center3Y - m_Center1Y) * (m_Center3Y - m_Center1Y));
    if (d >= 4 * m_Radius) return result; // cannot make CCC path
    auto theta = atan2(m_Center3Y - m_Center1Y, m_Center3X - m_Center1X) - acos(d / (4 * m_Radius));
    auto x3 = m_Center1X + 2 * m_Radius * cos(theta);
    auto y3 = m_Center1Y + 2 * m_Radius * sin(theta);
    auto pt1X = m_Center1X + (x3 - m_Center1X) / 2;
    auto pt1Y = m_Center1Y + (y3 - m_Center1Y) / 2;
    auto pt2X = m_Center3X + (x3 - m_Center1X) / 2;
    auto pt2Y = m_Center3Y + (y3 - m_Center1Y) / 2;
    result.length1 = arcLength(m_X1, m_Y1, pt1X, pt1Y, m_Radius, true);
    result.length2 = arcLength(pt1X, pt1Y, pt2X, pt2Y, m_Radius, false);
    result.length3 = arcLength(pt2X, pt2Y, m_X2, m_Y2, m_Radius, true);
    result.pathType = LRL;
    return result;
}

DubinsPath::DubinsParameters DubinsPath::computeLSL() const {
    auto tangents = std::move(computePossibleTangents(m_Center1X, m_Center1Y, m_Center3X, m_Center3Y, m_Radius));
    auto minLength = DBL_MAX;
    DubinsParameters result {0, 0, 0, INVALID};
    for (const auto& s : tangents) if (s.valid) {
        DubinsParameters parameters {
            arcLength(m_Center1X, m_Center1Y, s.x1, s.y1, m_Radius, true),
            s.length(),
            arcLength(s.x2, s.y2, m_Center3X, m_Center3Y, m_Radius, true),
            LSL
        };
        auto length = parameters.totalLength();
        if (length < minLength) {
            minLength = length;
            result = parameters;
        }
    }
    return result;
}

DubinsPath::DubinsParameters DubinsPath::computeLSR() const {
    auto tangents = std::move(computePossibleTangents(m_Center1X, m_Center1Y, m_Center4X, m_Center4Y, m_Radius));
    auto minLength = DBL_MAX;
    DubinsParameters result {0, 0, 0, INVALID};
    for (const auto& s : tangents) if (s.valid) {
            DubinsParameters parameters {
                    arcLength(m_Center1X, m_Center1Y, s.x1, s.y1, m_Radius, true),
                    s.length(),
                    arcLength(s.x2, s.y2, m_Center4X, m_Center4Y, m_Radius, false),
                    LSR
            };
            auto length = parameters.totalLength();
            if (length < minLength) {
                minLength = length;
                result = parameters;
            }
        }
    return result;
}

DubinsPath::DubinsParameters DubinsPath::computeRLR() const {
    DubinsParameters result {0, 0, 0, INVALID};
    auto d = sqrt((m_Center4X - m_Center2X) * (m_Center4X - m_Center2X) + (m_Center4Y - m_Center2Y) * (m_Center4Y - m_Center2Y));
    if (d >= 4 * m_Radius) return result; // cannot make CCC path
    auto theta = atan2(m_Center4Y - m_Center2Y, m_Center4X - m_Center2X) - acos(d / (4 * m_Radius));
    auto x3 = m_Center2X + 2 * m_Radius * cos(theta);
    auto y3 = m_Center2Y + 2 * m_Radius * sin(theta);
    auto pt1X = m_Center2X + (x3 - m_Center2X) / 2;
    auto pt1Y = m_Center2Y + (y3 - m_Center2Y) / 2;
    auto pt2X = m_Center4X + (x3 - m_Center2X) / 2;
    auto pt2Y = m_Center4Y + (y3 - m_Center2Y) / 2;
    result.length1 = arcLength(m_X1, m_Y1, pt1X, pt1Y, m_Radius, false);
    result.length2 = arcLength(pt1X, pt1Y, pt2X, pt2Y, m_Radius, true);
    result.length3 = arcLength(pt2X, pt2Y, m_X2, m_Y2, m_Radius, false);
    result.pathType = RLR;
    return result;
}

DubinsPath::DubinsParameters DubinsPath::computeRSL() const {
    auto tangents = std::move(computePossibleTangents(m_Center2X, m_Center2Y, m_Center3X, m_Center3Y, m_Radius));
    auto minLength = DBL_MAX;
    DubinsParameters result {0, 0, 0, INVALID};
    for (const auto& s : tangents) if (s.valid) {
            DubinsParameters parameters {
                    arcLength(m_Center2X, m_Center2Y, s.x1, s.y1, m_Radius, false),
                    s.length(),
                    arcLength(s.x2, s.y2, m_Center3X, m_Center3Y, m_Radius, true),
                    RSL
            };
            auto length = parameters.totalLength();
            if (length < minLength) {
                minLength = length;
                result = parameters;
            }
        }
    return result;
}

DubinsPath::DubinsParameters DubinsPath::computeRSR() const {
    auto tangents = std::move(computePossibleTangents(m_Center2X, m_Center2Y, m_Center4X, m_Center4Y, m_Radius));
    auto minLength = DBL_MAX;
    DubinsParameters result {0, 0, 0, INVALID};
    for (const auto& s : tangents) if (s.valid) {
            DubinsParameters parameters {
                    arcLength(m_Center2X, m_Center2Y, s.x1, s.y1, m_Radius, false),
                    s.length(),
                    arcLength(s.x2, s.y2, m_Center4X, m_Center4Y, m_Radius, false),
                    RSR
            };
            auto length = parameters.totalLength();
            if (length < minLength) {
                minLength = length;
                result = parameters;
            }
        }
    return result;
}

void DubinsPath::potentiallyReplaceParameters(const DubinsPath::DubinsParameters& parameters) {
    if (m_DubinsParameters.pathType == INVALID || m_DubinsParameters.totalLength() > parameters.totalLength()) {
        m_DubinsParameters = parameters;
    }
}

void DubinsPath::sample(double& x, double& y, double& theta, double& distance) const {
    auto t = distance / m_Radius;

}

DubinsPath::SegmentType DubinsPath::getSegmentType(int index, PathType pathType) {
    switch (pathType) {
        case LRL:
            switch (index) {
                case 0: return L;
                case 1: return R;
                case 2: return L;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case LSL:
            switch (index) {
                case 0: return L;
                case 1: return S;
                case 2: return L;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case LSR:
            switch (index) {
                case 0: return L;
                case 1: return S;
                case 2: return R;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case RLR:
            switch (index) {
                case 0: return R;
                case 1: return L;
                case 2: return R;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case RSL:
            switch (index) {
                case 0: return R;
                case 1: return S;
                case 2: return L;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case RSR:
            switch (index) {
                case 0: return R;
                case 1: return S;
                case 2: return R;
                default: throw std::runtime_error("Unknown Dubins segment");
            }
        case INVALID:
        default: throw std::runtime_error("Unknown Dubins segment");
    }
}
