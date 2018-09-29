#ifndef __GETRGEOMATH_H
#define __GETRGEOMATH_H

#include <array>
#include <vector>
#include <math.h>

struct SPoint
    {
    SPoint() = default;
    SPoint(double const & x, double const & y) : x(x), y(y) {}
    template<typename T>
	SPoint(T const & t) : x(t.x), y(t.y) {}
    double x{0}, y{0};

    SPoint operator -= (SPoint const & p)
        {
	x-=p.x;
	y-=p.y;
        return *this;
        }
    };

inline SPoint operator - (SPoint const & p1, SPoint const & p2)
    {
    return {p1.x-p2.x, p1.y-p2.y};
    }

inline SPoint operator + (SPoint const & p1, SPoint const & p2)
    {
    return {p2.x+p1.x, p2.y+p1.y};
    }

inline SPoint operator / (SPoint const & p, double const & d)
    {
    return {p.x/d, p.y/d};
    }

inline SPoint operator * (SPoint const & p, double const & d)
    {
    return {p.x*d, p.y*d};
    }


struct SPointT
    {
    SPointT() = default;
    SPointT(double const & x, double const & y, bool bSplit) : x(x), y(y), bCSplit(bSplit) {}
    template<typename T>
	SPointT(T const & t, bool bSplit) : x(t.x), y(t.y), bCSplit(bSplit) {}
    double x{0}, y{0};
    bool   bCSplit{false};
    };

struct SEbene
    {
    SEbene() = default;
//    SEbene(SPoint const & p1, SPoint const & p2) : x1(p1.x), y1(p1.y), x2(p2.x), y2(p2.y) {}
    SPoint M() const { return { (x1+x2)/2,(y1+y2)/2 }; }
    double x1{0}, y1{0};
    double x2{0}, y2{0};
    };

using SLine = SEbene;

struct SCollision
    {
    enum class EWhat
	{
	none,		// there was no collision
	A,		// move point A
	Ebene,		// move a Ebenenlage
	Grundpunkt,	// move a Grundpunkt
	}eWhat {EWhat::none};
    int nIndex {0};	// 1,2,3 = E1, E2, E3 ; 1,2 =G1, G2
    int nSubIx {0};	// 1,2,3 = P1, P2, Pm
    };

struct SUmkreis
    {
    double Radius;
    SPoint MidPnt;
    };

using VEbenenLagen = std::vector<SEbene>;
using VPolDreieck  = std::vector<SPoint>;
using VGelenke     = std::vector<SPoint>;
using A3Gelenke    = std::array<SPoint, 3>;


SEbene FixedLenLine(SLine & roL, double const & crnLen, bool const & crbFirst = true);
SPoint Intersection(SLine const & E1, SLine const & E2);
SEbene Perpendicle(SLine const & croLine);
SPoint PointMirror(SPoint const & croPoint, SLine const & croMirror);
SPoint CalcPolpunkt(SEbene const & E1, SEbene const & E2);
SUmkreis Umkreis( SPoint const & P1, SPoint const & P2, SPoint const & P3 );

// __GETRGEOMATH_H
#endif
