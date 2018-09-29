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
    };

struct SEbene
    {
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


SEbene FixedLenLine(SEbene & roL, double const & crnLenEbene, bool const & crbFirst = true);


// __GETRGEOMATH_H
#endif
