#include "getrgeomath.h"

using VEbenenLagen = std::vector<SEbene>;
using VPolDreieck  = std::vector<SPoint>;
using VGelenke     = std::vector<SPoint>;
using A3Gelenke    = std::array<SPoint, 3>;


SEbene FixedLenLine(SLine & roL, double const & crnLen, bool const & crbFirst)
    {
    auto const dx   { roL.x1 - roL.x2 };
    auto const dy   { roL.y1 - roL.y2 };
    auto const nLen { sqrt(dx*dx + dy*dy) };
    auto const q    { (double)crnLen / ((nLen!=0)?nLen:1) };
    if (crbFirst)
        {
        roL.x2 = roL.x1 - dx*q;
        roL.y2 = roL.y1 - dy*q;
        }
    else
        {
        roL.x1 = roL.x2 + dx*q;
        roL.y1 = roL.y2 + dy*q;
        }
    return roL;
    } // void FixedLenLine(...



SPoint Intersection(SLine const & E1, SLine const & E2)
    {
    auto const dx1 { E1.x2 - E1.x1 };
    auto const dx2 { E2.x2 - E2.x1 };

    auto const m1 { (E1.y2 - E1.y1) / dx1 }; // Steigungen ermitteln
    auto const m2 { (E2.y2 - E2.y1) / dx2 };

    // if (ROUND(m1,MAX_ACCURACY)==ROUND(m2,MAX_ACCURACY)) return false; // Geraden sind parallel

    auto const n1 { E1.y1 - (m1*E1.x1) }; // Abstände von X-Achse ermitteln
    auto const n2 { E2.y1 - (m2*E2.x1) };

    auto const  x { (n2-n1)/(m1-m2) }; // Schnittpunktkoordinate berechnen
    auto const  y { m1*x+n1 };

    return { x, y };
  } // Intersection


SEbene Perpendicle(SLine const & croLine)
    {
    auto const dx = (croLine.x2 - croLine.x1)/2.0;
    auto const dy = (croLine.y2 - croLine.y1)/2.0;

    SEbene I{ croLine.x2 - dy - dx,
	      croLine.y2 + dx - dy,
	      croLine.x2 + dy - dx,
	      croLine.y2 - dx - dy};

    return std::move(I);
    }

SPoint PointMirror(SPoint const & croPoint, SLine const & croMirror)
    {
    auto const dx = (croMirror.x2 - croMirror.x1)/2.0;
    auto const dy = (croMirror.y2 - croMirror.y1)/2.0;

    SEbene I{ croPoint.x - dy,
	      croPoint.y + dx,
	      croPoint.x + dy,
	      croPoint.y - dx };

    auto const S  { Intersection(croMirror, I) };
    auto const mx { (croMirror.x2 + croMirror.x1)/2.0 };
    auto const my { (croMirror.y2 + croMirror.y1)/2.0 };

    return { S.x + (S.x - croPoint.x), S.y + (S.y - croPoint.y) };
    }

SPoint CalcPolpunkt(SEbene const & E1, SEbene const & E2)
    {
    SEbene const L1{ E1.x1, E1.y1, E2.x1, E2.y1 };
    auto   const La{ Perpendicle(L1) };

    SEbene const L2{ E1.x2, E1.y2, E2.x2, E2.y2 };
    auto   const Lb{ Perpendicle(L2) };

    return std::move(Intersection(La, Lb));
    }

SUmkreis Umkreis( SPoint const & P1, SPoint const & P2, SPoint const & P3 )
    {
    SEbene const L1{ P1.x, P1.y, P2.x, P2.y };
    auto   const La{ Perpendicle(L1) };

    SEbene const L2{ P1.x, P1.y, P3.x, P3.y };
    auto   const Lb{ Perpendicle(L2) };

    auto const M{Intersection(La, Lb)};
    auto const R{sqrt( (double)(M.x - L1.x1)*(M.x - L1.x1) + (M.y - L1.y1)*(M.y - L1.y1) )};

    return { R, M };
    }
