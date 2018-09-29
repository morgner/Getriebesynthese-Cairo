#include "getrgeomath.h"

using VEbenenLagen = std::vector<SEbene>;
using VPolDreieck  = std::vector<SPoint>;
using VGelenke     = std::vector<SPoint>;
using A3Gelenke    = std::array<SPoint, 3>;



SEbene FixedLenLine(SEbene & roL, double const & crnLenEbene, bool const & crbFirst)
    {
    auto const dx   { roL.x1 - roL.x2 };
    auto const dy   { roL.y1 - roL.y2 };
    auto const nLen { sqrt(dx*dx + dy*dy) };
    auto const q    { (double)crnLenEbene / ((nLen!=0)?nLen:1) };
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


