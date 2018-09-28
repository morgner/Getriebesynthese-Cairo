#ifndef __GRUNDPUNKT_H
#define __GRUNDPUNKT_H


#include <vector>
#include <array>
#include <math.h>
#include "te.h"

using namespace std::string_literals;

using TRenderer = Cairo::RefPtr<Cairo::Context>;


struct SPoint
    {
    SPoint() = default;
    SPoint(double const & x, double const & y) : x(x), y(y) {}
    double x{0}, y{0};
    };

struct SEbene
    {
    SPoint M() const { return { (x1+x2)/2,(y1+y2)/2 }; }
    double x1{0}, y1{0};
    double x2{0}, y2{0};
    };

using SLine = SEbene;

struct SUmkreis
    {
    double Radius;
    SPoint MidPnt;
    };

using VEbenenLagen = std::vector<SEbene>;
using VPolDreieck  = std::vector<SPoint>;
using VGelenke     = std::vector<SPoint>;
using A3Gelenke    = std::array<SPoint, 3>;

SEbene FixedLenLine(SEbene & roL, double const & crnLenEbene, bool const & crbFirst = true)
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

SPoint Intersection(SEbene const & E1, SEbene const & E2)
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


SEbene Perpendicle(SEbene const & croLine)
    {
    auto const dx = (croLine.x2 - croLine.x1)/2.0;
    auto const dy = (croLine.y2 - croLine.y1)/2.0;

    SEbene I{ croLine.x2 - dy - dx,
	      croLine.y2 + dx - dy,
	      croLine.x2 + dy - dx,
	      croLine.y2 - dx - dy};

    return std::move(I);
    }

SPoint PointMirror(TRenderer const & cr, SPoint const & croPoint, SEbene const & croMirror)
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

/*
 * CGrundpunkt
 *
 *
 */
class CGrundpunkt
    {
    protected:

	double    m_dX{};
	double    m_dY{};
	SUmkreis  m_tUK{};
	A3Gelenke m_a3Gelenke{}; // Gelenkpunkte
	bool      m_bFixed{true};


    public:

	CGrundpunkt(SPoint const & P123) : m_dX(P123.x), m_dY(P123.y) {}
	CGrundpunkt(CGrundpunkt const & src) = default;

	void FixIt( bool const bFixit = true ) { m_bFixed = bFixit; }
	bool Isfix() const { return m_bFixed; }

	SPoint const   P123() const { return {m_dX, m_dY}; }
	SPoint const & G0() const { return m_tUK.MidPnt; }
	SPoint const & GPoint( int const & i ) const { return m_a3Gelenke[i]; }

	void UpdateAndShow( TRenderer const & cr, SPoint const & P123, VPolDreieck const & Poldreieck)
	    {
	    Update(P123);
	    Show(cr, Poldreieck);
	    }

	void Update(SPoint const & P123)
	    {
	    m_dX = P123.x;
	    m_dY = P123.y;
	    }

	void Show( TRenderer const & cr, VPolDreieck const & Poldreieck)
	    {
	    cr->set_line_width(1); //!!

	    // P123
	    cr->set_source_rgba(0, .5, .5, 0.75);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->fill();
	    cr->set_source_rgb (0, 0, 0);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->stroke();

	    A3Gelenke & G = m_a3Gelenke;
	    for ( int n=0, i=0, j=0; n < 3; ++n )
		{
		// Lines: P12-P13, P12-P23, P13-P23
		switch (n)
		    {
		    case 0: i=0; j=1; break;
		    case 1: i=0; j=2; break;
		    case 2: i=1; j=2; break;
		    }
		SLine const PL{ (double)Poldreieck[i].x, (double)Poldreieck[i].y,
			        (double)Poldreieck[j].x, (double)Poldreieck[j].y};
		G[n] = PointMirror( cr, { m_dX, m_dY }, PL );
		cr->set_source_rgba(1, 0, 0, 0.25);
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->fill();
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->set_source_rgb(0, 0, 0);
		cr->stroke();
		} // for (...

	    m_tUK = Umkreis( G[0], G[1], G[2] );
/*
	    cr->set_source_rgba(.5, .5, 0, 0.5);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->fill();
	    cr->set_line_width(1);
	    cr->set_source_rgb(0, 0, 0);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->stroke();
*/
	    }
    }; // class CGrundpunkt

using VGrundpunkte  = std::vector<CGrundpunkt>;
using A2Grundpunkte = std::array<CGrundpunkt, 2>;

// __GRUNDPUNKT_H
#endif
