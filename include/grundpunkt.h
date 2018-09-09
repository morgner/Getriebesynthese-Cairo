#ifndef __GRUNDPUNKT_H
#define __GRUNDPUNKT_H

/*
#include <fstream>
#include <iostream>     // std::cout
#include <algorithm>
#include <regex>
#include <map>
#include <utility>      // std::pair

#include <string>       // std::string
#include <streambuf>
#include <sstream>      // std::ostringstream
*/

#include <vector>
#include <array>
#include "te.h"

#include <gtkmm/drawingarea.h> // Cairo

using namespace std::string_literals; 
using TRenderer=Cairo::RefPtr<Cairo::Context>;


struct SPoint
    {
    SPoint() = default;
    template<typename T>
	SPoint(T x, T y) : x(x), y(y) {}

    double x{0}, y{0};
    };

struct SEbene
    {
    SPoint M() const { return { (x1+x2)/2,(y1+y2)/2 }; }
    double x1{0}, y1{0};
    double x2{0}, y2{0};
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

SEbene FixedLenLine(SEbene & roL, double const & crnLenEbene, bool const & crbFirst = true)
    {
    double const dx   = roL.x1 - roL.x2;
    double const dy   = roL.y1 - roL.y2;
    double const nLen = sqrt(dx*dx + dy*dy);
    double const q    = (double)crnLenEbene / ((nLen!=0)?nLen:1);
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
    double dx1,dx2,m1,n1,m2,n2;

    dx1 = E1.x2 - E1.x1;
    dx2 = E2.x2 - E2.x1;

    m1 = (E1.y2 - E1.y1) / dx1; // Steigungen ermitteln
    m2 = (E2.y2 - E2.y1) / dx2;

    // if (ROUND(m1,MAX_ACCURACY)==ROUND(m2,MAX_ACCURACY)) return false; // Geraden sind parallel

    n1 = E1.y1 - (m1*E1.x1); // Abstände von X-Achse ermitteln
    n2 = E2.y1 - (m2*E2.x1);

    double x = (n2-n1)/(m1-m2); // Schnittpunktkoordinate berechnen
    double y = m1*x+n1;

    return { x, y };
  } // Intersection


SEbene Perpendicle(SEbene const & croLine)
    {
    auto const dx = (croLine.x2 - croLine.x1)/2.0;
    auto const dy = (croLine.y2 - croLine.y1)/2.0;

    SEbene I;
    I.x1 = croLine.x2 - dy - dx;
    I.y1 = croLine.y2 + dx - dy;
    I.x2 = croLine.x2 + dy - dx;
    I.y2 = croLine.y2 - dx - dy;

    return std::move(I);
    }

SPoint PointMirror(TRenderer const & cr, SPoint const & croPoint, SEbene const & croMirror)
    {
    auto const dx = (croMirror.x2 - croMirror.x1)/2.0;
    auto const dy = (croMirror.y2 - croMirror.y1)/2.0;

    SEbene I;
    I.x1 = croPoint.x - dy - 0*dx;
    I.y1 = croPoint.y + dx - 0*dy;
    I.x2 = croPoint.x + dy - 0*dx;
    I.y2 = croPoint.y - dx - 0*dy;

    auto const S  = Intersection(croMirror, I);
    auto const mx = (croMirror.x2 + croMirror.x1)/2.0;
    auto const my = (croMirror.y2 + croMirror.y1)/2.0;

    return { S.x + (S.x - croPoint.x), S.y + (S.y - croPoint.y) };
    }

SPoint CalcPolpunkt(SEbene const & E1, SEbene const & E2)
    {
    SEbene const L1{ E1.x1, E1.y1, E2.x1, E2.y1 };
    auto La = Perpendicle(L1);

    SEbene const L2{ E1.x2, E1.y2, E2.x2, E2.y2 };
    auto Lb = Perpendicle(L2);

    return Intersection( La, Lb );
    }

SUmkreis RenderUmkreis( TRenderer const & cr, SPoint const & P1, SPoint const & P2, SPoint const & P3 )
    {
    SEbene const L1{ (double)P1.x, (double)P1.y, (double)P2.x, (double)P2.y };
    auto La = Perpendicle(L1);
    SEbene const L2{ (double)P1.x, (double)P1.y, (double)P3.x, (double)P3.y };
    auto Lb = Perpendicle(L2);

    auto M = Intersection( La, Lb );
    auto R = sqrt( (double)(M.x - L1.x1)*(M.x - L1.x1) + (M.y - L1.y1)*(M.y - L1.y1) );

    cr->set_source_rgba(1, 0, 0,0.5);
    cr->arc(M.x, M.y, R, 0, 2*M_PI);
    cr->fill();
    cr->stroke();

    return { R, M };
    }

SUmkreis Umkreis( SPoint const & P1, SPoint const & P2, SPoint const & P3 )
    {
    SEbene const L1{ (double)P1.x, (double)P1.y, (double)P2.x, (double)P2.y };
    auto La = Perpendicle(L1);
    SEbene const L2{ (double)P1.x, (double)P1.y, (double)P3.x, (double)P3.y };
    auto Lb = Perpendicle(L2);

    auto M = Intersection( La, Lb );
    auto R = sqrt( (double)(M.x - L1.x1)*(M.x - L1.x1) + (M.y - L1.y1)*(M.y - L1.y1) );

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
	SUmkreis m_tUK{};
	A3Gelenke m_a3Gelenke{}; // Gelenkpunkte
	bool      m_bFixed{true};


    public:

	CGrundpunkt(SPoint const & P123) : m_dX(P123.x), m_dY(P123.y) {}
	CGrundpunkt(CGrundpunkt const & src) = default;

	void FixIt( bool bFixit = true ) { m_bFixed = bFixit; }
	bool Isfix() { return m_bFixed; }

	SPoint P123() const { return {m_dX, m_dY}; }
	SPoint const & G0() const { return m_tUK.MidPnt; }
	SPoint const & GPoint( int i ) const { return m_a3Gelenke[i]; }

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
	    cr->set_line_width(2); //!!

	    // P123
	    cr->set_source_rgba(0, .5, .5, 0.75);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->fill();
	    cr->set_source_rgb (0, 0, 0);
	    cr->arc(m_dX, m_dY, 8, 0, 2*M_PI);
	    cr->stroke();

	    A3Gelenke G;
	    for ( int n=0, i=0, j=0; n < 3; ++n )
		{
		// Lines: P12-P13, P12-P23, P13-P23
		if (n == 0) { i=0; j=1; }
		if (n == 1) { i=0; j=2; }
		if (n == 2) { i=1; j=2; }
		SEbene PL{ (double)Poldreieck[i].x, (double)Poldreieck[i].y,
			    (double)Poldreieck[j].x, (double)Poldreieck[j].y};
		G[n] = PointMirror( cr, { m_dX, m_dY }, PL );
		cr->set_source_rgba(1, 0, 0, 0.25);
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->fill();
		cr->arc(m_dX, m_dY, 20, 0, 2*M_PI);
		cr->set_source_rgb(0, 0, 0);
		cr->stroke();
		} // for (...
	    m_a3Gelenke = G;

	    m_tUK = Umkreis( G[0], G[1], G[2] );
	    cr->set_source_rgba(.5, .5, 0, 0.5);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->fill();
	    cr->set_line_width(2);
	    cr->set_source_rgb(0, 0, 0);
	    cr->arc(m_tUK.MidPnt.x, m_tUK.MidPnt.y, 36, 0, 2*M_PI);
	    cr->stroke();
	    }
    }; // class CGrundpunkt

using VGrundpunkte  = std::vector<CGrundpunkt>;
using A2Grundpunkte = std::array<CGrundpunkt, 2>;

// __GRUNDPUNKT_H
#endif
